import os
import numpy as np
import subprocess
from collections import OrderedDict, Counter, defaultdict
from cpython cimport PyBytes_FromStringAndSize

######################
### Local Assembly ###
######################
cdef int getMinEdge(int numSegRaw):
    if numSegRaw < 5:
        return 1
    elif numSegRaw < 10:
        return 2
    elif numSegRaw < 40:
        return 3
    elif numSegRaw < 70:
        return 4
    else:
        return 5


cpdef assembleCluster(Cluster[::1] cltView, dict allCltData, int startIdx, int taskSize, object cmdArgs):
    cdef int numThread = cmdArgs.numThread
    cdef int nodeLen = cmdArgs.nodeLen
    cdef int i, endIdx, minEdge, rounds
    cdef str cmd, prefix

    endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]

    for i in range(startIdx, endIdx):
        if isSomaClt(&cltView[i]):
            continue

        minEdge = getMinEdge(cltView[i].numSegRaw)
        if cmdArgs.minEdge > 0:
            minEdge = cmdArgs.minEdge

        # 1. Primary assembling
        prefix = "tmp_assm/{}_{}".format(cltView[i].tid, cltView[i].idx)
        rounds = 0
        while rounds < 5:
            rounds += 1
            cmd = "wtdbg2 -p 5 -k 15 -l 256 -e {0} -S 1 -A --rescue-low-cov-edges --node-len {1} --ctg-min-length {1} " \
                "--ctg-min-nodes 1 -q -t {2} -i {3}.fa -fo {3}".format(minEdge, nodeLen, numThread, prefix)
            subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
            if os.path.isfile(f"{prefix}.ctg.lay.gz") == False:
                continue
            
            cmd = "wtpoa-cns -q -c 1 -t {0} -i {1}.ctg.lay.gz -fo {1}_assm.fa".format(numThread, prefix)
            subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
            if os.path.isfile(f"{prefix}_assm.fa") and (os.path.getsize(f"{prefix}_assm.fa") != 0):
                break
        
        if (not os.path.isfile(f"{prefix}_assm.fa")) or (os.path.getsize(f"{prefix}_assm.fa") == 0):
            if not outputReadAsAssm(cltView, allCltData, cmdArgs, i):
                continue

        # 2. First round polishing
        cmd = "minimap2 -aY {0}_assm.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
            "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 3. Second round polishing
        cmd = "minimap2 -aY {0}_assembled.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
            "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # 4. Recalibration
        recalibration(prefix, cmdArgs)

        # 5. Remove N-bases
        cmd = (
            "sed 's/N//g' {0}_assembled.fa | "
            "awk '{{if($1~/^>/){{ctg=$1; a[ctg]=\"\"}} else{{a[ctg]=a[ctg]\"\"$1}}}} "
            "END{{for(ctg in a){{if(a[ctg]!=\"\"){{print ctg; print a[ctg]}}}}}}' > {0}_tmp.fa && "
            "mv {0}_tmp.fa {0}_assembled.fa"
        ).format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


#####################
### Recalibration ###
#####################
cdef recalibration(str prefix, object cmdArgs):
    # 1. Map raw reads to polished seqs
    os.rename(f"{prefix}_assembled.fa", f"{prefix}_polished.fa")
    cmd = "minimap2 -aY {0}_polished.fa {0}.fa | " \
        "samtools sort | samtools view -bhS -F 3332 -o {0}_RawToPolish.bam && " \
        "samtools index {0}_RawToPolish.bam".format(prefix)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
    if not os.path.exists(f"{prefix}_RawToPolish.bam.bai"):
        os.rename(f"{prefix}_polished.fa", f"{prefix}_assembled.fa")
        return
    
    # 2. Load files
    cdef bytes inputFn = f"{prefix}_polished.fa".encode("utf-8")
    cdef faidx_t *inputFa = fai_load(inputFn)
    cdef BamFile inputBam = BamFile(f"{prefix}_RawToPolish.bam", "rb")
    cdef object outputFa = open(f"{prefix}_assembled.fa", "w")

    # 3. Initilize
    cdef int tid, numRef = faidx_nseq(inputFa)
    cdef int refLen, maxLen = 0
    for tid in range(numRef):
        refLen = faidx_seq_len(inputFa, faidx_iseq(inputFa, tid))
        if refLen > maxLen:
            maxLen = refLen

    # 4. Recalibrate
    cdef object queryArr = np.zeros(5*maxLen, dtype=np.int32)
    cdef object refArr = np.zeros(5*maxLen, dtype=np.int32)
    cdef int refStart, refEnd, skipNext
    cdef const char *cRefName
    cdef char *cRefSeq
    cdef str refSeq, recalibratedSeq

    for tid in range(numRef):
        # 5. Load sequence
        cRefName = faidx_iseq(inputFa, tid)
        refLen = faidx_seq_len(inputFa, cRefName)
        cRefSeq = faidx_fetch_seq(inputFa, cRefName, 0, refLen, &refLen)
        refSeq = cRefSeq.decode('utf-8')
        
        # 6. Define homopolymer regions
        polymerRegions = getPolymerRegions(refSeq, refLen)
        
        # 7. Collect query homo-polymers
        iterator = Iterator(inputBam, tid)
        queryPolymers = getQueryPolymers(iterator, queryArr, refArr, polymerRegions, refSeq)
        
        # 8. Output recalibrated sequence
        outputFa.write(">" + cRefName.decode('utf-8') + "\n")
        skipNext = 0
        for refStart in polymerRegions:
            refEnd = polymerRegions[refStart]

            # Ignore long homopolymer, un-covered region
            if (refEnd - refStart + 1 > 20) or (refStart not in queryPolymers):
                outputFa.write(refSeq[refStart:refEnd+1])
                continue
            
            # If previous base/homo-polymer has been recalibrated, skip this one
            if skipNext:
                outputFa.write(refSeq[refStart:refEnd+1])
                skipNext = 0
                continue
            
            recalibratedSeq, skipNext = getRecalibratedSeq(queryPolymers, refSeq, refStart, refEnd)
            outputFa.write(recalibratedSeq)
        
        outputFa.write("\n")
        del iterator

    # 9. Close files
    fai_destroy(inputFa)
    outputFa.close()
    inputBam.close()


cdef object getPolymerRegions(str refSeq, int refLen):
    cdef int start = 0
    cdef int end = 1
    cdef object regions = OrderedDict()

    while end < refLen:
        if refSeq[end] == refSeq[start]:
            end += 1
            continue
        
        regions[start] = end - 1
        start = end
        end += 1
    
    # Final homo-polymer
    regions[start] = end - 1
    return regions


cdef object getQueryPolymers(Iterator iterator, int[::1] queryArr, int[::1] refArr, object polymerRegions, str refSeq):
    cdef int refPos, refStart, refEnd
    cdef int queryPos, queryStart, queryEnd
    cdef int i, retValue, numPairs
    cdef int refLen = len(refSeq)
    cdef object queryPolymers = OrderedDict()
    cdef str readSeq, querySeq
    
    while True:
        # 1. Load read
        retValue = iterator.cnext1()
        if retValue < 0:
            return queryPolymers

        # 2. Get aligned-pairs
        readSeq = getReadSeq(iterator)
        numPairs = getAlignedPairs(iterator.bamRcord, &queryArr[0], &refArr[0])
        refEnd = -2

        # 3. Collect homo-polymer lengths
        for i in range(numPairs):
            queryPos = queryArr[i]
            refPos = refArr[i]

            # 4. Read cover homo-polymer start
            if refPos in polymerRegions:
                queryStart = queryPos
                refStart = refPos
                refEnd = polymerRegions[refPos]

                # Single base
                if refStart == refEnd:
                    if queryStart < 0:
                        queryStart = -1
                        queryEnd = -2
                        querySeq = ""
                    else:
                        # Find real start
                        queryPos = checkLeftSide(queryArr, refArr, i-1, 0)
                        if queryPos >= 0:
                            queryStart = queryPos + 1 
                        # Find real end
                        queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, refLen)
                        if queryPos >= 0:
                            queryEnd = queryPos - 1

                        querySeq = readSeq[queryStart:queryEnd+1]
                    
                    if queryStart == queryEnd:
                        extraLen = getExtraLen(queryStart, readSeq, refSeq[refStart])
                    else:
                        extraLen = 0
                    
                    if refStart in queryPolymers:
                        queryPolymers[refStart].append((querySeq, queryEnd-queryStart+1, extraLen))
                    else:
                        queryPolymers[refStart] = [(querySeq, queryEnd-queryStart+1, extraLen)]
                    continue
                
                # Homo-polymer region
                if queryStart < 0:
                    # Correction for boundary deletion
                    queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, refEnd)
                    if queryPos < 0:
                        refEnd = -2
                    else:
                        queryStart = queryPos
                else:
                    # Correction for boundary insertion
                    queryPos = checkLeftSide(queryArr, refArr, i-1, 0)
                    if queryPos >= 0:
                        queryStart = queryPos + 1            
            
            # 5. Read cover homo-polymer end
            elif refPos == refEnd:
                queryEnd = queryPos
                if queryEnd < 0:
                    # Correction for boundary deletion
                    queryPos = checkLeftSide(queryArr, refArr, i-1, refStart)
                    if queryPos < 0:
                        refEnd = -2
                        continue
                    else:
                        queryEnd = queryPos
                else:
                    queryPos = checkRightSide(queryArr, refArr, numPairs, i+1, refLen)
                    if queryPos >= 0:
                        queryEnd = queryPos - 1
                
                if queryEnd - queryStart == refEnd - refStart:
                    extraLen = getExtraLen(queryStart, readSeq, refSeq[refStart])
                else:
                    extraLen = 0
                
                # 6. Collect query homo-polymers
                if refStart in queryPolymers:
                    queryPolymers[refStart].append((readSeq[queryStart:queryEnd+1], queryEnd - queryStart + 1, extraLen))
                else:
                    queryPolymers[refStart] = [(readSeq[queryStart:queryEnd+1], queryEnd - queryStart + 1, extraLen)]
                
                refEnd = -2


cdef str getReadSeq(Iterator iterator):
    cdef int i, readLen = iterator.bamRcord.core.l_qseq
    cdef bytes bSeq = PyBytes_FromStringAndSize(NULL, readLen)
    cdef char *cSeq = <char*>bSeq

    for i in range(readLen):
        cSeq[i] = seq_nt16_str[bam_seqi(bam_get_seq(iterator.bamRcord), i)]

    return bSeq.decode('utf-8')


cdef int checkLeftSide(int[::1] queryArr, int[::1] refArr,  int i, int leftMost):
    while i >= 0:
        queryPos = queryArr[i]
        refPos = refArr[i]
        if (queryPos < 0) or (refPos < 0):
            i -= 1
            continue
        if refPos < leftMost:
            return -1
        return queryPos

    return -1


cdef int checkRightSide(int[::1] queryArr, int[::1] refArr, int numPairs,  int i, int rightMost):
    while i < numPairs:
        queryPos = queryArr[i]
        refPos = refArr[i]
        if (queryPos < 0) or (refPos < 0):
            i += 1
            continue
        if refPos > rightMost:
            return -1
        return queryPos

    return -1


cdef int getExtraLen(int queryStart, str readSeq, str polymerBase):
    if queryStart == 0:
        return 0    
    if readSeq[queryStart-1] == polymerBase:
        return 1
    else:
        return 0


cdef int getMostCommonLen(object queryPolymers, int refStart, int refEnd):
    cdef object counter = Counter([p[1] + p[2] for p in queryPolymers[refStart]])
    cdef int refLen = refEnd - refStart + 1
    cdef int refCount = counter[refLen]
    cdef int altLen, altCount

    for altLen, altCount in counter.most_common(2):
        if altLen != refLen:
            break
    
    if altCount <= refCount:
        return refLen
    else:
        return altLen


cdef tuple getRecalibratedSeq(object queryPolymers, str refSeq, int refStart, int refEnd):
    # 1. Get most common length
    cdef int mostCommonLen = getMostCommonLen(queryPolymers, refStart, refEnd)
    if mostCommonLen == 0:
        return '', 1

    # 2. Initialize counter for each position
    cdef list positionCounters = [defaultdict(int) for _ in range(mostCommonLen)]
    cdef object polymer
    cdef int i, queryLen
    cdef str querySeq, base
    
    # 3. Count for each position
    for polymer in queryPolymers[refStart]:
        queryLen = polymer[1]
        if queryLen != mostCommonLen:
            continue
        
        querySeq = polymer[0]
        for i in range(mostCommonLen):
            base = querySeq[i]
            positionCounters[i][base] += 1
    
    # 4. Get most common base for each position
    cdef list result = []
    cdef object counter
    cdef str refBase = refSeq[refStart]
    cdef int refCount
    if len(positionCounters[0]) == 0:
        return refSeq[refStart:refEnd+1], 0
    else:
        for counter in positionCounters:
            base = max(counter, key=counter.get)
            refCount = counter[refBase] # 0 if refBase not exists
            if refCount >= counter[base]:
                result.append(refBase)
            else:
                result.append(base)

        if (refEnd - refStart + 1) == mostCommonLen:
            return ''.join(result), 0
        else:
            return ''.join(result), 1