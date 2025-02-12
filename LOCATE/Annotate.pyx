import os
import numpy as np
import subprocess

##################
### Data Types ###
##################
AnnoDt = np.dtype([
    ('idx',         np.int32),
    ('cltTid',      np.int32),
    ('cltIdx',      np.int32),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
    ('strand',      np.uint8),
    ('tid',         np.int32),
    ('refStart',    np.int32),
    ('refEnd',      np.int32),
    ('flag',        np.uint32),
])


#########################
### Annotate Assembly ###
#########################
cdef annotateAssm(Cluster[::1] cltView, int startIdx, int endIdx, object cmdArgs):
    cdef int i
    for i in range(startIdx, endIdx):
        mapFlankToAssm(cltView[i].tid, cltView[i].idx)
        extractIns(&cltView[i])
        mapAssmFlankToLocal(cltView[i].tid, cltView[i].idx)
        if os.path.isfile("tmp_anno/{}_{}_AssmFlankToLocal.bam".format(cltView[i].tid, cltView[i].idx)) == True:
            reExtractIns(&cltView[i])
        mapInsToTE(cltView[i].tid, cltView[i].idx, cmdArgs)
        

cdef mapFlankToAssm(int tid, int idx):
    cdef str targetFn = "tmp_assm/{}_{}_assembled.fa".format(tid, idx)
    cdef str queryFn = "tmp_anno/{}_{}_flank.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_FlankToAssm.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(targetFn, queryFn, outFn)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef mapAssmFlankToLocal(int tid, int idx):
    cdef str targetFn = "tmp_anno/{}_{}_local.fa".format(tid, idx)
    cdef str queryFn = "tmp_anno/{}_{}_assmFlank.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_AssmFlankToLocal.bam".format(tid, idx)
    cdef str cmd = "minimap2 -x sr --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(targetFn, queryFn, outFn)
    if os.path.isfile(queryFn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef mapInsToTE(int tid, int idx, object cmdArgs):
    cdef str queryFn = "tmp_anno/{}_{}_insertion.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_InsToTE.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(cmdArgs.teFn, queryFn, outFn)
    if os.path.isfile(queryFn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


##########################
### Annotate Insertion ###
##########################
cdef object getClassArr(object cmdArgs):
    cdef int numTE = 0, maxNum = 190
    cdef object classArr = np.zeros(200, dtype=np.uint32)
    cdef uint32_t[::1] classView = classArr

    for l in open(cmdArgs.classFn, "r"):
        l = l.strip().split()
        if len(l) < 2:
            continue
        
        if numTE >= maxNum:
            maxNum = classView.shape[0] + 200
            classArr.resize((maxNum,), refcheck=False)
            classView = classArr
            maxNum -= 10
        
        if l[1] == "DNA":
            classView[numTE] = CLT_DNA
        elif l[1] == "LTR":
            classView[numTE] = CLT_LTR
        elif l[1] == "LINE":
            classView[numTE] = CLT_LINE
        elif l[1] == "SINE":
            classView[numTE] = CLT_SINE
        elif l[1] == "RETROPOSON":
            classView[numTE] = CLT_RETROPOSON
        
        numTE += 1
    
    classArr.resize((numTE,), refcheck=False)
    return classArr


cdef object getSizeArr(object cmdArgs):
    cdef int numTE = 0, maxNum = 190
    cdef object sizeArr = np.zeros(200, dtype=np.int32)
    cdef int32_t[::1] sizeView = sizeArr
    cdef indexFn = cmdArgs.teFn + ".fai"

    for l in open(indexFn, "r"):
        l = l.strip().split()
        if len(l) < 1:
            continue
        
        if numTE >= maxNum:
            maxNum = sizeView.shape[0] + 200
            sizeArr.resize((maxNum,), refcheck=False)
            sizeView = sizeArr
            maxNum -= 10
        
        sizeView[numTE] = int(l[1])
        numTE += 1
    
    sizeArr.resize((numTE,), refcheck=False)
    return sizeArr


cdef object getLtrArr():
    cdef int numTE = 0, maxNum = 190
    cdef object ltrArr = np.zeros(200, dtype=np.int32)
    cdef int32_t[::1] ltrView = ltrArr

    for l in open("tmp_anno/ltrSize.txt", "r"):
        l = l.strip().split()
        if len(l) < 1:
            continue
        
        if numTE >= maxNum:
            maxNum = ltrView.shape[0] + 200
            ltrArr.resize((maxNum,), refcheck=False)
            ltrView = ltrArr
            maxNum -= 10
        
        ltrView[numTE] = int(l[0])
        numTE += 1
    
    ltrArr.resize((numTE,), refcheck=False)
    return ltrArr


cdef object annotateIns(Cluster[::1] cltView, int startIdx, int endIdx, object cmdArgs):
    cdef int i, numTmp, numAnno = 0, maxNum = 2900
    cdef object annoArr = np.zeros(3000, dtype=AnnoDt)
    cdef object classArr = getClassArr(cmdArgs)
    cdef object sizeArr = getSizeArr(cmdArgs)
    cdef object ltrArr = getLtrArr()
    cdef Annotation[::1] annoView = annoArr
    cdef uint32_t[::1] classView = classArr
    cdef int[::1] sizeView = sizeArr
    cdef int[::1] ltrView = ltrArr
    cdef str bamFn, assmFn

    for i in range(startIdx, endIdx):
        # 1. Check if cluster has assembled insSeq
        assmFn = "tmp_assm/{}_{}_polished.fa".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(assmFn) != 0:
                cltView[i].flag |= CLT_ASSEMBLED
        
        bamFn = "tmp_anno/{}_{}_InsToTE.bam".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(bamFn) == False:
            continue

        if numAnno >= maxNum:
            maxNum = annoView.shape[0] + 3000
            annoArr.resize((maxNum,), refcheck=False)
            annoView = annoArr
            maxNum -= 100

        # 2. Annotate TE fragment, polyA/T for insSeq
        numTmp = fillAnnoArr(&cltView[i], &annoView[numAnno], &classView[0], i)

        # 3. Annotate TSD
        mapTsdToLocal(cltView[i].tid, cltView[i].idx)
        bamFn = "tmp_anno/{}_{}_TsdToLocal.bam".format(cltView[i].tid, cltView[i].idx)
        if os.path.isfile(bamFn) == True:
            annoTsd(&cltView[i], &annoView[numAnno], numTmp)

        # 4. Define insSeq structure
        setInsStruc(&cltView[i], &annoView[numAnno], numTmp, &classView[0], &sizeView[0], &ltrView[0])

        # 5. Perform post-filtering
        postFilter(&cltView[i])
        
        numAnno += numTmp
    
    annoArr.resize((numAnno,), refcheck=False)
    annoArr.sort(order=['idx', 'queryStart', 'queryEnd'])
    return annoArr


cdef mapTsdToLocal(int tid, int idx):
    cdef str targetFn = "tmp_anno/{}_{}_local.fa".format(tid, idx)
    cdef str queryFn = "tmp_anno/{}_{}_tsd.fa".format(tid, idx)
    cdef str outFn = "tmp_anno/{}_{}_TsdToLocal.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(targetFn, queryFn, outFn)
    if os.path.isfile(queryFn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotateCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]
    
    # 1. Define insertion seq from assembled contig(s)
    annotateAssm(cltView, startIdx, endIdx, cmdArgs)

    # 2. Annotate TE-fragment, PolyA/T, TSD and structure for insSeq
    #    Also perform post-filtering
    cdef object annoArr = annotateIns(cltView, startIdx, endIdx, cmdArgs)
    
    # 3. Output formated cluster and annotation records
    cdef Annotation[::1] annoView = annoArr
    cdef bytes teFn = cmdArgs.teFn.encode()
    cdef bytes refFn = cmdArgs.refFn.encode()
    if (annoView.shape[0] > 0):
        outputAnno(&annoView[0], annoView.shape[0], startIdx, teFn)
    
    if (cltView.shape[0] > 0):
        outputClt(&cltView[0], startIdx, endIdx, refFn, teFn)

    # Output cltArr and annoArr for post-anno-filtering
    cdef object cltArr = np.asarray(cltView)
    annoArr.tofile('tmp_anno/{}_annoArr.dat'.format(startIdx))
    cltArr[startIdx:endIdx].tofile('tmp_anno/{}_cltArr.dat'.format(startIdx))