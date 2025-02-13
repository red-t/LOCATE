import os
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from .FileIO import outputSomaCltSeqs, outputRefFlank, mergeOutput
from .Cluster import buildCluster
from .Assemble import assembleCluster
from .Annotate import annotateCluster


###########################
### Get Background Info ###
###########################
cdef dict getBackgroundInfo(str genomeBamFn, int numThread):
    cdef BamFile genomeBam = BamFile(genomeBamFn, "rb", numThread)
    cdef int i, tid = 0, maxChromLen = 0

    for i in range(genomeBam.header.n_targets):
        if maxChromLen < sam_hdr_tid2len(genomeBam.header, i):
            maxChromLen = sam_hdr_tid2len(genomeBam.header, i)
            tid = i
    
    cdef Iterator iterator = Iterator(genomeBam, tid)
    cdef int numAln = 0, maxNumAln = 499900
    cdef object readLenArr = np.zeros(500000, dtype=np.int32)
    cdef object divArr = np.zeros(500000, dtype=np.float32)
    cdef float[::1] divView = divArr
    cdef int[::1] readLenView = readLenArr
    cdef int retValue, alnLen
    cdef int64_t sumAlnLen = 0
    cdef float divergence
    cdef dict bgInfo = {}

    while True:
        retValue = iterator.cnext1()
        if retValue < 0:
            bgInfo["numChrom"] = genomeBam.header.n_targets
            bgInfo["bgDiv"] = np.mean(divArr[:numAln])
            bgInfo["bgDepth"] = float(sumAlnLen) / maxChromLen
            bgInfo["bgReadLen"] = np.median(readLenArr[:numAln])
            genomeBam.close(); del iterator; del genomeBam; return bgInfo

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numAln >= maxNumAln:
            maxNumAln = divArr.shape[0] + 500000
            divArr.resize((maxNumAln,), refcheck=False)
            readLenArr.resize((maxNumAln,), refcheck=False)
            divView = divArr
            readLenView = readLenArr
            maxNumAln -= 100
            
        getMapLenAndDiv(&alnLen, &divergence, iterator.bamRcord)
        sumAlnLen += alnLen
        divView[numAln] = divergence
        readLenView[numAln] = iterator.bamRcord.core.l_qseq
        numAln += 1


cdef buildTERef(object cmdArgs):
    with open("tmp_build/tmp.fa", "w") as outFa:
        outFa.write('>0\nAAAAAAAAAAAAAA\n')
    
    cdef str queryFn = "tmp_build/tmp.fa"
    cdef str outFn = "tmp_build/tmp.bam"
    cdef str cmd = "minimap2 -t 1 -aY {} {} | samtools view -bhS -o {} -".format(cmdArgs.teFn, queryFn, outFn)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


#######################
### Parallel Module ###
#######################
cdef tuple divideTask(int numTask, int poolSize):
    cdef int i, taskSize = int(numTask / poolSize)
    cdef list startList = []

    if taskSize == 0:
        taskSize = 1
    
    for i in range(0, numTask, taskSize):
        startList.append(i)
    
    return taskSize, startList


cpdef object runInParallel(object cmdArgs):
    cdef set subProcTup
    cdef dict allCltData = {}, tidToCltData, bgInfo
    cdef object subProc, assembleArr, retValue
    cdef int startIdx, taskSize
    cdef list startList

    # 1. Get Background Info
    bgInfo = getBackgroundInfo(cmdArgs.genomeBamFn, cmdArgs.numThread)
    print("Background divergence: {}\nBackground depth: {}\nBackground read length: {}" \
          "".format(bgInfo["bgDiv"], bgInfo["bgDepth"], bgInfo["bgReadLen"]))
    
    # 2. Define LTR size
    cdef bytes teFn = cmdArgs.teFn.encode()
    cdef bytes classFn = cmdArgs.classFn.encode()
    defineLTR(teFn, classFn)
    
    # 3. Build TE reference
    buildTERef(cmdArgs)

    with ProcessPoolExecutor(max_workers=cmdArgs.numProcess) as executor:
        # 4. Build Cluster
        subProcTup = set([executor.submit(buildCluster, \
                                          bgInfo["bgDiv"], \
                                          bgInfo["bgDepth"], \
                                          bgInfo["bgReadLen"], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
        for subProc in as_completed(subProcTup):
            tidToCltData = subProc.result()
            allCltData = {**allCltData, **tidToCltData}
        
        # 5. Local Assembly
        highQualArr = getHighQualClts(allCltData)
        taskSize, startList = divideTask(highQualArr.shape[0], cmdArgs.numProcess)
        subProcTup = set([executor.submit(assembleCluster, \
                                          highQualArr, \
                                          allCltData, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 6. Output sequence for clusters without assembly
        subProcTup = set([executor.submit(outputSomaCltSeqs, \
                                          allCltData[tid][0], \
                                          allCltData[tid][1], \
                                          cmdArgs, tid) for tid in range(bgInfo["numChrom"])])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 7. Output reference flank sequence for high-qual clusters
        subProcTup = set([executor.submit(outputRefFlank, \
                                          highQualArr, \
                                          allCltData, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 8. Annotate clusters
        subProcTup = set([executor.submit(annotateCluster, \
                                          highQualArr, \
                                          startIdx, \
                                          taskSize, \
                                          cmdArgs) for startIdx in startList])
        for subProc in as_completed(subProcTup):
            retValue = subProc.result()
        
        # 9. Merge Output
        mergeOutput()
    
    return allCltData