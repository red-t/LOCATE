import os
import numpy as np
import pandas as pd
import subprocess
from autogluon.tabular import TabularPredictor

##################
### Data Types ###
##################
SegmentDt = np.dtype([
    ('mapQual',         np.uint8),
    ('queryStart',      np.int32),
    ('queryEnd',        np.int32),
    ('refPosition',     np.int32),
    ('segType',         np.uint8),
    ('alnType',         np.uint8),
    ('fileOffset',      np.int64),
    ('alnRefStart',     np.int32),
    ('alnRefEnd',       np.int32),
    ('order',           np.uint8),
    ('numSeg',          np.uint8),
    ('overhang',        np.int32),
    ('matchLen',        np.int32),
    ('readLen',         np.int32),
    ('alnLocationType', np.uint8),
    ('numTeAlignment',  np.uint8),
    ('sumQueryMapLen',  np.int32),
    ('sumAlnScore',     np.float32),
    ('sumDivergence',   np.float32),
    ('startIdx',        np.int32),
    ('endIdx',          np.int32),
])

TeAlignmentDt = np.dtype([
    ('segIdx',      np.int32),
    ('AlnScore',    np.int32),
    ('queryStart',  np.int32),
    ('queryEnd',    np.int32),
    ('mapLen',      np.int32),
    ('divergence',  np.float32),
])

ClusterDt = np.dtype([
    ('tid',                 np.int32),
    ('idx',                 np.int32),
    ('refStart',            np.int32),
    ('refEnd',              np.int32),
    ('startIdx',            np.int32),
    ('endIdx',              np.int32),
    ('numSeg',              np.float32),
    ('cltType',             np.uint8),
    ('locationType',        np.uint8),
    ('numSegType',          np.uint8),
    ('entropy',             np.float32),
    ('balanceRatio',        np.float32),
    ('lowMapQualFrac',      np.float32),
    ('dualClipFrac',        np.float32),
    ('alnFrac1',            np.float32),
    ('alnFrac2',            np.float32),
    ('alnFrac4',            np.float32),
    ('alnFrac8',            np.float32),
    ('alnFrac16',           np.float32),
    ('meanMapQual',         np.float32),
    ('meanAlnScore',        np.float32),
    ('meanQueryMapFrac',    np.float32),
    ('meanDivergence',      np.float32),
    ('bgDiv',               np.float32),
    ('bgDepth',             np.float32),
    ('bgReadLen',           np.float32),
    ('teAlignedFrac',       np.float32),
    ('isInBlacklist',       np.uint8),
    ('probability',         np.float32),
    ('flag',                np.uint32),
    ('numSegRaw',           np.int32),
    ('numLeft',             np.int32),
    ('numMiddle',           np.int32),
    ('numRight',            np.int32),
    ('tid1',                np.int32),
    ('leftMost',            np.int32),
    ('tid2',                np.int32),
    ('rightMost',           np.int32),
    ('insLen',              np.int32),
    ('repTid',              np.int32),
])


########################
### Construct SegArr ###
########################
cdef object getSegArr(BamFile genomeBam, Args args):
    cdef int retValue, numSeg = 0, maxNum = 9900
    cdef object segArr = np.zeros(10000, dtype=SegmentDt)
    cdef Segment[::1] segView = segArr
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    
    while True:
        retValue = iterator.cnext1()
        if retValue < 0:
            del iterator
            segArr.resize((numSeg,), refcheck=False)
            return segArr

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numSeg >= maxNum:
            maxNum = segView.shape[0] + 10000
            segArr.resize((maxNum,), refcheck=False)
            segView = segArr
            maxNum -= 100

        retValue = fillSegArr(iterator.bamRcord, &segView[numSeg], iterator.offset, args.minSegLen)
        numSeg += retValue


cdef updateSegArr(Segment[::1] segView, Args args):
    cdef int i
    for i in range(segView.shape[0]):
        updateSegment(&segView[i], args.repeatAiList, args.gapAiList)


cdef updateSegArrByTe(Segment[::1] segView, Args args):
    cdef BamFile teBam = BamFile("tmp_build/all_seg_{}.bam".format(args.tid), "rb", args.numThread)
    cdef Iterator iterator = Iterator(teBam)
    cdef object teArr = getTeArr(iterator)
    cdef TeAlignment[::1] teView = teArr
    cdef int i

    teArr.sort(order=['segIdx', 'queryStart'])
    for i in range(teView.shape[0]):
        updateSegByTeArr(&segView[0], &teView[0], i)

    del iterator; teBam.close(); del teBam


#########################
### Construct TeArr ###
#########################
cdef mapSegToTE(str teFn, Args args):
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t {0} -aY {1} tmp_build/all_seg_{2}.fa | " \
                   "samtools view -@ {0} -bhS -o tmp_build/all_seg_{2}.bam -".format(args.numThread, teFn, args.tid)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef object getTeArr(Iterator iterator):
    cdef int retValue, numTeAln=0, maxNum=9900
    cdef object teArr = np.zeros(10000, dtype=TeAlignmentDt)
    cdef TeAlignment[::1] teView = teArr
    
    while True:
        retValue = iterator.cnext2()
        if retValue < 0:
            teArr.resize((numTeAln,), refcheck=False)
            return teArr

        if bamIsInvalid(iterator.bamRcord):
            continue

        if numTeAln >= maxNum:
            maxNum = teArr.shape[0] + 10000
            teArr.resize((maxNum,), refcheck=False)
            teView = teArr
            maxNum -= 100

        fillTeArr(iterator.bamRcord, &teView[numTeAln])
        numTeAln += 1


##########################
### Construct CltArr ###
##########################
cdef object getCltArr(Segment[::1] segView, Args args):
    cdef int startIdx=0, endIdx, numClt=0, maxNum=9900
    cdef object cltArr = np.zeros(10000, dtype=ClusterDt)
    cdef Cluster[::1] cltView = cltArr
    
    while startIdx < segView.shape[0]:
        if overhangIsShort(&segView[startIdx], args.minOverhang):
            startIdx += 1; continue

        if numClt > maxNum:
            maxNum = cltArr.shape[0] + 10000
            cltArr.resize((maxNum,), refcheck=False)
            cltView = cltArr
            maxNum -= 100
        
        # Initialize numClt-th cluster
        cltView[numClt].tid = args.tid
        cltView[numClt].refStart = segView[startIdx].refPosition - 1
        cltView[numClt].refEnd = segView[startIdx].refPosition + args.maxDistance
        cltView[numClt].idx = numClt
        cltView[numClt].startIdx = startIdx

        endIdx = startIdx + 1
        while endIdx < segView.shape[0]:
            if overhangIsShort(&segView[endIdx], args.minOverhang):
                endIdx += 1; continue
            if segView[endIdx].refPosition > cltView[numClt].refEnd:
                break

            cltView[numClt].refEnd = segView[endIdx].refPosition + args.maxDistance
            endIdx += 1
        
        cltView[numClt].endIdx = endIdx
        cltView[numClt].refEnd = cltView[numClt].refEnd - args.maxDistance
        startIdx = endIdx; numClt += 1
    
    cltArr.resize((numClt,), refcheck=False)
    return cltArr


cdef updateCltArr(Cluster[::1] cltView, Segment[::1] segView, BamFile genomeBam, Args args):
    cdef int i
    
    args.genomeBam = genomeBam.htsFile
    args.firstBamRecord = bam_init1()
    args.secondBamRecord = bam_init1()

    for i in range(cltView.shape[0]):
        updateCluster(&cltView[i], &segView[0], args)
    
    bam_destroy1(args.firstBamRecord); bam_destroy1(args.secondBamRecord)
    destroyAiList(args.repeatAiList); destroyAiList(args.gapAiList)


######################
### Filter Cluster ###
######################
cdef filterByBlacklist(Cluster[::1] cltView, Args args):
    cdef int i
    for i in range(cltView.shape[0]):
        intersectBlackList(&cltView[i], args)
    
    destroyAiList(args.blackAiList)


cdef object filterByModel(object cltArr, object cmdArgs):
    cdef object cltDf = pd.DataFrame(cltArr)
    
    filterGermByModel(cltDf, cmdArgs.highFreqModel)
    filterSomaByModel(cltDf, cmdArgs.lowFreqModel)

    return cltDf.to_records(index=False)

cdef filterGermByModel(object cltDf, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object germDf = cltDf.loc[(cltDf['cltType']==0) & (cltDf['teAlignedFrac']>=0.8) & (cltDf['isInBlacklist']==0)]

    probability = predictor.predict_proba(germDf)
    probability.columns = ['0', 'probability']
    cltDf.update(probability)

cdef filterSomaByModel(object cltDf, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object somaDf = cltDf.loc[(cltDf['cltType']>0) & (cltDf['teAlignedFrac']>=0.8) & (cltDf['isInBlacklist']==0)]

    probability = predictor.predict_proba(somaDf)
    probability.columns = ['0', 'probability']
    cltDf.update(probability)


#####################
### Build Cluster ###
#####################
cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid):
    # 1. Construct segments
    cdef Args args = newArgs(tid, bgDiv, bgDepth, bgReadLen, cmdArgs)
    cdef BamFile genomeBam = BamFile(cmdArgs.genomeBamFn, "rb", cmdArgs.numThread)
    cdef object segArr = getSegArr(genomeBam, args)

    # 2. Compute segment features
    cdef const char *chrom = sam_hdr_tid2name(genomeBam.header, tid)
    args.repeatAiList = newAiList(cmdArgs.repeatFn, chrom)
    args.gapAiList = newAiList(cmdArgs.gapFn, chrom)

    updateSegArr(segArr, args)
    segArr.sort(order='refPosition')
    ouputAllSegSeqs(segArr, genomeBam, args)

    mapSegToTE(cmdArgs.teFn, args)
    updateSegArrByTe(segArr, args)
    
    # 3. Construct cluster
    cdef tidToCltData = {}
    cdef object cltArr = getCltArr(segArr, args)

    # 4. Compute cluster features
    updateCltArr(cltArr, segArr, genomeBam, args)

    # 5. Filter clusters
    args.blackAiList = newAiList(cmdArgs.blacklistFn, chrom)
    filterByBlacklist(cltArr, args)
    cltArr = filterByModel(cltArr, cmdArgs)

    # 6. Output cluster seqs for assembling
    outputGermCltSeqs(cltArr, segArr, genomeBam, args)

    tidToCltData[tid] = (cltArr, segArr)
    genomeBam.close(); del genomeBam
    return tidToCltData


##############################
### Get High-Qual Clusters ###
##############################
cdef object getHighQualClts(dict allCltData):
    cdef int i, tid, numClt = 0, maxNum = 1900
    cdef object highQualArr = np.zeros(2000, dtype=ClusterDt)
    cdef Cluster[::1] highQualView = highQualArr
    cdef Cluster[::1] cltView

    for tid in range(len(allCltData)):
        cltView = allCltData[tid][0]
        if cltView.shape[0] == 0:
            continue
        
        for i in range(cltView.shape[0]):
            if isLowQualClt(&cltView[i]):
                continue
            
            if numClt >= maxNum:
                maxNum = highQualArr.shape[0] + 2000
                highQualArr.resize((maxNum,), refcheck=False)
                highQualView = highQualArr
                maxNum -= 100
            
            highQualView[numClt] = cltView[i]
            numClt += 1

    highQualArr.resize((numClt,), refcheck=False)
    return highQualArr