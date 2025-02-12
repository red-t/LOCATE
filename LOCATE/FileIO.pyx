import os
import subprocess
from collections import Counter

########################
### Class Definition ###
########################
cdef int MAX_POS = (1 << 31) - 1
cdef class BamFile:
    def __cinit__(self, str filePath, str mode, int numThread=1, BamFile template=None):
        cdef bytes filePathBytes = os.fsencode(filePath)
        cdef bytes indexPathBytes = os.fsencode(filePath + ".bai")
        cdef bytes modeBytes = mode.encode()

        self.filePath = filePathBytes
        self.indexFilePath = indexPathBytes
        self.mode = modeBytes
        self.numThread = numThread

        self.openFile(template)
    
    cdef void openFile(self, BamFile template=None):
        self.htsFile = self.openHtsFile()
        if self.htsFile == NULL:
            raise IOError("could not open BAM file `{}`".format(self.filePath.decode()))

        if self.mode == b'rb':
            with nogil:
                self.header = sam_hdr_read(self.htsFile)
            if self.header == NULL:
                raise IOError("file does not have a valid header, is it BAM format?")
            
            if os.path.exists(self.indexFilePath.decode()):
                with nogil:
                    self.index = sam_index_load2(self.htsFile, self.filePath, self.indexFilePath)
                if not self.index:
                    raise IOError('unable to open index file `{}`'.format(self.indexFilePath.decode()))

        elif (self.mode == b'wb') or (self.mode == b'wF'):
            if template:
                self.header = sam_hdr_dup(template.header)
            else:
                raise ValueError("need a template for copying header")

            if self.mode == b'wb':
                with nogil:
                    sam_hdr_write(self.htsFile, self.header)
    
    cdef htsFile *openHtsFile(self) except? NULL:
        '''open file in 'rb/wb/wF' mode'''
        cdef htsFile *htsFile
        with nogil:
            htsFile = sam_open(self.filePath, self.mode)
            if htsFile != NULL:
                hts_set_threads(htsFile, self.numThread)
            return htsFile
    
    cdef void write(self, bam1_t *bam):
        cdef int retValue
        with nogil:
            retValue = sam_write1(self.htsFile, self.header, bam)
        if retValue < 0:
            raise IOError("sam_write1 failed with error code {}".format(retValue))
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def close(self):
        if self.htsFile:
            sam_close(self.htsFile)
            self.htsFile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.header:
            sam_hdr_destroy(self.header)
            self.header = NULL

    def __dealloc__(self):
        if self.htsFile:
            sam_close(self.htsFile)
            self.htsFile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.header:
            sam_hdr_destroy(self.header)
            self.header = NULL


cdef class Iterator:
    def __cinit__(self):
        self.bamRcord = bam_init1()
        if self.bamRcord == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self, BamFile bamFile, int tid=-1):
        self.htsFile = bamFile.htsFile
        if tid >= 0:
            with nogil:
                self.iter = sam_itr_queryi(bamFile.index, tid, 0, MAX_POS)
    
    def __dealloc__(self):
        if self.bamRcord:
            bam_destroy1(self.bamRcord)
            self.bamRcord = NULL
        if self.iter:
            sam_itr_destroy(self.iter)
            self.iter = NULL
    
    def __iter__(self):
        return self

    cdef int cnext1(self):
        '''cversion of iterator. retValue>=0 if success.'''
        cdef int retValue

        with nogil:
            if self.iter.curr_off == 0:
                if self.iter.n_off > 0:
                    self.offset = self.iter.off[0].u
                else:
                    self.offset = self.iter.curr_off
            else:
                self.offset = self.iter.curr_off

            retValue = hts_itr_next(self.htsFile.fp.bgzf, self.iter, self.bamRcord, self.htsFile)

        return retValue
    
    cdef int cnext2(self):
        '''directly read a alignment record'''
        cdef int retValue

        with nogil:
            retValue = bam_read1(self.htsFile.fp.bgzf, self.bamRcord)
        
        return retValue

    cdef int cnext3(self, int64_t offset):
        '''read a alignment record with specified offset'''
        cdef int retValue

        with nogil:
            bgzf_seek(self.htsFile.fp.bgzf, offset, SEEK_SET)
            retValue = bam_read1(self.htsFile.fp.bgzf, self.bamRcord)
        
        return retValue


######################
### Construct Args ###
######################
cdef Args newArgs(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs):
    cdef Args args

    args.tid = tid
    args.bgDiv = bgDiv
    args.bgDepth = bgDepth
    args.bgReadLen = bgReadLen
    args.numThread = cmdArgs.numThread
    args.minSegLen = cmdArgs.minSegLen
    args.maxDistance = cmdArgs.maxDistance
    args.minOverhang = cmdArgs.minOverhang
    return args


#######################
### Construc AiList ###
#######################
cdef AiList* newAiList(str bedFn, const char *chrom):
    cdef bytes bedFnBytes = bedFn.encode()
    cdef AiList *aiList = initAiList()

    # Load intervals only if bedFn exists
    if bedFn:
        readBED(aiList, bedFnBytes, chrom)

    constructAiList(aiList, 20)
    return aiList


###########################
### Segment Sequence IO ###
###########################
cdef ouputAllSegSeqs(Segment[::1] segView, BamFile genomeBam, Args args):
    cdef str outputFn = "tmp_build/all_seg_{}.fa".format(args.tid)
    cdef BamFile outputFa = BamFile(outputFn, "wF", args.numThread, genomeBam)
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, retValue

    for i in range(segView.shape[0]):
        retValue = iterator.cnext3(segView[i].fileOffset)
        trimSegment(iterator.bamRcord, destRecord, i, segView[i].queryStart, segView[i].queryEnd)
        outputFa.write(destRecord)
    
    bam_destroy1(destRecord); outputFa.close(); del outputFa; del iterator


cdef outputGermCltSeqs(Cluster[::1] cltView, Segment[::1] segView, BamFile genomeBam, Args args):
    cdef str outputFn
    cdef BamFile outputFa
    cdef Iterator iterator = Iterator(genomeBam, args.tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef int i, j

    for i in range(cltView.shape[0]):
        if isLowQualClt(&cltView[i]) or isSomaClt(&cltView[i]):
            continue

        outputFn = "tmp_assm/{}_{}.fa".format(args.tid, i)
        outputFa = BamFile(outputFn, "wF", args.numThread, genomeBam)
        for j in range(cltView[i].startIdx, cltView[i].endIdx):
            if overhangIsShort(&segView[j], args.minOverhang):
                continue
            outputSingleSeq(segView, outputFa, iterator, destRecord, j)
        
        outputFa.close()

    bam_destroy1(destRecord); del iterator


cpdef outputSomaCltSeqs(Cluster[::1] cltView, Segment[::1] segView, object cmdArgs, int tid):
    cdef int i, j
    cdef Args args
    cdef str outputFn
    cdef BamFile outputFa
    cdef BamFile genomeBam = BamFile(cmdArgs.genomeBamFn, "rb", cmdArgs.numThread)
    cdef Iterator iterator = Iterator(genomeBam, tid)
    cdef bam1_t *destRecord = bam_init1()

    args.minOverhang = cmdArgs.minOverhang
    for i in range(cltView.shape[0]):
        if isLowQualClt(&cltView[i]):
            continue
        
        # Skip successfully assembled clusters
        outputFn = "tmp_assm/{}_{}_assembled.fa".format(tid, i)
        if os.path.isfile(outputFn):
            if os.path.getsize(outputFn) != 0:
                continue
        
        outputFa = BamFile(outputFn, "wF", cmdArgs.numThread, genomeBam)
        j = getOuputSegIdx(&cltView[i], &segView[0], args)
        outputSingleSeq(segView, outputFa, iterator, destRecord, j)
        outputFa.close()

    bam_destroy1(destRecord); del iterator
    genomeBam.close(); del genomeBam


cpdef int outputReadAsAssm(Cluster[::1] cltView, dict allCltData, object cmdArgs, int i):
    if cltView[i].numSegRaw < 3:
        return 0

    cdef Segment[::1] segView = allCltData[cltView[i].tid][1]
    cdef BamFile genomeBam = BamFile(cmdArgs.genomeBamFn, "rb", cmdArgs.numThread)
    cdef Iterator iterator = Iterator(genomeBam, cltView[i].tid)
    cdef bam1_t *destRecord = bam_init1()
    cdef str outputFn = "tmp_assm/{}_{}_assm.fa".format(cltView[i].tid, cltView[i].idx)
    cdef BamFile outputFa = BamFile(outputFn, "wF", cmdArgs.numThread, genomeBam)
    cdef Args args

    args.minOverhang = cmdArgs.minOverhang
    j = getOuputSegIdx(&cltView[i], &segView[0], args)
    outputSingleSeq(segView, outputFa, iterator, destRecord, j)
    
    outputFa.close()
    bam_destroy1(destRecord); del iterator
    genomeBam.close(); del genomeBam

    return 1


cdef outputSingleSeq(Segment[::1] segView, BamFile outputFa, Iterator iterator, bam1_t *destRecord, int j, int flankSize=3000):
    cdef int start, end, retValue

    retValue = iterator.cnext3(segView[j].fileOffset)
    setTrimRegion(&segView[j], &start, &end, flankSize)
    trimSegment(iterator.bamRcord, destRecord, j, start, end)
    outputFa.write(destRecord)


#########################
### Flank Sequence IO ###
#########################
cpdef outputRefFlank(Cluster[::1] cltView, dict allCltData, int startIdx, int taskSize, object cmdArgs):
    cdef int endIdx = startIdx + taskSize
    if endIdx > cltView.shape[0]:
        endIdx = cltView.shape[0]
    
    # Reset breakpoint to most common site
    cdef object counter = Counter()
    cdef Segment[::1] segView
    for i in range(startIdx, endIdx):
        segView = allCltData[cltView[i].tid][1]
        for j in range(cltView[i].startIdx, cltView[i].endIdx):
            if overhangIsShort(&segView[j], cmdArgs.minOverhang):
                continue
            counter[segView[j].refPosition] += 1
        
        cltView[i].refStart = counter.most_common(1)[0][0]
        cltView[i].refEnd = cltView[i].refStart + 1
        counter.clear()
    
    cdef bytes refFn = cmdArgs.refFn.encode('utf-8')
    extractRefFlanks(refFn, &cltView[0], startIdx, endIdx)


####################
### Merge Output ###
####################
cpdef mergeOutput():
    cdef str cmd = "cat tmp_anno/*cltFormated.txt | sort -k1,1 > tmp_anno/sorted_cltFormated.txt && " \
          "cat tmp_anno/*annoFormated.txt | sort -k1,1 > tmp_anno/sorted_annoFormated.txt && " \
          "join -1 1 -2 1 -t $'\t' -o '1.2 1.3 1.4 2.3 1.5 2.2 2.1 2.4 2.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15' " \
          "tmp_anno/sorted_cltFormated.txt tmp_anno/sorted_annoFormated.txt > tmp_anno/result.txt"
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')