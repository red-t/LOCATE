from .HtslibExternal cimport *
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.errno  cimport errno
from libc.string cimport strerror

cdef class BamFile:
    """BamFile(str filePath, str mode, int numThread=1, BamFile template=None)"""

    cdef char       *filePath
    cdef char       *indexFilePath
    cdef char       *mode
    cdef int        numThread
    cdef htsFile    *htsFile
    cdef hts_idx_t  *index
    cdef sam_hdr_t  *header

    cdef void      openFile(self, BamFile template=*)
    cdef htsFile   *openHtsFile(self) except? NULL
    cdef void      write(self, bam1_t *bam)


cdef class Iterator:
    """Iterator(BamFile bamFile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef bam1_t     *bamRcord
    cdef htsFile    *htsFile
    cdef hts_itr_t  *iter
    cdef int64_t    offset

    cdef int cnext1(self)
    cdef int cnext2(self)
    cdef int cnext3(self, int64_t offset)


cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    ctypedef struct AiList:
        pass
    
    AiList *initAiList()
    void destroyAiList(AiList *ail)
    void readBED(AiList *ail, const char *bedFn, const char *targetChrom)
    void constructAiList(AiList *ail, int minCoverageLen)


cdef extern from "src/cluster_utils.h" nogil:
    ##############################
    ### Cluster related macros ###
    #############################
    uint32_t CLT_PASS
    uint32_t CLT_ASSEMBLED
    uint32_t CLT_LEFT_FLANK_MAP
    uint32_t CLT_RIGHT_FLANK_MAP
    uint32_t CLT_DIFF_FLANK_MAP
    uint32_t CLT_SAME_FLANK_MAP
    uint32_t CLT_TE_MAP
    uint32_t CLT_POLYA
    uint32_t CLT_TSD
    uint32_t CLT_5P_FULL
    uint32_t CLT_3P_FULL
    uint32_t CLT_SINGLE_TE
    uint32_t CLT_LARGE_GAP
    uint32_t CLT_SELF_TO_SELF
    uint32_t CLT_DNA
    uint32_t CLT_LTR
    uint32_t CLT_LINE
    uint32_t CLT_SINE
    uint32_t CLT_RETROPOSON

    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Cluster:
        int         tid
        int         idx
        int         refStart
        int         refEnd
        int         startIdx
        int         endIdx
        float       numSeg
        uint8_t     cltType
        uint8_t     locationType
        uint8_t     numSegType
        float       entropy
        float       balanceRatio
        float       lowMapQualFrac
        float       dualClipFrac
        float       alnFrac1
        float       alnFrac2
        float       alnFrac4
        float       alnFrac8
        float       alnFrac16
        float       meanMapQual
        float       meanAlnScore
        float       meanQueryMapFrac
        float       meanDivergence
        float       bgDiv
        float       bgDepth
        float       bgReadLen
        float       teAlignedFrac
        uint8_t     isInBlacklist
        float       probability
        uint32_t    flag
        int         numSegRaw
        int         numLeft
        int         numMiddle
        int         numRight
        int         tid1
        int         leftMost
        int         tid2
        int         rightMost
        int         insLen
        int         repTid
    
    ctypedef struct Args:
        int         numThread
        int         tid
        int         minSegLen
        int         maxDistance
        int         minOverhang
        float       bgDiv
        float       bgDepth
        float       bgReadLen
        htsFile     *genomeBam
        bam1_t      *firstBamRecord
        bam1_t      *secondBamRecord
        AiList      *repeatAiList
        AiList      *gapAiList
        AiList      *blackAiList
    
    ########################
    ### Define Candidate ###
    ########################
    int overhangIsShort(Segment *segment, int minOverhang)


cdef extern from "src/seg_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Segment:
        uint8_t     mapQual
        int         queryStart
        int         queryEnd
        int         refPosition
        uint8_t     segType
        uint8_t     alnType
        int64_t     fileOffset
        int         alnRefStart
        int         alnRefEnd
        uint8_t     order
        uint8_t     numSeg
        int         overhang
        int         matchLen
        int         readLen
        uint8_t     alnLocationType
        uint8_t     numTeAlignment
        int         sumQueryMapLen
        float       sumAlnScore
        float       sumDivergence
        int         startIdx
        int         endIdx
    
    ###############################
    ### Initialize TeAlignments ###
    ###############################
    int bamIsInvalid(bam1_t *bam)
    
    ####################
    ### Trim Segment ###
    ####################
    int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIdx, int sourceStart, int sourceEnd)

    #####################
    ### Alinged Pairs ###
    #####################
    int getAlignedPairs(bam1_t *read, int *queryArr, int *refArr)


cdef extern from "src/io_utils.h" nogil:
    ###########################
    ### Segment Sequence IO ###
    ###########################
    int isLowQualClt(Cluster *clt)
    int isSomaClt(Cluster *clt)

    int getOuputSegIdx(Cluster *clt, Segment *segArr, Args args)
    void setTrimRegion(Segment *segment, int *start, int *end, int flankSize)

    #########################
    ### Flank Sequence IO ###
    #########################
    void extractRefFlanks(char *refFn, Cluster *cltArr, int startIdx, int endIdx)

    #############################
    ### Insertion Sequence IO ###
    #############################
    void extractIns(Cluster *clt)
    void reExtractIns(Cluster *clt)


cdef Args newArgs(int tid, float bgDiv, float bgDepth, float bgReadLen, object cmdArgs)
cdef AiList* newAiList(str bedFn, const char *chrom)
cdef ouputAllSegSeqs(Segment[::1] segView, BamFile genomeBam, Args args)
cdef outputGermCltSeqs(Cluster[::1] cltView, Segment[::1] segView, BamFile genomeBam, Args args)
cpdef outputSomaCltSeqs(Cluster[::1] cltView, Segment[::1] segView, object cmdArgs, int tid)
cpdef int outputReadAsAssm(Cluster[::1] cltView, dict allCltData, object cmdArgs, int i)
cpdef outputRefFlank(Cluster[::1] cltView, dict allCltData, int startIdx, int taskSize, object cmdArgs)
cpdef mergeOutput()