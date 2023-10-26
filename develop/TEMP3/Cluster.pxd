from .AlignmentFileIO cimport BamFile, Iterator
from .htslib_external cimport *

cdef extern from "src/cluster_utils.h" nogil:
    ctypedef packed struct Cluster:
        int         refStart
        int         refEnd
        int         startIndex
        int         endIndex
        float       numSeg
        uint16_t    directionFlag
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
        float       meanDe
        float       bgDiv
        float       bgDe
        float       bgDepth
        float       bgReadLen
        float       teAlignedFrac
        int         teTid
    
    ctypedef struct Args:
        int         numThread
        int         tid
        int         minSegLen
        int         maxDistance
        int         numTeTid
        int         *teTidCountTable
        int         minOverhang
        float       bgDiv
        float       bgDe
        float       bgDepth
        float       bgReadLen
        htsFile     *genomeBamFile
        bam1_t      *firstBamRecord
        bam1_t      *secondBamRecord
        AiList      *repeatAiList
        AiList      *gapAiList
    
    int isValidCandidate(Cluster *cluster)
    int overhangIsShort(Segment *segment, int minOverhang)
    Args initArgs(int numThread, int tid, int minSegLen, int maxDistance, int minOverhang, float bgDiv, float bgDe, float bgDepth, float bgReadLen)
    void updateCluster(Cluster *cltArray, Segment *segArray, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct Segment:
        uint16_t    flag
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
        uint8_t     alnLocationType
        uint8_t     numTeAlignment
        int         sumQueryMapLen
        float       sumAlnScore
        float       sumDivergence
        float       sumDe
        uint16_t    directionFlag
        int         startIndex
        int         endIndex
        int         teTid

    ctypedef packed struct TeAlignment:
        int     segIndex
        int     AlnScore
        int     queryStart
        int     queryEnd
        int     mapLen
        float   divergence
        float   de
        int16_t flag
        int     teTid
    
    int fillSegmentArray(bam1_t *bamRecord, Segment *segmentArray, int64_t fileOffset, int minSegmentLength)
    void updateSegment(Segment *segArray, AiList *repeatAiList, AiList *gapAiList)
    void updateSegByTeArray(Segment *segArray, TeAlignment *teArray, int teIndex)
    void countTeTids(Segment *segArray, TeAlignment *teArray, int *teTidCountTable, int numTeTid)

    #########################
    ### Alignment records ###
    #########################
    int bamIsInvalid(bam1_t *bamRecord)
    void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bamRecord);
    float getDe(bam1_t *b)
    int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIndex, int sourceStart, int sourceEnd)
    void fillTeArray(bam1_t *bamRecord, TeAlignment *teArray)


cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    ctypedef struct AiList:
        pass

    AiList *initAiList()
    void destroyAiList(AiList *ail)
    void readBED(AiList *ail, const char* bedFileName, const char* targetChrom)
    void constructAiList(AiList *ail, int minCoverageLen)


cpdef dict buildCluster(
    str genomeBamFilePath,
    str repeatPath,
    str gapPath,
    str referenceTe,
    int numThread,
    int tid,
    int minSegLen,
    int maxDistance,
    float bgDiv,
    float bgDe,
    float bgDepth,
    float bgReadLen)