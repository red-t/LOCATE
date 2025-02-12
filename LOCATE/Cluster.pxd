from .FileIO cimport *

cdef extern from "src/cluster_utils.h" nogil:
    void updateCluster(Cluster *cltArr, Segment *segArr, Args args)
    void intersectBlackList(Cluster *clt, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct TeAlignment:
        int     segIdx
        int     AlnScore
        int     queryStart
        int     queryEnd
        int     mapLen
        float   divergence
    
    int fillSegArr(bam1_t *bam, Segment *segArr, int64_t fileOffset, int minSegmentLength)
    void updateSegment(Segment *segArr, AiList *repeatAiList, AiList *gapAiList)
    void updateSegByTeArr(Segment *segArr, TeAlignment *teArr, int teIdx)

    #########################
    ### Alignment records ###
    #########################
    void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bam)
    void fillTeArr(bam1_t *bam, TeAlignment *teArr)
    

cpdef dict buildCluster(float bgDiv, float bgDepth, float bgReadLen, object cmdArgs, int tid)
cdef object getHighQualClts(dict allCltData)