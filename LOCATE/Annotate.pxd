from .FileIO cimport *

cdef extern from "src/io_utils.h":
    ###################
    ### Cluster I/O ###
    ###################
    void outputClt(Cluster *cltArr, int startIdx, int endIdx, const char *refFn, const char *teFn)

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Annotation:
        int         idx
        int         cltTid
        int         cltIdx
        int         queryStart
        int         queryEnd
        uint8_t     strand
        int         tid
        int         refStart
        int         refEnd
        uint32_t    flag
    
    ###################################
    ### Annotate Insertion sequence ###
    ###################################
    int fillAnnoArr(Cluster *clt, Annotation *annoArr, uint32_t *classArr, int idx)
    void annoTsd(Cluster *clt, Annotation *annoArr, int numAnno)
    void setInsStruc(Cluster *clt, Annotation *annoArr, int numAnno, uint32_t *classArr, int *sizeArr, int *ltrArr)

    ######################
    ### Annotation I/O ###
    ######################
    void outputAnno(Annotation *annoArr, int numAnno, int startIdx, const char *teFn)

cdef extern from "src/post_filter.h":
    ###################
    ### Cluster I/O ###
    ###################
    void postFilter(Cluster *clt)

cpdef annotateCluster(Cluster[::1] cltView, int startIdx, int taskSize, object cmdArgs)