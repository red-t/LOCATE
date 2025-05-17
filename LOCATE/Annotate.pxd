from .FileIO cimport *

cdef extern from "src/io_utils.h":
    ###################
    ### Cluster I/O ###
    ###################
    void output_clusters(Cluster *clt_arr, int start_idx, int end_idx, const char *ref_fn, const char *te_fn)

cdef extern from "src/anno_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Annotation:
        int         idx
        int         cltTid
        int         cltIdx
        int         query_start
        int         query_end
        uint8_t     strand
        int         tid
        int         ref_start
        int         ref_end
        uint32_t    flag
        int         extra
    
    ###################################
    ### Annotate Insertion sequence ###
    ###################################
    int fill_anno_arr(Cluster *clt, Annotation *anno_arr, uint32_t *class_arr, int idx)
    void annotate_tsd(Cluster *clt, Annotation *anno_arr, int num_anno)
    void set_ins_structure(Cluster *clt, Annotation *anno_arr, int num_anno, uint32_t *class_arr, int *size_arr, int *ltr_arr)

    ######################
    ### Annotation I/O ###
    ######################
    void output_annotations(Annotation *anno_arr, int num_anno, int start_idx, const char *te_fn)

cdef extern from "src/post_filter.h":
    ###################
    ### Cluster I/O ###
    ###################
    void post_filter(Cluster *clt)

cpdef annotate_cluster(Cluster[::1] clt_view, tuple block, object cmd_args)