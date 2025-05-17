from .FileIO cimport *

cdef extern from "src/cluster_utils.h" nogil:
    void update_cluster(Cluster *clt_arr, Segment *seg_arr, Args args)
    void intersect_black_list(Cluster *clt, Args args)


cdef extern from "src/seg_utils.h" nogil:
    #######################
    ### Segment records ###
    #######################
    ctypedef packed struct TeAlignment:
        int     seg_idx
        int     aln_score
        int     query_start
        int     query_end
        int     map_len
        float   divergence
    
    int fill_seg_arr(bam1_t *bam, Segment *seg_arr, int64_t file_offset, int minSegmentLength)
    void update_segment(Segment *seg_arr, AiList *repeat_ailist, AiList *gap_ailist)
    void update_seg_by_te_arr(Segment *seg_arr, TeAlignment *te_arr, int teIdx)

    #########################
    ### Alignment records ###
    #########################
    void get_mapping_length_and_divergence(int *mapLenPtr, float *divergencePtr, bam1_t *bam)
    void fill_te_arr(bam1_t *bam, TeAlignment *te_arr)
    

cpdef dict build_cluster(float bg_div, float bg_depth, float bg_read_len, object cmd_args, int tid, int extra_thread)
cdef object get_high_quality_clusters(dict cluster_data_by_tid)