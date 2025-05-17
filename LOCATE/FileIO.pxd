from .HtslibExternal cimport *
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.errno  cimport errno
from libc.string cimport strerror

cdef class BamFile:
    """
    BamFile(str file_path, str mode, int num_thread=1, BamFile template=None)
    
    A class for reading/writing BAM/FASTA files.
    """

    cdef char       *file_path
    cdef char       *index_file_path
    cdef char       *mode
    cdef int        num_thread
    cdef htsFile    *hts_file
    cdef hts_idx_t  *index
    cdef sam_hdr_t  *header

    cdef void      open_bam_file(self, BamFile template=*)
    cdef htsFile   *_open_hts_file(self) except? NULL
    cdef void      write(self, bam1_t *bam_record)


cdef class Iterator:
    """
    Iterator(BamFile bamFile, int tid)

    A class for iterating over mapped reads in single chromosome.
    """
    cdef bam1_t     *bam_record
    cdef htsFile    *hts_file
    cdef hts_itr_t  *hts_iter
    cdef int64_t    offset

    cdef int next_record_by_tid(self)
    cdef int next_record(self)
    cdef int next_record_by_offset(self, int64_t offset)


cdef extern from "src/AIList.h" nogil:
    ##############
    ### AIList ###
    ##############
    ctypedef struct AiList:
        pass
    
    AiList *initAiList()
    void destroyAiList(AiList *ail)
    void readBED(AiList *ail, const char *bed_fn, const char *targetChrom)
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
        int         ref_start
        int         ref_end
        int         start_idx
        int         end_idx
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
        int         num_thread
        int         tid
        int         min_seg_len
        int         max_dist
        int         overhang
        float       bg_div
        float       bg_depth
        float       bg_read_len
        htsFile     *genome_bam
        bam1_t      *first_bam_record
        bam1_t      *second_bam_record
        AiList      *repeat_ailist
        AiList      *gap_ailist
        AiList      *blacklist_ailist
    
    ########################
    ### Define Candidate ###
    ########################
    int overhang_too_short(Segment *segment, int overhang)


cdef extern from "src/seg_utils.h" nogil:
    ##################
    ### Structures ###
    ##################
    ctypedef packed struct Segment:
        uint8_t     map_qual
        int         query_start
        int         query_end
        int         ref_position
        uint8_t     seg_type
        uint8_t     aln_type
        int64_t     file_offset
        int         aln_refstart
        int         aln_refend
        uint8_t     order
        uint8_t     numSeg
        int         overhang
        int         match_len
        int         read_len
        uint8_t     aln_location_type
        uint8_t     num_te_aln
        int         sum_query_maplen
        float       sum_aln_score
        float       sum_divergence
        int         start_idx
        int         end_idx
    
    ###############################
    ### Initialize TeAlignments ###
    ###############################
    int bam_is_invalid(bam1_t *bam)
    
    ####################
    ### Trim Segment ###
    ####################
    int trim_segment(bam1_t *source_record, bam1_t *dest_record, int seg_idx, int source_start, int source_end)

    #####################
    ### Alinged Pairs ###
    #####################
    int get_aligned_pairs(bam1_t *read, int *query_arr, int *ref_arr)


cdef extern from "src/io_utils.h" nogil:
    ###########################
    ### Segment Sequence IO ###
    ###########################
    int is_lowqual_clt(Cluster *clt)
    int is_lowfreq_clt(Cluster *clt)

    int get_ouput_segidx(Cluster *clt, Segment *seg_arr, Args args)
    void setTrimRegion(Segment *segment, int *start, int *end, int flank_size)

    #########################
    ### Flank Sequence IO ###
    #########################
    void extract_ref_flankseq(char *ref_fn, Cluster *clt_arr, int start_idx, int end_idx)

    #############################
    ### Insertion Sequence IO ###
    #############################
    void extract_ins_seq(Cluster *clt)
    void re_extract_ins_seq(Cluster *clt)


cdef Args new_args(int tid, float bg_div, float bg_depth, float bg_read_len, object cmd_args)
cdef AiList* new_ailist(str bed_fn, const char *chrom)
cdef ouput_seg_seqs(Segment[::1] seg_view, BamFile genome_bam, Args args)
cdef output_highfreq_clusters_seqs(Cluster[::1] clt_view, Segment[::1] seg_view, BamFile genome_bam, Args args)
cpdef output_lowfreq_clusters_seq(Cluster[::1] clt_view, Segment[::1] seg_view, object cmd_args, int tid, int extra_thread)
cpdef int output_read_as_assmbly(Cluster[::1] clt_view, dict cluster_data_by_tid, object cmd_args, int i)
cpdef output_reference_flank(Cluster[::1] clt_view, dict cluster_data_by_tid, tuple block, object cmd_args)
cpdef merge_output()