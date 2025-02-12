from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
from cpython cimport PyBytes_Check
from posix.types cimport off_t


########################
### htslib/kstring.h ###
########################
cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s


######################
### htslib/hfile.h ###
######################
cdef extern from "htslib/hfile.h" nogil:
    ctypedef struct hFILE


#####################
### htslib/bgzf.h ###
#####################
cdef extern from "htslib/bgzf.h" nogil:
    ctypedef struct bgzf_mtaux_t
    ctypedef struct bgzidx_t
    ctypedef struct z_stream

    ctypedef struct BGZF:
        unsigned           errcode
        unsigned           is_write
        int           is_be
        int           compress_level
        int           is_compressed
        int           is_gzip
        int           cache_size
        int64_t       block_address
        int64_t       uncompressed_address
        void         *uncompressed_block
        void         *compressed_block
        void         *cache
        hFILE        *fp
        bgzf_mtaux_t *mt
        bgzidx_t     *idx
        int           idx_build_otf
        z_stream     *gz_stream

    int SEEK_SET
    int64_t bgzf_tell(BGZF *fp)
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence)


######################
### htslib/hts.h ###
######################
cdef extern from "htslib/hts.h" nogil:
    ctypedef struct cram_fd

    union FilePointerUnion:
        BGZF    *bgzf
        cram_fd *cram
        hFILE   *hfile
        void    *voidp

    enum htsFormatCategory:
        unknown_category
        sequence_data    # Sequence data -- SAM, BAM, CRAM, etc
        variant_data     # Variant calling data -- VCF, BCF, etc
        index_file       # Index file associated with some data file
        region_list      # Coordinate intervals or regions -- BED, etc
        category_maximum

    enum htsExactFormat:
        unknown_format
        binary_format
        text_format
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed
        format_maximum

    enum htsCompression:
        no_compression, gzip, bgzf, custom
        compression_maximum

    ctypedef struct htsVersion:
        short major, minor

    ctypedef struct htsFormat:
        htsFormatCategory category
        htsExactFormat    format
        htsVersion        version
        htsCompression    compression
        short             compression_level
        void              *specific  

    ctypedef struct htsFile:
        uint8_t  is_bin
        uint8_t  is_write
        uint8_t  is_be
        uint8_t  is_cram
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        FilePointerUnion fp
        htsFormat format

    const char *seq_nt16_str

    int hts_set_threads(htsFile *fp, int n)

    ctypedef struct hts_idx_t

    ctypedef struct hts_pair64_max_t:
        uint64_t u, v
        uint64_t max

    ctypedef struct hts_itr_t:
        int n_off
        hts_pair64_max_t *off
        uint64_t curr_off

    void hts_idx_destroy(hts_idx_t *idx)
    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)


####################
### htslib/sam.h ###
####################
cdef extern from "htslib/sam.h" nogil:
    #**********************
    #*** SAM/BAM header ***
    #**********************

    ctypedef struct sam_hrecs_t

    ctypedef struct sam_hdr_t:
         int32_t n_targets, ignore_sam_err
         uint32_t l_text
         uint32_t *target_len
         uint8_t *cigar_tab
         char **target_name
         char *text
         void *sdict
         sam_hrecs_t *hrecs
         uint32_t ref_count
    
    ctypedef sam_hdr_t bam_hdr_t

    #****************************
    #*** CIGAR related macros ***
    #****************************

    int BAM_CMATCH
    int BAM_CINS
    int BAM_CDEL
    int BAM_CREF_SKIP
    int BAM_CSOFT_CLIP
    int BAM_CHARD_CLIP
    int BAM_CPAD
    int BAM_CEQUAL
    int BAM_CDIFF
    int BAM_CBACK

    char    *BAM_CIGAR_STR
    int      BAM_CIGAR_SHIFT
    uint32_t BAM_CIGAR_MASK
    uint32_t BAM_CIGAR_TYPE

    char bam_cigar_op(uint32_t c)
    uint32_t bam_cigar_oplen(uint32_t c)
    char bam_cigar_opchr(uint32_t)
    uint32_t bam_cigar_gen(char, uint32_t)
    int bam_cigar_type(char o)

    # @abstract the read is paired in sequencing, no matter whether it is mapped in a pair
    int BAM_FPAIRED
    # @abstract the read is mapped in a proper pair
    int BAM_FPROPER_PAIR
    # @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
    int BAM_FUNMAP
    # @abstract the mate is unmapped
    int BAM_FMUNMAP
    # @abstract the read is mapped to the reverse strand
    int BAM_FREVERSE
    # @abstract the mate is mapped to the reverse strand
    int BAM_FMREVERSE
    # @abstract this is read1
    int BAM_FREAD1
    # @abstract this is read2
    int BAM_FREAD2
    # @abstract not primary alignment
    int BAM_FSECONDARY
    # @abstract QC failure
    int BAM_FQCFAIL
    # @abstract optical or PCR duplicate
    int BAM_FDUP
    # @abstract supplementary alignment
    int BAM_FSUPPLEMENTARY

    #*************************
    #*** Alignment records ***
    #*************************

    ctypedef struct bam1_core_t:
        int32_t     pos
        int32_t     tid
        uint16_t    bin
        uint8_t     qual
        uint8_t     l_extranul
        uint16_t    flag
        uint8_t     l_qname
        uint32_t    n_cigar
        int32_t     l_qseq
        int32_t     mtid
        int32_t     mpos
        int32_t     isize

    ctypedef struct bam1_t:
        bam1_core_t core
        uint64_t    id
        uint8_t     *data
        int         l_data
        uint32_t    m_data
        uint32_t    mempolicy

    int bam_is_rev(bam1_t *b)
    char *bam_get_qname(bam1_t *b)
    uint32_t *bam_get_cigar(bam1_t *b)
    char *bam_get_seq(bam1_t *b)
    uint8_t *bam_get_qual(bam1_t *b)
    char bam_seqi(char *s, int i)

    #**************************
    #*** Exported functions ***
    #**************************

    #***************
    #*** BAM I/O ***
    #***************

    void sam_hdr_destroy(sam_hdr_t *h)
    sam_hdr_t* sam_hdr_dup(const sam_hdr_t *h0)
    ctypedef htsFile samFile
    sam_hdr_t *sam_hdr_read(samFile *fp)
    int sam_hdr_write(samFile *fp, const sam_hdr_t *h)
    const char *sam_hdr_tid2name(const sam_hdr_t *h, int tid)
    int sam_hdr_tid2len(const sam_hdr_t *h, int tid)

    ## /* Alignment */ ##
    bam1_t *bam_init1()
    void bam_destroy1(bam1_t *b)
    int bam_read1(BGZF *fp, bam1_t *b)

    #*************************
    #*** BAM/CRAM indexing ***
    #*************************

    hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx)
    void sam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)

    #***************
    #*** SAM I/O ***
    #***************

    htsFile *sam_open(const char *fn, const char *mode)
    int sam_close(htsFile *fp)
    int sam_read1(htsFile *fp, sam_hdr_t *h, bam1_t *b)
    int sam_write1(samFile *fp, const sam_hdr_t *h, const bam1_t *b)

    #*************************************
    #*** Manipulating auxiliary fields ***
    #*************************************

    uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
    int64_t  bam_aux2i(const uint8_t *s)
    double   bam_aux2f(const uint8_t *s)


######################
### htslib/faidx.h ###
######################
cdef extern from "htslib/faidx.h" nogil:
    ctypedef struct faidx_t

    faidx_t *fai_load(const char *fn)
    void fai_destroy(faidx_t *fai)
    int faidx_nseq(const faidx_t *fai)
    const char *faidx_iseq(const faidx_t *fai, int i)
    char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
    int faidx_seq_len(const faidx_t *fai, const char *seq)


