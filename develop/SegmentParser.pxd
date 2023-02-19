from htslib_external cimport *
from pysam.libcutils cimport force_bytes
from cpython cimport PyBytes_FromStringAndSize

cdef class InsertSegment:
    """A class for insert/clip segment"""
    cdef bam1_t *_delegate
    cdef public:
        int32_t qstart, qend, rpos  # [qstart, qend) region of query sequence; rpos=="rend", rstart=rpos-1
        uint8_t orient              # left-clip: 0, mid-insert: 1, right-clip: 2


cdef InsertSegment makeInsertSegment(bam1_t *src, int32_t qstart, int32_t qend,
                                     int32_t rpos, uint8_t orient)


cdef int parse_cigar(bam1_t *src, list tmp_segl, uint8_t minl=*)