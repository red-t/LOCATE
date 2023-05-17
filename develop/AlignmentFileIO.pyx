import os
import numpy as np


SEG_DTYPE = np.dtype([
    ('flag',        np.uint16),
    ('mapq',        np.uint8),
    ('qst',         np.int32),
    ('qed',         np.int32),
    ('rpos',        np.int32),
    ('sflag',       np.uint8),
    ('rflag',       np.uint8),
    ('offset',      np.int64),
    ('refst',       np.int32),
    ('refed',       np.int32),
    ('ith',         np.uint8),
    ('nseg',        np.uint8),
    ('overhang',    np.int32),
    ('nmatch',      np.int32),
    ('loc_flag1',   np.uint8),
    ('loc_flag2',   np.uint8),
])
#
# ---------------------------------------------------------------
#
cdef class BamFile:
    def __cinit__(self,
                  str     filepath,
                  int32_t nthreads = 1,
                  str     mode     = 'rb',
                  BamFile template = None):
        cdef bytes bfilepath    = os.fsencode(filepath)
        cdef bytes bidx         = os.fsencode(filepath+".bai")
        cdef bytes bmode        = mode.encode(TEXT_ENCODING, ERROR_HANDLER)

        self.threads    = nthreads
        self.mode       = bmode
        self.filename   = bfilepath
        self.index_filename = bidx
        self._open(template)
        
    def close(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            bam_hdr_destroy(self.hdr)
            self.hdr = NULL

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.hdr:
            bam_hdr_destroy(self.hdr)
            self.hdr = NULL

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    cdef void _open(self, BamFile template=None):
        '''open BAM file for reading/writting'''
        # open file as htsFile 
        self.htsfile = self._open_htsfile()
        if self.htsfile == NULL:
            if errno:
                raise IOError(errno, "could not open alignment file `{}`: {}".format(
                    self.filename.decode(TEXT_ENCODING, ERROR_HANDLER),
                    strerror(errno)))
            else:
                raise ValueError("could not open alignment file `{}`".format(
                    self.filename.decode(TEXT_ENCODING, ERROR_HANDLER)))

        # for reading
        if self.mode == b'rb':
            # bam files require a valid header
            with nogil:
                self.hdr = sam_hdr_read(self.htsfile)
            if self.hdr == NULL:
                raise ValueError("file does not have a valid header, is it BAM format?")
            
            # open BAM index file to enable random access
            with nogil:
                self.index = sam_index_load2(self.htsfile, self.filename, self.index_filename)
            if not self.index:
                if errno:
                    raise IOError(errno, strerror(errno))
                else:
                    raise IOError('unable to open index file `{}`'.format(
                        self.index_filename.decode(TEXT_ENCODING, ERROR_HANDLER)))
        # for writing
        elif self.mode == b'wb':
            # copy header from template
            if template:
                self.hdr = bam_hdr_dup(template.hdr)
            else:
                raise ValueError("not enough information to construct header. Please provide template")
            
            # write header to htsfile
            with nogil:
                sam_hdr_write(self.htsfile, self.hdr)
    
    cdef htsFile *_open_htsfile(self) except? NULL:
        '''open file in 'rb/wb' mode, return htsFile object if success.'''
        cdef int32_t threads = self.threads
        cdef htsFile *htsfile

        if isinstance(self.filename, bytes):
            with nogil:
                htsfile = hts_open(self.filename, self.mode)
                if htsfile != NULL:
                    hts_set_threads(htsfile, threads)
                return htsfile

    cpdef object extract_seg(self,
                             BamFile wbf,
                             int tid,
                             int minl=50):
        '''extract segments from chromosome tid
        
        Usage: extract_seg(0, 50) to extract segments with length >= 50
            from chromosome 0.

        Parameters
        ----------
        wbf: BamFile
            BamFile object opened for writting.

        tid: int
            tid of the target chromosome.

        minl: int
            minimum length of the segment. segment with cigar_length < minl
            will not be used to create a new InsertSegment.
        
        Returns
        -------
        segs: object
		    strctured numpy array with `dtype=SEG_DTYPE`.
        '''
        cdef Iterator ite

        ite = Iterator(self, wbf, tid, minl)
        for i in ite:
            continue

        return ite.segs
    
    cdef void write(self, bam1_t *src):
        '''write a single alignment to disk.

        Parameters
        ----------
        src: bam1_t*
            bam1_t type pointer, which point to the source memory address
            of the alignment loaded into memory.
        '''
        cdef int ret

        with nogil:
            ret = sam_write1(self.htsfile, self.hdr, src)
        if ret < 0:
            raise IOError(
            "sam_write1 failed with error code {}".format(ret))
#
# ---------------------------------------------------------------
#
cdef class Iterator:
    def __cinit__(self):
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))
        if self.b == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

    def __init__(self,
                 BamFile bamfile,
                 BamFile wbf,
                 int tid,
                 int minl):

        self.bamfile = bamfile
        self.htsfile = bamfile.htsfile
        self.index   = bamfile.index
        self.wbf     = wbf
        self.tid     = tid
        self.minl    = minl
        with nogil:
            self.iter = sam_itr_queryi(self.index,
                                       tid,
                                       0,
                                       MAX_POS)
    
    def __dealloc__(self):
        bam_destroy1(self.b)
        hts_itr_destroy(self.iter)
    
    def __iter__(self):
        return self

    cdef int cnext(self):
        '''cversion of iterator. retval>=0 if success.'''
        cdef int        retval
        cdef int64_t    offset

        with nogil:
            offset = bgzf_tell(self.htsfile.fp.bgzf)
            retval = hts_itr_next(self.htsfile.fp.bgzf,
                                  self.iter,
                                  self.b,
                                  self.htsfile)
        self.offset = offset
        return retval
    
    def __next__(self):
        '''fetch a record and parse it's CIGAR, store results in SEG_DICT'''
        cdef int32_t    n, retval, M
        cdef int32_t    N = 0
        cdef seg_dtype_struct[::1] segs_view

        self.segs = np.zeros(10000, dtype=SEG_DTYPE)
        template  = np.zeros(10000, dtype=SEG_DTYPE)
        segs_view = self.segs
        M = segs_view.shape[0] - 20

        while 1:
            retval = self.cnext()
            # If current iterator is not exhausted, parse the alignment
            if retval > 0:
                n = self.b.core.n_cigar
                if n==0:
                    continue
                if N == M:
                    self.segs = np.concatenate((self.segs, template))
                    segs_view = self.segs
                    M = segs_view.shape[0] - 20

                retval = parse_cigar(self.b, segs_view, N, self.offset, self.minl)
                if retval > 0:
                    # total number of extracted segments
                    N += retval
                    self.wbf.write(self.b)
                continue

            self.segs = self.segs[:N,]
            del template
            raise StopIteration