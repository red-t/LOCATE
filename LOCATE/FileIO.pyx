import os
import pandas as pd
import logging
import subprocess
from collections import Counter

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

########################
### Constants ###
########################
cdef int MAX_POSITION = (1 << 31) - 1

########################
### BamFile Class ###
########################
cdef class BamFile:
    """
    A class for handling BAM file operations.
    """
    def __cinit__(self, str file_path, str mode, int num_thread=1, BamFile template=None):
        cdef bytes file_path_bytes = os.fsencode(file_path)
        cdef bytes index_path_bytes = os.fsencode(file_path + ".bai")
        cdef bytes mode_bytes = mode.encode()

        self.file_path = file_path_bytes
        self.index_file_path = index_path_bytes
        self.mode = mode_bytes
        self.num_thread = num_thread

        self.open_bam_file(template)
    
    cdef void open_bam_file(self, BamFile template=None):
        """
        Open a BAM file in the specified mode.
        """
        try:
            self.hts_file = self._open_hts_file()
            if self.hts_file == NULL:
                raise IOError(f"Could not open BAM file `{self.file_path.decode()}`")

            if self.mode == b'rb':
                with nogil:
                    self.header = sam_hdr_read(self.hts_file)
                if self.header == NULL:
                    raise IOError("File does not have a valid header, is it BAM format?")
                
                if os.path.exists(self.index_file_path.decode()):
                    with nogil:
                        self.index = sam_index_load2(self.hts_file, self.file_path, self.index_file_path)
                    if not self.index:
                        raise IOError(f"Unable to open index file `{self.index_file_path.decode()}`")

            elif self.mode in (b'wb', b'wF'):
                if template:
                    self.header = sam_hdr_dup(template.header)
                else:
                    raise ValueError("Need a template for copying header")

                if self.mode == b'wb':
                    with nogil:
                        sam_hdr_write(self.hts_file, self.header)
        except Exception as e:
            logging.error(f"Error opening BAM file: {e}")
            raise
    
    cdef htsFile *_open_hts_file(self) except? NULL:
        """
        Open a BAM file in 'rb/wb/wF' mode.
        """
        cdef htsFile *hts_file
        with nogil:
            hts_file = sam_open(self.file_path, self.mode)
            if hts_file != NULL:
                hts_set_threads(hts_file, self.num_thread)
            return hts_file
    
    cdef void write(self, bam1_t *bam_record):
        """
        Write a BAM record to the file.
        """
        cdef int return_value
        with nogil:
            return_value = sam_write1(self.hts_file, self.header, bam_record)
        if return_value < 0:
            raise IOError(f"sam_write1 failed with error code {return_value}")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def close(self):
        """
        Close the BAM file and release resources.
        """
        if self.hts_file:
            sam_close(self.hts_file)
            self.hts_file = NULL
        if self.index:
            hts_idx_destroy(self.index)
            self.index = NULL
        if self.header:
            sam_hdr_destroy(self.header)
            self.header = NULL

    def __dealloc__(self):
        self.close()


########################
### Iterator Class ###
########################
cdef class Iterator:
    """
    A class for iterating over BAM records.
    """
    def __cinit__(self):
        self.bam_record = bam_init1()
        if self.bam_record == NULL:
            raise MemoryError(f"Could not allocate memory of size {sizeof(bam1_t)}")

    def __init__(self, BamFile bam_file, int tid=-1, int beg=0, int end=MAX_POSITION):
        self.hts_file = bam_file.hts_file
        if tid >= 0:
            with nogil:
                self.hts_iter = sam_itr_queryi(bam_file.index, tid, beg, end)
    
    def __dealloc__(self):
        if self.bam_record:
            bam_destroy1(self.bam_record)
            self.bam_record = NULL
        if self.hts_iter:
            sam_itr_destroy(self.hts_iter)
            self.hts_iter = NULL
    
    def __iter__(self):
        return self

    cdef int next_record_by_tid(self):
        """
        Read next alignment record on specified chromosome.
        """
        cdef int return_value
        with nogil:
            if self.hts_iter.curr_off == 0:
                if self.hts_iter.n_off > 0:
                    self.offset = self.hts_iter.off[0].u
                else:
                    self.offset = self.hts_iter.curr_off
            else:
                self.offset = self.hts_iter.curr_off

            return_value = hts_itr_next(self.hts_file.fp.bgzf, self.hts_iter, self.bam_record, self.hts_file)
        return return_value
    
    cdef int next_record(self):
        """
        Directly read next alignment record.
        """
        cdef int return_value
        with nogil:
            return_value = bam_read1(self.hts_file.fp.bgzf, self.bam_record)
        return return_value

    cdef int next_record_by_offset(self, int64_t offset):
        """
        Read an alignment record with a specified offset.
        """
        cdef int return_value
        with nogil:
            bgzf_seek(self.hts_file.fp.bgzf, offset, SEEK_SET)
            return_value = bam_read1(self.hts_file.fp.bgzf, self.bam_record)
        return return_value


########################
### Utility Functions ###
########################
cdef Args new_args(int tid, float bg_div, float bg_depth, float bg_read_len, object cmd_args):
    """
    Construct an Args object with the given parameters.
    """
    cdef Args args
    args.tid = tid
    args.bg_div = bg_div
    args.bg_depth = bg_depth
    args.bg_read_len = bg_read_len
    args.num_thread = cmd_args.num_thread
    args.min_seg_len = cmd_args.min_seg_len
    args.max_dist = cmd_args.max_dist
    args.overhang = cmd_args.overhang
    return args


cdef AiList* new_ailist(str bed_fn, const char *chrom):
    """
    Construct a new AiList object.
    """
    cdef bytes bed_fn_bytes = bed_fn.encode()
    cdef AiList *ai_list = initAiList()

    if bed_fn:
        readBED(ai_list, bed_fn_bytes, chrom)

    constructAiList(ai_list, 20)
    return ai_list


########################
### Output Functions ###
########################
cdef ouput_seg_seqs(Segment[::1] seg_view, BamFile genome_bam, Args args):
    """
    Output all segments' sequences to a file.
    """
    cdef str output_fn = "tmp_build/all_seg_{}.fa".format(args.tid)
    cdef BamFile output_fa = BamFile(output_fn, "wF", args.num_thread, genome_bam)
    cdef Iterator iterator = Iterator(genome_bam, args.tid)
    cdef bam1_t *dest_record = bam_init1()
    cdef int i, return_value

    try:
        for i in range(seg_view.shape[0]):
            return_value = iterator.next_record_by_offset(seg_view[i].file_offset)
            if return_value < 0:
                raise IOError(f"Failed to read record at offset {seg_view[i].file_offset}")
            
            trim_segment(
                iterator.bam_record,
                dest_record, i,
                seg_view[i].query_start,
                seg_view[i].query_end
            )
            output_fa.write(dest_record)
    finally:
        bam_destroy1(dest_record)
        output_fa.close()
        del output_fa
        del iterator


cdef output_highfreq_clusters_seqs(Cluster[::1] clt_view, Segment[::1] seg_view, BamFile genome_bam, Args args):
    """
    Output segment sequences of each high-frequency cluster for assembly.
    """
    cdef str output_fn
    cdef BamFile output_fa
    cdef Iterator iterator = Iterator(genome_bam, args.tid)
    cdef bam1_t *dest_record = bam_init1()
    cdef int i, j

    try:
        for i in range(clt_view.shape[0]):
            if is_lowqual_clt(&clt_view[i]) or is_lowfreq_clt(&clt_view[i]):
                continue

            output_fn = "tmp_assm/{}_{}.fa".format(args.tid, i)
            output_fa = BamFile(output_fn, "wF", args.num_thread, genome_bam)

            for j in range(clt_view[i].start_idx, clt_view[i].end_idx):
                if overhang_too_short(&seg_view[j], args.overhang):
                    continue
                output_single_seq(seg_view, output_fa, iterator, dest_record, j)

            output_fa.close()
    finally:
        bam_destroy1(dest_record)
        del iterator


cpdef output_lowfreq_clusters_seq(Cluster[::1] clt_view, Segment[::1] seg_view, object cmd_args, int tid, int extra_thread):
    """
    Output one segment sequence for each low-frequency cluster as assembly.
    """
    cdef int i, j
    cdef int num_thread = 1 + extra_thread
    cdef Args args
    cdef str output_fn
    cdef BamFile output_fa
    cdef BamFile genome_bam = BamFile(cmd_args.genome_bam_fn, "rb", num_thread)
    cdef Iterator iterator = Iterator(genome_bam, tid)
    cdef bam1_t *dest_record = bam_init1()

    try:
        args.overhang = cmd_args.overhang
        for i in range(clt_view.shape[0]):
            if is_lowqual_clt(&clt_view[i]):
                continue

            # Skip successfully assembled clusters
            output_fn = "tmp_assm/{}_{}_assembled.fa".format(tid, i)
            if os.path.isfile(output_fn) and os.path.getsize(output_fn) != 0:
                continue

            output_fa = BamFile(output_fn, "wF", num_thread, genome_bam)
            j = get_ouput_segidx(&clt_view[i], &seg_view[0], args)
            output_single_seq(seg_view, output_fa, iterator, dest_record, j)
            output_fa.close()
    finally:
        bam_destroy1(dest_record)
        del iterator
        genome_bam.close()
        del genome_bam


cpdef int output_read_as_assmbly(Cluster[::1] clt_view, dict cluster_data_by_tid, object cmd_args, int i):
    """
    Output one segment sequence as assembly for high-frequency cluster that failed to be assembled.
    """
    if clt_view[i].numSegRaw < 3:
        return 0

    cdef Segment[::1] seg_view = cluster_data_by_tid[clt_view[i].tid][1]
    cdef BamFile genome_bam = BamFile(cmd_args.genome_bam_fn, "rb", cmd_args.num_thread)
    cdef Iterator iterator = Iterator(genome_bam, clt_view[i].tid)
    cdef bam1_t *dest_record = bam_init1()
    cdef str output_fn = "tmp_assm/{}_{}_assm.fa".format(clt_view[i].tid, clt_view[i].idx)
    cdef BamFile output_fa = BamFile(output_fn, "wF", cmd_args.num_thread, genome_bam)
    cdef Args args

    try:
        args.overhang = cmd_args.overhang
        j = get_ouput_segidx(&clt_view[i], &seg_view[0], args)
        output_single_seq(seg_view, output_fa, iterator, dest_record, j)
    finally:
        output_fa.close()
        bam_destroy1(dest_record)
        del iterator
        genome_bam.close()
        del genome_bam

    return 1


cdef output_single_seq(Segment[::1] seg_view, BamFile output_fa, Iterator iterator, bam1_t *dest_record, int j, int flank_size=3000):
    """
    Output a single sequence to the output file.
    """
    cdef int start, end, return_value
    return_value = iterator.next_record_by_offset(seg_view[j].file_offset)
    if return_value < 0:
        raise IOError(f"Failed to read record at offset {seg_view[j].file_offset}")

    setTrimRegion(&seg_view[j], &start, &end, flank_size)
    trim_segment(iterator.bam_record, dest_record, j, start, end)
    output_fa.write(dest_record)


cpdef output_reference_flank(Cluster[::1] clt_view, dict cluster_data_by_tid, tuple block, object cmd_args):
    """
    Output reference flank sequences for a range of clusters.

    Parameters:
    - clt_view: Array of Cluster objects.
    - cluster_data_by_tid: Dictionary mapping TID to cluster data.
    - block: Tuple (start_idx, end_idx) specifying the range of clusters to process.
    - cmd_args: Command-line arguments object containing configuration.

    This function resets the breakpoint of each cluster to the most common site
    within its segments and extracts the reference flank sequences.
    """
    # Reset breakpoint to most common site
    start_idx, end_idx = block
    cdef object counter = Counter()
    cdef Segment[::1] seg_view
    for i in range(start_idx, end_idx):
        seg_view = cluster_data_by_tid[clt_view[i].tid][1]
        for j in range(clt_view[i].start_idx, clt_view[i].end_idx):
            if overhang_too_short(&seg_view[j], cmd_args.overhang):
                continue
            counter[seg_view[j].ref_position] += 1
        
        clt_view[i].ref_start = counter.most_common(1)[0][0]
        clt_view[i].ref_end = clt_view[i].ref_start + 1
        counter.clear()
    
    # Extract reference flank sequences for the processed clusters
    cdef bytes ref_fn = cmd_args.ref_fn.encode('utf-8')
    extract_ref_flankseq(ref_fn, &clt_view[0], start_idx, end_idx)


cpdef merge_output():
    """
    Merge output files into a single result, with additional flag parsing.
    """
    clt_files = [os.path.join("tmp_anno", f) for f in os.listdir("tmp_anno") if f.endswith("cltFormated.txt")]
    anno_files = [os.path.join("tmp_anno", f) for f in os.listdir("tmp_anno") if f.endswith("annoFormated.txt")]

    # Load and merge clusters
    clt_dfs = [pd.read_csv(f, sep="\t", header=None) for f in clt_files]
    clt_df = pd.concat(clt_dfs, ignore_index=True)
    clt_df.columns = [
        "insertion_id", "chrom", "start", "end", "prob", "total_support_reads",
        "leftclip_reads", "spanning_reads", "rightclip_reads", "assembled",
        "tsd_seq", "insertion_seq", "left_flank_seq", "right_flank_seq", "flag"
        ]
    # Parse the flag field
    parse_flag(clt_df)

    # Load and merge annotations
    anno_dfs = [pd.read_csv(f, sep="\t", header=None) for f in anno_files]
    anno_df = pd.concat(anno_dfs, ignore_index=True)
    anno_df.columns = ["insertion_id", "strand", "te_family", "query_region", "annotation_region"]

    # Merge clt_df and anno_df
    result_df = pd.merge(clt_df, anno_df, on="insertion_id", how="outer")
    
    # Create the "extra_info" column
    result_df["extra_info"] = result_df.apply(generate_extra_info, axis=1)

    # Select necessary columns
    result_df = result_df[[
        "chrom", "start", "end", "te_family", "prob", "strand", "passed_filtering", "query_region", "annotation_region",
        "total_support_reads", "tsd_seq", "insertion_seq", "left_flank_seq", "right_flank_seq", "extra_info"
    ]]
    
    # Save the result to a file
    result_df.to_csv("result.tsv", sep="\t", index=False)


def parse_flag(df):
    df['passed_filtering'] = (df['flag'] & CLT_PASS) != 0
    df['assembled'] = (df['flag'] & CLT_ASSEMBLED) != 0
    df['has_polya'] = (df['flag'] & CLT_POLYA) != 0
    df['has_tsd'] = (df['flag'] & CLT_TSD) != 0
    df['singleton'] = (df['flag'] & CLT_SINGLE_TE) != 0
    df['self2self'] = (df['flag'] & CLT_SELF_TO_SELF) != 0
    df['solo_ltr'] = (df['flag'] & CLT_SOLO_LTR) != 0

    # Add reconstructed_ends column
    def define_reconstructed_ends(flag):
        if (flag & CLT_LEFT_FLANK_MAP) != 0:
            return "only_left"
        elif (flag & CLT_RIGHT_FLANK_MAP) != 0:
            return "only_right"
        elif (flag & (CLT_DIFF_FLANK_MAP | CLT_SAME_FLANK_MAP)) != 0:
            return "both_end"
        else:
            return "unknown"
    df['reconstructed_ends'] = df['flag'].apply(define_reconstructed_ends)

    # Add truncation column
    def define_truncation(flag):
        if ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_FULL) != 0):
            return "full"
        elif ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_UNKNOWN) != 0):
            return "3p_unknown"
        elif ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_UNKNOWN) == 0):
            return "3p_turncated"
        elif ((flag & CLT_3P_FULL) != 0) and ((flag & CLT_5P_UNKNOWN) != 0):
            return "5p_unknown"
        elif ((flag & CLT_3P_FULL) != 0) and ((flag & CLT_5P_UNKNOWN) == 0):
            return "5p_turncated"
        elif ((flag & CLT_3P_UNKNOWN) != 0) and ((flag & CLT_5P_UNKNOWN) != 0):
            return "5p3p_unknown"
        elif ((flag & CLT_5P_FULL) == 0) and ((flag & CLT_3P_FULL) == 0):
            return "5p3p_truncated"
        else:
            return "unknown"
    df['truncation'] = df['flag'].apply(define_truncation)

    # Add te_class column
    def define_te_class(flag):
        if (flag & CLT_DNA) != 0:
            return "DNA"
        elif (flag & CLT_LTR) != 0:
            return "LTR"
        elif (flag & CLT_LINE) != 0:
            return "LINE"
        elif (flag & CLT_SINE) != 0:
            return "SINE"
        elif (flag & CLT_RETROPOSON) != 0:
            return "Retroposon"
        else:
            return "unknown"
    df['te_class'] = df['flag'].apply(define_te_class)


def generate_extra_info(row):
    """
    Generate the extra_info string for a given row.
    """
    return (
        f"insID={row['insertion_id']},"
        f"leftClipReads={row['leftclip_reads']},"
        f"spanningReads={row['spanning_reads']},"
        f"rightClipReads={row['rightclip_reads']},"
        f"teClass={row['te_class']},"
        f"assembled={row['assembled']},"
        f"truncation={row['truncation']},"
        f"reconstructedEnds={row['reconstructed_ends']},"
        f"hasPolyA={row['has_polya']},"
        f"hasTSD={row['has_tsd']},"
        f"singleton={row['singleton']},"
        f"self2self={row['self2self']},"
        f"soloLTR={row['solo_ltr']}"
    )
