import os
import pandas as pd
import logging
import subprocess
from collections import Counter
from .genotyper import genotype_from_counts, define_genotype_threshold

logger = logging.getLogger(__name__)

########################
### Constants ###
########################
cdef int MAX_POSITION = (1 << 31) - 1


cdef inline bint _resize_if_needed(object arr, int size, int *capacity, int *threshold):
    if size < threshold[0]:
        return False
    capacity[0] = <int>(capacity[0] * 1.5)
    threshold[0] = <int>(capacity[0] * 0.9)
    arr.resize((capacity[0],), refcheck=False)
    return True


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
            logger.error(f"Error opening BAM file: {e}")
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
cdef output_seg_seqs(Segment[::1] seg_view, BamFile genome_bam, Args args):
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
            j = get_output_segidx(&clt_view[i], &seg_view[0], args)
            output_single_seq(seg_view, output_fa, iterator, dest_record, j)
            output_fa.close()
    finally:
        bam_destroy1(dest_record)
        del iterator
        genome_bam.close()
        del genome_bam


cpdef int output_read_as_assmbly(Cluster[::1] clt_view, dict cluster_data_by_tid, object cmd_args, int i, str output_fn):
    """
    Output one segment sequence as assembly for high-frequency cluster that failed to be assembled.
    """
    if clt_view[i].numSegRaw < 3:
        return 0

    cdef Segment[::1] seg_view = cluster_data_by_tid[clt_view[i].tid][1]
    cdef BamFile genome_bam = BamFile(cmd_args.genome_bam_fn, "rb", cmd_args.num_thread)
    cdef Iterator iterator = Iterator(genome_bam, clt_view[i].tid)
    cdef bam1_t *dest_record = bam_init1()
    cdef BamFile output_fa = BamFile(output_fn, "wF", cmd_args.num_thread, genome_bam)
    cdef Args args

    try:
        args.overhang = cmd_args.overhang
        j = get_output_segidx(&clt_view[i], &seg_view[0], args)
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


cpdef merge_output(genotyper='bayesian', output_format='both', sample_name='SAMPLE'):
    """
    Merge output files into a single result, with additional flag parsing.
    """
    clt_files = [os.path.join("tmp_anno", f) for f in os.listdir("tmp_anno") if f.endswith("cltFormated.txt")]
    anno_files = [os.path.join("tmp_anno", f) for f in os.listdir("tmp_anno") if f.endswith("annoFormated.txt")]

    # Load and merge clusters
    clt_dfs = [pd.read_csv(f, sep="\t", header=None) for f in clt_files]
    clt_df = pd.concat(clt_dfs, ignore_index=True)
    clt_df.columns = [
        "insertion_id", "chrom", "start", "end", "prob", "total_support",
        "leftclip_reads", "spanning_reads", "rightclip_reads", "num_ref", "assembled",
        "tsd_seq", "insertion_seq", "upstream_seq", "downstream_seq", "flag", "frequency"
        ]
    # Parse the flag field
    parse_flag(clt_df, genotyper=genotyper)

    # Load and merge annotations
    anno_dfs = [pd.read_csv(f, sep="\t", header=None) for f in anno_files]
    anno_df = pd.concat(anno_dfs, ignore_index=True)
    anno_df.columns = ["insertion_id", "strand", "family", "query_region", "target_region"]

    # Merge clt_df and anno_df
    result_df = pd.merge(clt_df, anno_df, on="insertion_id", how="outer")

    # Create the "extra_info" column
    result_df["extra_info"] = result_df.apply(generate_extra_info, axis=1)

    # Write VCF output (needs full DataFrame before column selection)
    if output_format in ('vcf', 'both'):
        vcf_path = os.path.abspath("result.vcf")
        write_vcf(result_df, vcf_path, sample_name)
        logger.info("VCF output: %s", vcf_path)

    # Select necessary columns for TSV
    output_columns = [
        "chrom", "start", "end", "family", "frequency", "strand", "genotype", "passed", "query_region",
        "target_region", "total_support", "tsd_seq", "insertion_seq", "upstream_seq", "downstream_seq", "extra_info"
    ]
    if 'genotype_quality' in result_df.columns:
        # Insert genotype_quality after genotype
        idx = output_columns.index("genotype")
        output_columns.insert(idx + 1, "genotype_quality")

    # Save the result as TSV
    if output_format in ('tsv', 'both'):
        tsv_path = os.path.abspath("result.tsv")
        result_df[output_columns].to_csv(tsv_path, sep="\t", index=False)
        logger.info("TSV output: %s", tsv_path)

    # Summary logging
    num_passed = result_df['passed'].sum()
    num_failed = len(result_df) - num_passed
    logger.info("Total insertions: %d (passed: %d, failed: %d)", len(result_df), num_passed, num_failed)


def _define_reconstructed_ends(flag):
    if (flag & CLT_LEFT_FLANK_MAP) != 0:
        return "only_left"
    elif (flag & CLT_RIGHT_FLANK_MAP) != 0:
        return "only_right"
    elif (flag & (CLT_DIFF_FLANK_MAP | CLT_SAME_FLANK_MAP)) != 0:
        return "both_end"
    else:
        return "unknown"


def _define_truncation(flag):
    if ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_FULL) != 0):
        return "full"
    elif ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_UNKNOWN) != 0):
        return "3p_unknown"
    elif ((flag & CLT_5P_FULL) != 0) and ((flag & CLT_3P_UNKNOWN) == 0):
        return "3p_truncated"
    elif ((flag & CLT_3P_FULL) != 0) and ((flag & CLT_5P_UNKNOWN) != 0):
        return "5p_unknown"
    elif ((flag & CLT_3P_FULL) != 0) and ((flag & CLT_5P_UNKNOWN) == 0):
        return "5p_truncated"
    elif ((flag & CLT_3P_UNKNOWN) != 0) and ((flag & CLT_5P_UNKNOWN) != 0):
        return "5p3p_unknown"
    elif ((flag & CLT_5P_FULL) == 0) and ((flag & CLT_3P_FULL) == 0):
        return "5p3p_truncated"
    else:
        return "unknown"


def _define_te_class(flag):
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


def parse_flag(df, genotyper='bayesian'):
    df['passed'] = (df['flag'] & CLT_PASS) != 0
    df['assembled'] = (df['flag'] & CLT_ASSEMBLED) != 0
    df['has_polya'] = (df['flag'] & CLT_POLYA) != 0
    df['has_tsd'] = (df['flag'] & CLT_TSD) != 0
    df['singleton'] = (df['flag'] & CLT_SINGLE_TE) != 0
    df['self2self'] = (df['flag'] & CLT_SELF_TO_SELF) != 0
    df['solo_ltr'] = (df['flag'] & CLT_SOLO_LTR) != 0

    df['reconstructed_ends'] = df['flag'].apply(_define_reconstructed_ends)
    df['truncation'] = df['flag'].apply(_define_truncation)
    df['te_class'] = df['flag'].apply(_define_te_class)

    if genotyper == 'threshold':
        df['genotype'] = df['frequency'].apply(define_genotype_threshold)
    else:
        genotypes = [genotype_from_counts(
            r.leftclip_reads, r.spanning_reads, r.rightclip_reads, r.num_ref
        ) for _, r in df.iterrows()]
        df['genotype'] = [g[0] for g in genotypes]
        df['genotype_quality'] = [g[1] for g in genotypes]


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


def _vcf_escape(value):
    """Escape characters reserved in VCF INFO field values."""
    if not isinstance(value, str):
        value = str(value) if value is not None else "."
    return (value.replace(";", "%3B")
                .replace("=", "%3D")
                .replace(",", "%2C")
                .replace("\n", " ")
                .replace("\t", " "))


def write_vcf(df, output_path, sample_name="SAMPLE"):
    """
    Write insertion results in VCF 4.3 format.

    Parameters
    ----------
    df : pandas.DataFrame
        Result DataFrame with all columns from merge_output.
    output_path : str
        Path for the output .vcf file.
    sample_name : str
        Name for the single sample column.
    """
    has_gq = "genotype_quality" in df.columns
    chromosomes = sorted(df["chrom"].unique())

    with open(output_path, "w") as f:
        # --- Header lines ---
        f.write("##fileformat=VCFv4.3\n")
        f.write("##source=LOCATE\n")

        for chrom in chromosomes:
            if pd.isna(chrom):
                continue
            f.write(f"##contig=<ID={chrom}>\n")

        # INFO definitions
        info_defs = {
            "END": ("Integer", "End position of the structural variant"),
            "SVTYPE": ("String", "Type of structural variant"),
            "SVLEN": ("Integer", "Insertion length"),
            "FAMILY": ("String", "TE family name(s)"),
            "AF": ("Float", "Allele frequency"),
            "STRAND": ("String", "Insertion strand orientation (+/-)"),
            "TSD": ("String", "Target site duplication sequence"),
            "TE_CLASS": ("String", "TE class (DNA/LTR/LINE/SINE/Retroposon/unknown)"),
            "TRUNCATION": ("String", "Truncation status"),
            "RECONSTRUCTED_ENDS": ("String", "Reconstructed ends status"),
            "HAS_POLYA": ("Integer", "Has polyA tail"),
            "HAS_TSD": ("Integer", "Has target site duplication"),
            "ASSEMBLED": ("Integer", "Insertion was assembled"),
            "SINGLETON": ("Integer", "Singleton insertion"),
            "SELF2SELF": ("Integer", "Self-to-self insertion"),
            "SOLO_LTR": ("Integer", "Solo LTR"),
            "SUPPORT": ("Integer", "Total supporting reads"),
            "LEFT_CLIP": ("Integer", "Left-clipped reads"),
            "SPANNING": ("Integer", "Spanning/mid-insert reads"),
            "RIGHT_CLIP": ("Integer", "Right-clipped reads"),
            "QV": ("Float", "ML model probability"),
        }
        for tag, (vtype, desc) in info_defs.items():
            f.write(f'##INFO=<ID={tag},Number=1,Type={vtype},Description="{desc}">\n')

        # FORMAT definitions
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')

        # FILTER definitions
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=FAIL,Description="Failed post-filtering">\n')

        # Column header
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        # --- Data lines ---
        for _, row in df.iterrows():
            chrom = str(row["chrom"]) if not pd.isna(row.get("chrom")) else "."
            pos = int(row["start"]) + 1 if not pd.isna(row.get("start")) else 1
            vid = _vcf_escape(str(row["insertion_id"])) if not pd.isna(row.get("insertion_id")) else "."
            qual = str(int(row["genotype_quality"])) if has_gq and not pd.isna(row.get("genotype_quality")) else "."

            # FILTER
            passed = row.get("passed", False)
            filt = "PASS" if passed else "FAIL"

            # Build INFO
            svlen = len(str(row.get("insertion_seq", "")) or "") if not pd.isna(row.get("insertion_seq")) else 0
            end = int(row["end"]) if not pd.isna(row.get("end")) else pos

            info_parts = [
                f"END={end}",
                "SVTYPE=INS",
                f"SVLEN={svlen}",
                f"FAMILY={_vcf_escape(row.get('family', '.'))}",
                f"AF={float(row.get('frequency', 0))}",
                f"STRAND={_vcf_escape(row.get('strand', '.'))}",
                f"TSD={_vcf_escape(row.get('tsd_seq', '.'))}",
                f"TE_CLASS={_vcf_escape(row.get('te_class', '.'))}",
                f"TRUNCATION={_vcf_escape(row.get('truncation', '.'))}",
                f"RECONSTRUCTED_ENDS={_vcf_escape(row.get('reconstructed_ends', '.'))}",
                f"HAS_POLYA={1 if row.get('has_polya') else 0}",
                f"HAS_TSD={1 if row.get('has_tsd') else 0}",
                f"ASSEMBLED={1 if row.get('assembled') else 0}",
                f"SINGLETON={1 if row.get('singleton') else 0}",
                f"SELF2SELF={1 if row.get('self2self') else 0}",
                f"SOLO_LTR={1 if row.get('solo_ltr') else 0}",
                f"SUPPORT={int(row.get('total_support', 0))}",
                f"LEFT_CLIP={int(row.get('leftclip_reads', 0))}",
                f"SPANNING={int(row.get('spanning_reads', 0))}",
                f"RIGHT_CLIP={int(row.get('rightclip_reads', 0))}",
                f"QV={float(row.get('prob', 0))}",
            ]
            info = ";".join(info_parts)

            # FORMAT and sample
            genotype = row.get("genotype", "./.")
            gt_val = str(genotype) if not pd.isna(genotype) else "./."
            if has_gq and not pd.isna(row.get("genotype_quality")):
                fmt = "GT:GQ"
                sample_val = f"{gt_val}:{int(row['genotype_quality'])}"
            else:
                fmt = "GT"
                sample_val = gt_val

            f.write(f"{chrom}\t{pos}\t{vid}\tN\t<INS>\t{qual}\t{filt}\t{info}\t{fmt}\t{sample_val}\n")
