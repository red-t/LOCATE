import os
import numpy as np
import pandas as pd
import subprocess
from autogluon.tabular import TabularPredictor

##################
### Data Types ###
##################
SegmentDt = np.dtype([
    ('map_qual',            np.uint8),
    ('query_start',         np.int32),
    ('query_end',           np.int32),
    ('ref_position',        np.int32),
    ('seg_type',            np.uint8),
    ('aln_type',            np.uint8),
    ('file_offset',         np.int64),
    ('aln_refstart',        np.int32),
    ('aln_refend',          np.int32),
    ('order',               np.uint8),
    ('numSeg',              np.uint8),
    ('overhang',            np.int32),
    ('match_len',           np.int32),
    ('read_len',            np.int32),
    ('aln_location_type',   np.uint8),
    ('num_te_aln',          np.uint8),
    ('sum_query_maplen',    np.int32),
    ('sum_aln_score',       np.float32),
    ('sum_divergence',      np.float32),
    ('start_idx',           np.int32),
    ('end_idx',             np.int32),
])

TeAlignmentDt = np.dtype([
    ('seg_idx',      np.int32),
    ('aln_score',    np.int32),
    ('query_start',  np.int32),
    ('query_end',    np.int32),
    ('map_len',      np.int32),
    ('divergence',   np.float32),
])

ClusterDt = np.dtype([
    ('tid',                 np.int32),
    ('idx',                 np.int32),
    ('ref_start',           np.int32),
    ('ref_end',             np.int32),
    ('start_idx',           np.int32),
    ('end_idx',             np.int32),
    ('numSeg',              np.float32),
    ('cltType',             np.uint8),
    ('locationType',        np.uint8),
    ('numSegType',          np.uint8),
    ('entropy',             np.float32),
    ('balanceRatio',        np.float32),
    ('lowMapQualFrac',      np.float32),
    ('dualClipFrac',        np.float32),
    ('alnFrac1',            np.float32),
    ('alnFrac2',            np.float32),
    ('alnFrac4',            np.float32),
    ('alnFrac8',            np.float32),
    ('alnFrac16',           np.float32),
    ('meanMapQual',         np.float32),
    ('meanAlnScore',        np.float32),
    ('meanQueryMapFrac',    np.float32),
    ('meanDivergence',      np.float32),
    ('bgDiv',               np.float32),
    ('bgDepth',             np.float32),
    ('bgReadLen',           np.float32),
    ('teAlignedFrac',       np.float32),
    ('isInBlacklist',       np.uint8),
    ('probability',         np.float32),
    ('flag',                np.uint32),
    ('numSegRaw',           np.int32),
    ('numLeft',             np.int32),
    ('numMiddle',           np.int32),
    ('numRight',            np.int32),
    ('tid1',                np.int32),
    ('leftMost',            np.int32),
    ('tid2',                np.int32),
    ('rightMost',           np.int32),
    ('insLen',              np.int32),
    ('repTid',              np.int32),
])


########################
### Construct SegArr ###
########################
cdef object get_seg_arr(BamFile genome_bam, Args args):
    """
    Constructs an array of segments (seg_arr) from a BAM file.

    Parameters:
        genome_bam (BamFile): Input BAM file object.
        args (Args): Object containing parameters such as tid and minimum segment length.

    Returns:
        object: A NumPy array containing segment information with dtype SegmentDt.
    """
    cdef int return_value, size = 0, capacity = 10000, threshold = int(capacity * 0.9)
    cdef object seg_arr = np.zeros(capacity, dtype=SegmentDt)
    cdef Segment[::1] seg_view = seg_arr
    cdef Iterator iterator = Iterator(genome_bam, args.tid)
    
    try:
        while True:
            return_value = iterator.next_record_by_tid()
            if return_value < 0:
                seg_arr.resize((size,), refcheck=False)
                return seg_arr

            if bam_is_invalid(iterator.bam_record):
                continue

            if size >= threshold:
                capacity = int(1.5 * capacity)
                threshold = int(capacity * 0.9)
                seg_arr.resize((capacity,), refcheck=False)
                seg_view = seg_arr

            return_value = fill_seg_arr(iterator.bam_record, &seg_view[size], iterator.offset, args.min_seg_len)
            size += return_value
    finally:
        del iterator


cdef update_seg_arr(Segment[::1] seg_view, Args args):
    """
    Updates each segment in the segment array with additional information.

    Parameters:
        seg_view (Segment[::1]): Memory view of the segment array.
        args (Args): Object containing parameters such as repeat and gap region lists.
    """
    cdef int i
    for i in range(seg_view.shape[0]):
        update_segment(&seg_view[i], args.repeat_ailist, args.gap_ailist)


cdef update_seg_arr_by_te(Segment[::1] seg_view, Args args):
    """
    Updates the segment array based on transposable element (TE) alignment information.

    Parameters:
        seg_view (Segment[::1]): Memory view of the segment array.
        args (Args): Object containing parameters such as thread count and tid.
    """
    cdef BamFile te_bam = BamFile("tmp_build/all_seg_{}.bam".format(args.tid), "rb", args.num_thread)
    cdef Iterator iterator = Iterator(te_bam)
    cdef object te_arr = get_te_arr(iterator)
    cdef TeAlignment[::1] te_view = te_arr
    cdef int i

    try:
        te_arr.sort(order=['seg_idx', 'query_start'])
        for i in range(te_view.shape[0]):
            update_seg_by_te_arr(&seg_view[0], &te_view[0], i)
    finally:
        del iterator
        te_bam.close()
        del te_bam


#########################
### Construct TeArr ###
#########################
cdef map_seg_to_te(str te_fn, Args args):
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t {0} -aY {1} tmp_build/all_seg_{2}.fa | " \
                   "samtools view -@ {0} -bhS -o tmp_build/all_seg_{2}.bam -".format(args.num_thread, te_fn, args.tid)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef object get_te_arr(Iterator iterator):
    """
    Constructs an array of transposon alignments (te_arr) from a BAM file.

    Parameters:
        iterator (Iterator): An iterator linked to the transposon BAM file.

    Returns:
        object: A NumPy array containing transposon alignments information with dtype TeAlignmentDt.
    """
    cdef int return_value, size = 0, capacity = 10000, threshold = int(capacity * 0.9)
    cdef object te_arr = np.zeros(capacity, dtype=TeAlignmentDt)
    cdef TeAlignment[::1] te_view = te_arr

    try:
        while True:
            return_value = iterator.next_record()
            if return_value < 0:
                break

            if bam_is_invalid(iterator.bam_record):
                continue

            if size >= threshold:
                capacity = int(capacity * 1.5)
                threshold = int(capacity * 0.9)
                te_arr.resize((capacity,), refcheck=False)
                te_view = te_arr

            fill_te_arr(iterator.bam_record, &te_view[size])
            size += 1
    finally:
        te_arr.resize((size,), refcheck=False)  # 调整到实际大小
        return te_arr


##########################
### Construct CltArr ###
##########################
cdef object get_clt_arr(Segment[::1] seg_view, Args args):
    """
    Constructs an array of clusters (clt_arr) from the given segment array (seg_view).

    Parameters:
        seg_view (Segment[::1]): Memory view of the segment array.
        args (Args): Object containing parameters such as overhang threshold and max distance.

    Returns:
        object: A NumPy array containing cluster information with dtype ClusterDt.
    """
    cdef int start_idx = 0, end_idx, size = 0, capacity = 10000, threshold = int(capacity * 0.9)
    cdef object clt_arr = np.zeros(capacity, dtype=ClusterDt)
    cdef Cluster[::1] clt_view = clt_arr

    while start_idx < seg_view.shape[0]:
        if overhang_too_short(&seg_view[start_idx], args.overhang):
            start_idx += 1
            continue

        if size >= threshold:
            capacity = int(capacity * 1.5)
            threshold = int(capacity * 0.9)
            clt_arr.resize((capacity,), refcheck=False)
            clt_view = clt_arr

        # Initialize current cluster
        clt_view[size].tid = args.tid
        clt_view[size].ref_start = seg_view[start_idx].ref_position - 1
        clt_view[size].ref_end = seg_view[start_idx].ref_position + args.max_dist
        clt_view[size].idx = size
        clt_view[size].start_idx = start_idx

        end_idx = start_idx + 1
        while end_idx < seg_view.shape[0]:
            if overhang_too_short(&seg_view[end_idx], args.overhang):
                end_idx += 1
                continue
            if seg_view[end_idx].ref_position > clt_view[size].ref_end:
                break

            clt_view[size].ref_end = seg_view[end_idx].ref_position + args.max_dist
            end_idx += 1

        clt_view[size].end_idx = end_idx
        clt_view[size].ref_end -= args.max_dist
        start_idx = end_idx
        size += 1

    clt_arr.resize((size,), refcheck=False)
    return clt_arr


cdef update_clt_arr(Cluster[::1] clt_view, Segment[::1] seg_view, BamFile genome_bam, Args args):
    """
    Updates the features of each cluster in the cluster array (clt_view).

    Parameters:
        clt_view (Cluster[::1]): Memory view of the cluster array.
        seg_view (Segment[::1]): Memory view of the segment array.
        genome_bam (BamFile): BAM file object containing genomic alignments.
        args (Args): Object containing parameters and auxiliary data for updating clusters.
    """
    cdef int i

    try:
        args.genome_bam = genome_bam.hts_file
        args.first_bam_record = bam_init1()
        args.second_bam_record = bam_init1()
        
        # Update cluster features
        for i in range(clt_view.shape[0]):
            update_cluster(&clt_view[i], &seg_view[0], args)
    finally:
        bam_destroy1(args.first_bam_record)
        bam_destroy1(args.second_bam_record)
        destroyAiList(args.repeat_ailist)
        destroyAiList(args.gap_ailist)


######################
### Filter Cluster ###
######################
cdef filter_by_blacklist(Cluster[::1] clt_view, Args args):
    cdef int i
    try:
        for i in range(clt_view.shape[0]):
            intersect_black_list(&clt_view[i], args)
    finally:
        destroyAiList(args.blacklist_ailist)


cdef object filter_by_model(object clt_arr, object cmd_args):
    cdef object clt_df = pd.DataFrame(clt_arr)
    try:
        filter_high_freq_clusters(clt_df, cmd_args.high_freq_model)
        filter_low_freq_clusters(clt_df, cmd_args.low_freq_model)
    finally:
        return clt_df.to_records(index=False)


cdef filter_high_freq_clusters(object clt_df, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object high_freq_df = clt_df.loc[(clt_df['cltType']==0) & (clt_df['teAlignedFrac']>=0.8) & (clt_df['isInBlacklist']==0)]

    probability = predictor.predict_proba(high_freq_df)
    probability.columns = ['0', 'probability']
    clt_df.update(probability)


cdef filter_low_freq_clusters(object clt_df, str modelPath):
    cdef object predictor = TabularPredictor.load(modelPath)
    cdef object low_freq_df = clt_df.loc[(clt_df['cltType']>0) & (clt_df['teAlignedFrac']>=0.8) & (clt_df['isInBlacklist']==0)]

    probability = predictor.predict_proba(low_freq_df)
    probability.columns = ['0', 'probability']
    clt_df.update(probability)


#####################
### Build Cluster ###
#####################
cpdef dict build_cluster(float bg_div, float bg_depth, float bg_read_len, object cmd_args, int tid, int extra_thread):
    """
    Builds clusters from genomic data.

    Parameters:
        bg_div (float): Background divergence value.
        bg_depth (float): Background depth value.
        bg_read_len (float): Background read length.
        cmd_args (object): Command-line arguments or configuration object containing file paths and parameters.
        tid (int): Target ID (chromosome or contig ID).
        extra_thread (int): Number of additional threads to use.

    Returns:
        dict: A dictionary where the key is the target ID (tid) and the value is a tuple containing:
              - clt_arr: Array of clusters.
              - seg_arr: Array of segments.
    """

    # Step 1: Construct segments
    cmd_args.num_thread = 1 + extra_thread
    cdef Args args = new_args(tid, bg_div, bg_depth, bg_read_len, cmd_args)
    cdef BamFile genome_bam = BamFile(cmd_args.genome_bam_fn, "rb", cmd_args.num_thread)
    cdef object seg_arr = get_seg_arr(genome_bam, args)

    # Step 2: Compute segment features
    cdef const char *chrom = sam_hdr_tid2name(genome_bam.header, tid)
    args.repeat_ailist = new_ailist(cmd_args.repeat_fn, chrom)
    args.gap_ailist = new_ailist(cmd_args.gap_fn, chrom)

    update_seg_arr(seg_arr, args)
    seg_arr.sort(order='ref_position')
    ouput_seg_seqs(seg_arr, genome_bam, args)

    map_seg_to_te(cmd_args.te_fn, args)
    update_seg_arr_by_te(seg_arr, args)
    
    # Step 3: Construct clusters
    cdef cluster_data = {}
    cdef object clt_arr = get_clt_arr(seg_arr, args)

    # Step 4: Compute cluster features
    update_clt_arr(clt_arr, seg_arr, genome_bam, args)

    # Step 5: Filter clusters
    args.blacklist_ailist = new_ailist(cmd_args.blacklist_fn, chrom)
    filter_by_blacklist(clt_arr, args)
    clt_arr = filter_by_model(clt_arr, cmd_args)

    # Step 6: Output cluster sequences
    output_highfreq_clusters_seqs(clt_arr, seg_arr, genome_bam, args)

    cluster_data[tid] = (clt_arr, seg_arr)
    
    genome_bam.close()
    del genome_bam
    return cluster_data


##############################
### Get High-Qual Clusters ###
##############################
cdef object get_high_quality_clusters(dict cluster_data_by_tid):
    """
    Extracts high-quality clusters from the given cluster data.

    Parameters:
        cluster_data_by_tid (dict): A dictionary where the key is the target ID (tid) and the value is a tuple containing:
                                    - clt_arr: Array of clusters.
                                    - seg_arr: Array of segments.

    Returns:
        object: A NumPy array containing high-quality clusters with dtype ClusterDt.
    """
    cdef int i, tid, size = 0, capacity = 2000, threshold = int(capacity * 0.9)
    cdef object high_quality_arr = np.zeros(capacity, dtype=ClusterDt)
    cdef Cluster[::1] high_quality_view = high_quality_arr
    cdef Cluster[::1] clt_view

    for tid in range(len(cluster_data_by_tid)):
        clt_view = cluster_data_by_tid[tid][0]
        if clt_view.shape[0] == 0:
            continue

        for i in range(clt_view.shape[0]):
            if is_lowqual_clt(&clt_view[i]):
                continue

            # Resize array if needed
            if size >= threshold:
                capacity = int(capacity * 1.5)
                threshold = int(capacity * 0.9)
                high_quality_arr.resize((capacity,), refcheck=False)
                high_quality_view = high_quality_arr

            # Add high-quality cluster to the array
            high_quality_view[size] = clt_view[i]
            size += 1

    # Resize to the final size
    high_quality_arr.resize((size,), refcheck=False)
    return high_quality_arr