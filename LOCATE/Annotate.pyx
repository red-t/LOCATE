import os
import numpy as np
import subprocess

##################
### Data Types ###
##################
AnnoDt = np.dtype([
    ('idx',         np.int32),
    ('cltTid',      np.int32),
    ('cltIdx',      np.int32),
    ('query_start', np.int32),
    ('query_end',   np.int32),
    ('strand',      np.uint8),
    ('tid',         np.int32),
    ('ref_start',   np.int32),
    ('ref_end',     np.int32),
    ('flag',        np.uint32),
    ('extra',       np.int32),
])


#########################
### Annotate Assembly ###
#########################
cdef annotate_assembly(Cluster[::1] clt_view, int start_idx, int end_idx, object cmd_args):
    cdef int i
    for i in range(start_idx, end_idx):
        map_flank_to_assembly(clt_view[i].tid, clt_view[i].idx)
        extract_ins_seq(&clt_view[i])
        map_assm_flank_to_local(clt_view[i].tid, clt_view[i].idx)
        if os.path.isfile("tmp_anno/{}_{}_AssmFlankToLocal.bam".format(clt_view[i].tid, clt_view[i].idx)) == True:
            re_extract_ins_seq(&clt_view[i])
        map_ins_seq_to_te(clt_view[i].tid, clt_view[i].idx, cmd_args)
        

cdef map_flank_to_assembly(int tid, int idx):
    cdef str target_fn = "tmp_assm/{}_{}_assembled.fa".format(tid, idx)
    cdef str query_fn = "tmp_anno/{}_{}_flank.fa".format(tid, idx)
    cdef str output_fn = "tmp_anno/{}_{}_FlankToAssm.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(target_fn, query_fn, output_fn)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef map_assm_flank_to_local(int tid, int idx):
    cdef str target_fn = "tmp_anno/{}_{}_local.fa".format(tid, idx)
    cdef str query_fn = "tmp_anno/{}_{}_assmFlank.fa".format(tid, idx)
    cdef str output_fn = "tmp_anno/{}_{}_AssmFlankToLocal.bam".format(tid, idx)
    cdef str cmd = "minimap2 -x sr --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(target_fn, query_fn, output_fn)
    if os.path.isfile(query_fn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


cdef map_ins_seq_to_te(int tid, int idx, object cmd_args):
    cdef str query_fn = "tmp_anno/{}_{}_insertion.fa".format(tid, idx)
    cdef str output_fn = "tmp_anno/{}_{}_InsToTE.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(cmd_args.te_fn, query_fn, output_fn)
    if os.path.isfile(query_fn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


##########################
### Annotate Insertion ###
##########################
cdef object get_class_arr(object cmd_args):
    cdef int size = 0, capacity = 200, threshold = int(0.9 * capacity)
    cdef object class_arr = np.zeros(capacity, dtype=np.uint32)
    cdef uint32_t[::1] class_view = class_arr

    for l in open(cmd_args.class_fn, "r"):
        l = l.strip().split()
        if len(l) < 2:
            continue
        
        if size >= threshold:
            capacity = int(1.5 * capacity)
            threshold = int(0.9 * capacity)
            class_arr.resize((capacity,), refcheck=False)
            class_view = class_arr
        
        if l[1] == "DNA":
            class_view[size] = CLT_DNA
        elif l[1] == "LTR":
            class_view[size] = CLT_LTR
        elif l[1] == "LINE":
            class_view[size] = CLT_LINE
        elif l[1] == "SINE":
            class_view[size] = CLT_SINE
        elif l[1] == "RETROPOSON":
            class_view[size] = CLT_RETROPOSON
        
        size += 1
    
    class_arr.resize((size,), refcheck=False)
    return class_arr


cdef object get_size_arr(object cmd_args):
    cdef int size = 0, capacity = 200, threshold = int(0.9 * capacity)
    cdef object size_arr = np.zeros(capacity, dtype=np.int32)
    cdef int32_t[::1] size_view = size_arr
    cdef indexFn = cmd_args.te_fn + ".fai"

    for l in open(indexFn, "r"):
        l = l.strip().split()
        if len(l) < 1:
            continue
        
        if size >= threshold:
            capacity = int(1.5 * capacity)
            threshold = int(0.9 * capacity)
            size_arr.resize((capacity,), refcheck=False)
            size_view = size_arr
        
        size_view[size] = int(l[1])
        size += 1
    
    size_arr.resize((size,), refcheck=False)
    return size_arr


cdef object get_ltr_arr():
    cdef int num_te = 0, max_num = 190
    cdef object ltr_arr = np.zeros(200, dtype=np.int32)
    cdef int32_t[::1] ltr_view = ltr_arr

    for l in open("tmp_anno/ltrSize.txt", "r"):
        l = l.strip().split()
        if len(l) < 1:
            continue
        
        if num_te >= max_num:
            max_num = ltr_view.shape[0] + 200
            ltr_arr.resize((max_num,), refcheck=False)
            ltr_view = ltr_arr
            max_num -= 10
        
        ltr_view[num_te] = int(l[0])
        num_te += 1
    
    ltr_arr.resize((num_te,), refcheck=False)
    return ltr_arr


cdef object annotate_ins_seq(Cluster[::1] clt_view, int start_idx, int end_idx, object cmd_args):
    cdef int i, num_anno, size = 0, capacity = 3000, threshold = int(0.9 * capacity)
    cdef object anno_arr = np.zeros(capacity, dtype=AnnoDt)
    cdef object class_arr = get_class_arr(cmd_args)
    cdef object size_arr = get_size_arr(cmd_args)
    cdef object ltr_arr = get_ltr_arr()
    cdef Annotation[::1] anno_view = anno_arr
    cdef uint32_t[::1] class_view = class_arr
    cdef int[::1] size_view = size_arr
    cdef int[::1] ltr_view = ltr_arr
    cdef str bam_fn, assembly_fn

    for i in range(start_idx, end_idx):
        # 1. Check if cluster has assembled insSeq
        assembly_fn = "tmp_assm/{}_{}_polished.fa".format(clt_view[i].tid, clt_view[i].idx)
        if os.path.isfile(assembly_fn) != 0:
                clt_view[i].flag |= CLT_ASSEMBLED
        
        bam_fn = "tmp_anno/{}_{}_InsToTE.bam".format(clt_view[i].tid, clt_view[i].idx)
        if os.path.isfile(bam_fn) == False:
            continue

        if size >= threshold:
            capacity = int(1.5 * capacity)
            threshold = int(0.9 * capacity)
            anno_arr.resize((capacity,), refcheck=False)
            anno_view = anno_arr

        # 2. Annotate TE fragment, polyA/T for insSeq
        num_anno = fill_anno_arr(&clt_view[i], &anno_view[size], &class_view[0], i)

        # 3. Annotate TSD
        map_tsd_to_local(clt_view[i].tid, clt_view[i].idx)
        bam_fn = "tmp_anno/{}_{}_TsdToLocal.bam".format(clt_view[i].tid, clt_view[i].idx)
        if os.path.isfile(bam_fn) == True:
            annotate_tsd(&clt_view[i], &anno_view[size], num_anno)

        # 4. Define insSeq structure
        set_ins_structure(&clt_view[i], &anno_view[size], num_anno, &class_view[0], &size_view[0], &ltr_view[0])

        # 5. Perform post-filtering
        post_filter(&clt_view[i])
        
        size += num_anno
    
    anno_arr.resize((size,), refcheck=False)
    anno_arr.sort(order=['idx', 'query_start', 'query_end'])
    return anno_arr


cdef map_tsd_to_local(int tid, int idx):
    cdef str target_fn = "tmp_anno/{}_{}_local.fa".format(tid, idx)
    cdef str query_fn = "tmp_anno/{}_{}_tsd.fa".format(tid, idx)
    cdef str output_fn = "tmp_anno/{}_{}_TsdToLocal.bam".format(tid, idx)
    cdef str cmd = "minimap2 -k11 -w5 --sr -O4,8 -n2 -m20 --secondary=no -t 1 -aY {} {} | " \
                   "samtools view -bhS -o {} -".format(target_fn, query_fn, output_fn)
    if os.path.isfile(query_fn) == False:
        return
    
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


##########################
### Genotype Insertion ###
##########################
cdef object compute_frequency(Cluster[::1] clt_view, tuple block, dict cluster_data_by_tid, object cmd_args, int extra_thread):
    cdef BamFile genome_bam = BamFile(cmd_args.genome_bam_fn, "rb", 1 + extra_thread)
    cdef Iterator iterator
    cdef Segment[::1] seg_view

    # Step 1: Process each clusters
    cdef i, j, beg, end, num_ref, return_value, start_idx, end_idx
    start_idx, end_idx = block
    for i in range(start_idx, end_idx):
        # Collect offsets of each support reads
        seg_view = cluster_data_by_tid[clt_view[i].tid][1]
        alt_offsets = set()
        for j in range(clt_view[i].start_idx, clt_view[i].end_idx):
            if overhang_too_short(&seg_view[j], cmd_args.overhang):
                continue
            alt_offsets.add(seg_view[j].file_offset)

        # Step 2: Initialize iterator for the corresponding region
        iterator = Iterator(genome_bam, clt_view[i].tid, clt_view[i].ref_start - 50000, clt_view[i].ref_start + 50000)

        # Step 3: Count ref_support reads
        try:
            num_ref = 0
            beg = clt_view[i].ref_start - 75
            end = clt_view[i].ref_start + 75
            while True:
                return_value = iterator.next_record_by_tid()
                if return_value < 0:
                    break
                if bam_is_invalid(iterator.bam_record):
                    continue
                if iterator.offset in alt_offsets:
                    continue
                elif (iterator.bam_record.core.pos < beg) and (bam_endpos(iterator.bam_record) > end):
                    num_ref += 1
        finally:
            del iterator

        # Step 4: Compute frequency
        clt_view[i].frequency = float(clt_view[i].numLeft + clt_view[i].numRight + 2*clt_view[i].numMiddle) / (clt_view[i].numLeft + clt_view[i].numRight + 2*clt_view[i].numMiddle + 2*num_ref)
        
    genome_bam.close()




###################################
### Annotate Insertion Sequence ###
###################################
cpdef annotate_cluster(Cluster[::1] clt_view, tuple block, dict cluster_data_by_tid, object cmd_args, int extra_thread):
    # Step 1: Define insertion seq from assembled contig(s)
    start_idx, end_idx = block
    annotate_assembly(clt_view, start_idx, end_idx, cmd_args)

    # Step 2: Annotate TE-fragment, PolyA/T, TSD and structure for insSeq
    #         Also perform post-filtering
    cdef object anno_arr = annotate_ins_seq(clt_view, start_idx, end_idx, cmd_args)

    # Step 3: Genotype clusters
    compute_frequency(clt_view, block, cluster_data_by_tid, cmd_args, extra_thread)
    
    # Step 4: Output formated cluster and annotation records
    cdef Annotation[::1] anno_view = anno_arr
    cdef bytes te_fn = cmd_args.te_fn.encode()
    cdef bytes ref_fn = cmd_args.ref_fn.encode()
    if (anno_view.shape[0] > 0):
        output_annotations(&anno_view[0], anno_view.shape[0], start_idx, te_fn)
    
    if (clt_view.shape[0] > 0):
        output_clusters(&clt_view[0], start_idx, end_idx, ref_fn, te_fn)
    