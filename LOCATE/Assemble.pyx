import os
import numpy as np
import subprocess
from collections import OrderedDict, Counter, defaultdict
from cpython cimport PyBytes_FromStringAndSize

######################
### Local Assembly ###
######################
cdef int get_min_edge(int num_seg_raw):
    # Dynamically set min_edge
    if num_seg_raw < 5:
        return 1
    elif num_seg_raw < 10:
        return 2
    elif num_seg_raw < 40:
        return 3
    elif num_seg_raw < 70:
        return 4
    else:
        return 5


cpdef assemble_cluster(Cluster[::1] clt_view, dict cluster_data_by_tid, tuple block, object cmd_args, int extra_thread):
    """
    Assembles clusters into sequences using multiple steps, including primary assembly, polishing, and recalibration.

    Parameters:
        clt_view (Cluster[::1]): Memory view of the cluster array.
        cluster_data_by_tid (dict): A dictionary where the key is the target ID (tid) and the value is a tuple containing:
                                    - clt_arr: Array of clusters.
                                    - seg_arr: Array of segments.
        block (tuple): A tuple (start_idx, end_idx) specifying the range of clusters to process.
        cmd_args (object): Command-line arguments or configuration object containing file paths and parameters.
        extra_thread (int): Number of additional threads to use.
    """

    cdef int num_thread = 1 + extra_thread
    cdef int node_len = cmd_args.node_len
    cdef int i, min_edge, rounds
    cdef str cmd, prefix
    
    start_idx, end_idx = block
    for i in range(start_idx, end_idx):
        if is_lowfreq_clt(&clt_view[i]):
            continue

        min_edge = get_min_edge(clt_view[i].numSegRaw)
        if cmd_args.min_edge > 0:
            min_edge = cmd_args.min_edge

        # Step 1: Primary assembling
        prefix = "tmp_assm/{}_{}".format(clt_view[i].tid, clt_view[i].idx)
        rounds = 0
        while rounds < 5:
            rounds += 1
            cmd = "wtdbg2 -p 5 -k 15 -l 256 -e {0} -S 1 -A --rescue-low-cov-edges --node-len {1} --ctg-min-length {1} " \
                "--ctg-min-nodes 1 -q -t {2} -i {3}.fa -fo {3}".format(min_edge, node_len, num_thread, prefix)
            subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
            if os.path.isfile(f"{prefix}.ctg.lay.gz") == False:
                continue
            
            cmd = "wtpoa-cns -q -c 1 -t {0} -i {1}.ctg.lay.gz -fo {1}_assm.fa".format(num_thread, prefix)
            subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
            if os.path.isfile(f"{prefix}_assm.fa") and (os.path.getsize(f"{prefix}_assm.fa") != 0):
                break
        
        if (not os.path.isfile(f"{prefix}_assm.fa")) or (os.path.getsize(f"{prefix}_assm.fa") == 0):
            if output_read_as_assmbly(clt_view, cluster_data_by_tid, cmd_args, i) == 0:
                continue

        # Step 2: First round polishing
        cmd = "minimap2 -aY {0}_assm.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
            "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # Step 3: Second round polishing
        cmd = "minimap2 -aY {0}_assembled.fa {0}.fa | samtools sort | samtools view -bhS -F 3332 -o {0}_RawToAssm.bam && " \
            "samtools consensus --ff 3332 -m simple -c 0 -d 1 -H 0.9 {0}_RawToAssm.bam -o {0}_assembled.fa".format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')

        # Step 4: Recalibration
        recalibration(prefix, cmd_args)

        # Step 5: Remove N-bases
        cmd = (
            "sed 's/N//g' {0}_assembled.fa | "
            "awk '{{if($1~/^>/){{ctg=$1; a[ctg]=\"\"}} else{{a[ctg]=a[ctg]\"\"$1}}}} "
            "END{{for(ctg in a){{if(a[ctg]!=\"\"){{print ctg; print a[ctg]}}}}}}' > {0}_tmp.fa && "
            "mv {0}_tmp.fa {0}_assembled.fa"
        ).format(prefix)
        subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')


#####################
### Recalibration ###
#####################
cdef recalibration(str prefix, object cmd_args):
    # Step 1: Map raw reads to polished seqs
    os.rename(f"{prefix}_assembled.fa", f"{prefix}_polished.fa")
    cmd = "minimap2 -aY {0}_polished.fa {0}.fa | " \
        "samtools sort | samtools view -bhS -F 3332 -o {0}_RawToPolish.bam && " \
        "samtools index {0}_RawToPolish.bam".format(prefix)
    subprocess.run(cmd, stderr=subprocess.DEVNULL, shell=True, executable='/bin/bash')
    if not os.path.exists(f"{prefix}_RawToPolish.bam.bai"):
        os.rename(f"{prefix}_polished.fa", f"{prefix}_assembled.fa")
        return
    
    # Step 2: Load files
    cdef bytes input_fn = f"{prefix}_polished.fa".encode("utf-8")
    cdef faidx_t *input_fa = fai_load(input_fn)
    cdef BamFile input_bam = BamFile(f"{prefix}_RawToPolish.bam", "rb")
    cdef object output_fa = open(f"{prefix}_assembled.fa", "w")

    # Step 3: Initilize
    cdef int tid, num_ref = faidx_nseq(input_fa)
    cdef int ref_len, max_len = 0
    for tid in range(num_ref):
        ref_len = faidx_seq_len(input_fa, faidx_iseq(input_fa, tid))
        if ref_len > max_len:
            max_len = ref_len

    # Step 4: Recalibrate
    cdef object query_arr = np.zeros(5 * max_len, dtype=np.int32)
    cdef object ref_arr = np.zeros(5 * max_len, dtype=np.int32)
    cdef int ref_start, ref_end, skip_next
    cdef const char *ref_name_char
    cdef char *ref_seq_char
    cdef str ref_seq, recalibrated_seq

    for tid in range(num_ref):
        # Step 5: Load sequence
        ref_name_char = faidx_iseq(input_fa, tid)
        ref_len = faidx_seq_len(input_fa, ref_name_char)
        ref_seq_char = faidx_fetch_seq(input_fa, ref_name_char, 0, ref_len, &ref_len)
        ref_seq = ref_seq_char.decode('utf-8')
        
        # Step 6: Define homopolymer regions
        polymer_regions = get_polymer_regions(ref_seq, ref_len)
        
        # Step 7: Collect query homo-polymers
        iterator = Iterator(input_bam, tid)
        query_polymers = get_query_polymers(iterator, query_arr, ref_arr, polymer_regions, ref_seq)
        
        # Step 8: Output recalibrated sequence
        output_fa.write(">" + ref_name_char.decode('utf-8') + "\n")
        skip_next = 0
        for ref_start in polymer_regions:
            ref_end = polymer_regions[ref_start]

            # Ignore long homopolymer, un-covered region
            if (ref_end - ref_start + 1 > 20) or (ref_start not in query_polymers):
                output_fa.write(ref_seq[ref_start:ref_end+1])
                continue
            
            # If previous base/homo-polymer has been recalibrated, skip this one
            if skip_next:
                output_fa.write(ref_seq[ref_start:ref_end+1])
                skip_next = 0
                continue
            
            recalibrated_seq, skip_next = get_recalibrated_seq(query_polymers, ref_seq, ref_start, ref_end)
            output_fa.write(recalibrated_seq)
        
        output_fa.write("\n")
        del iterator

    # 9. Close files
    fai_destroy(input_fa)
    output_fa.close()
    input_bam.close()


cdef object get_polymer_regions(str ref_seq, int ref_len):
    cdef int start = 0
    cdef int end = 1
    cdef object regions = OrderedDict()

    while end < ref_len:
        if ref_seq[end] == ref_seq[start]:
            end += 1
            continue
        
        regions[start] = end - 1
        start = end
        end += 1
    
    # Final homo-polymer
    regions[start] = end - 1
    return regions


cdef object get_query_polymers(Iterator iterator, int[::1] query_arr, int[::1] ref_arr, object polymer_regions, str ref_seq):
    cdef int ref_pos, ref_start, ref_end
    cdef int query_pos, query_start, query_end
    cdef int i, return_value, num_pairs
    cdef int ref_len = len(ref_seq)
    cdef object query_polymers = OrderedDict()
    cdef str read_seq, query_seq
    
    while True:
        # 1. Load read
        return_value = iterator.next_record_by_tid()
        if return_value < 0:
            return query_polymers

        # 2. Get aligned-pairs
        read_seq = get_read_seq(iterator)
        num_pairs = get_aligned_pairs(iterator.bam_record, &query_arr[0], &ref_arr[0])
        ref_end = -2

        # 3. Collect homo-polymer lengths
        for i in range(num_pairs):
            query_pos = query_arr[i]
            ref_pos = ref_arr[i]

            # 4. Read cover homo-polymer start
            if ref_pos in polymer_regions:
                query_start = query_pos
                ref_start = ref_pos
                ref_end = polymer_regions[ref_pos]

                # Single base
                if ref_start == ref_end:
                    if query_start < 0:
                        query_start = -1
                        query_end = -2
                        query_seq = ""
                    else:
                        # Find real start
                        query_pos = check_left_side(query_arr, ref_arr, i-1, 0)
                        if query_pos >= 0:
                            query_start = query_pos + 1 
                        # Find real end
                        query_pos = check_right_side(query_arr, ref_arr, num_pairs, i+1, ref_len)
                        if query_pos >= 0:
                            query_end = query_pos - 1

                        query_seq = read_seq[query_start:query_end+1]
                    
                    if query_start == query_end:
                        extra_len = get_extra_len(query_start, read_seq, ref_seq[ref_start])
                    else:
                        extra_len = 0
                    
                    if ref_start in query_polymers:
                        query_polymers[ref_start].append((query_seq, query_end-query_start+1, extra_len))
                    else:
                        query_polymers[ref_start] = [(query_seq, query_end-query_start+1, extra_len)]
                    continue
                
                # Homo-polymer region
                if query_start < 0:
                    # Correction for boundary deletion
                    query_pos = check_right_side(query_arr, ref_arr, num_pairs, i+1, ref_end)
                    if query_pos < 0:
                        ref_end = -2
                    else:
                        query_start = query_pos
                else:
                    # Correction for boundary insertion
                    query_pos = check_left_side(query_arr, ref_arr, i-1, 0)
                    if query_pos >= 0:
                        query_start = query_pos + 1            
            
            # 5. Read cover homo-polymer end
            elif ref_pos == ref_end:
                query_end = query_pos
                if query_end < 0:
                    # Correction for boundary deletion
                    query_pos = check_left_side(query_arr, ref_arr, i-1, ref_start)
                    if query_pos < 0:
                        ref_end = -2
                        continue
                    else:
                        query_end = query_pos
                else:
                    query_pos = check_right_side(query_arr, ref_arr, num_pairs, i+1, ref_len)
                    if query_pos >= 0:
                        query_end = query_pos - 1
                
                if query_end - query_start == ref_end - ref_start:
                    extra_len = get_extra_len(query_start, read_seq, ref_seq[ref_start])
                else:
                    extra_len = 0
                
                # 6. Collect query homo-polymers
                if ref_start in query_polymers:
                    query_polymers[ref_start].append((read_seq[query_start:query_end+1], query_end - query_start + 1, extra_len))
                else:
                    query_polymers[ref_start] = [(read_seq[query_start:query_end+1], query_end - query_start + 1, extra_len)]
                
                ref_end = -2


cdef str get_read_seq(Iterator iterator):
    cdef int i, read_len = iterator.bam_record.core.l_qseq
    cdef bytes bSeq = PyBytes_FromStringAndSize(NULL, read_len)
    cdef char *cSeq = <char*>bSeq

    for i in range(read_len):
        cSeq[i] = seq_nt16_str[bam_seqi(bam_get_seq(iterator.bam_record), i)]

    return bSeq.decode('utf-8')


cdef int check_left_side(int[::1] query_arr, int[::1] ref_arr,  int i, int leftMost):
    while i >= 0:
        query_pos = query_arr[i]
        ref_pos = ref_arr[i]
        if (query_pos < 0) or (ref_pos < 0):
            i -= 1
            continue
        if ref_pos < leftMost:
            return -1
        return query_pos

    return -1


cdef int check_right_side(int[::1] query_arr, int[::1] ref_arr, int num_pairs,  int i, int rightMost):
    while i < num_pairs:
        query_pos = query_arr[i]
        ref_pos = ref_arr[i]
        if (query_pos < 0) or (ref_pos < 0):
            i += 1
            continue
        if ref_pos > rightMost:
            return -1
        return query_pos

    return -1


cdef int get_extra_len(int query_start, str read_seq, str polymerBase):
    if query_start == 0:
        return 0    
    if read_seq[query_start-1] == polymerBase:
        return 1
    else:
        return 0


cdef int get_most_common_len(object query_polymers, int ref_start, int ref_end):
    cdef object counter = Counter([p[1] + p[2] for p in query_polymers[ref_start]])
    cdef int ref_len = ref_end - ref_start + 1
    cdef int ref_count = counter[ref_len]
    cdef int alt_len, alt_count

    for alt_len, alt_count in counter.most_common(2):
        if alt_len != ref_len:
            break
    
    if alt_count <= ref_count:
        return ref_len
    else:
        return alt_len


cdef tuple get_recalibrated_seq(object query_polymers, str ref_seq, int ref_start, int ref_end):
    # 1. Get most common length
    cdef int most_common_len = get_most_common_len(query_polymers, ref_start, ref_end)
    if most_common_len == 0:
        return '', 1

    # 2. Initialize counter for each position
    cdef list position_counters = [defaultdict(int) for _ in range(most_common_len)]
    cdef object polymer
    cdef int i, query_len
    cdef str query_seq, base
    
    # 3. Count for each position
    for polymer in query_polymers[ref_start]:
        query_len = polymer[1]
        if query_len != most_common_len:
            continue
        
        query_seq = polymer[0]
        for i in range(most_common_len):
            base = query_seq[i]
            position_counters[i][base] += 1
    
    # 4. Get most common base for each position
    cdef list result = []
    cdef object counter
    cdef str ref_base = ref_seq[ref_start]
    cdef int ref_count
    if len(position_counters[0]) == 0:
        return ref_seq[ref_start:ref_end+1], 0
    else:
        for counter in position_counters:
            base = max(counter, key=counter.get)
            ref_count = counter[ref_base] # 0 if ref_base not exists
            if ref_count >= counter[base]:
                result.append(ref_base)
            else:
                result.append(base)

        if (ref_end - ref_start + 1) == most_common_len:
            return ''.join(result), 0
        else:
            return ''.join(result), 1