import os
import logging
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from .FileIO import output_lowfreq_clusters_seq, output_reference_flank, merge_output
from .Cluster import build_cluster
from .Assemble import assemble_cluster
from .Annotate import annotate_cluster

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

###########################
### Helper Functions ###
###########################
cdef list divide_tasks(int total_tasks, int max_workers):
    """
    Divide tasks into blocks for parallel processing.

    Parameters:
        total_tasks (int): Total number of tasks.
        max_workers (int): Maximum number of workers.

    Returns:
        list[tuple[int, int]]: A list of task blocks represented as (start, end) tuples.
    """
    cdef int start = 0
    cdef int end = 0
    cdef list task_blocks = []
    cdef int base_tasks = total_tasks // max_workers
    cdef int remaining_tasks = total_tasks % max_workers

    while end < total_tasks:
        # Assign one extra task to workers while there are remaining tasks
        num_tasks = base_tasks + 1 if remaining_tasks > 0 else base_tasks
        if remaining_tasks > 0:
            remaining_tasks -= 1

        end = start + num_tasks
        task_blocks.append((start, end))
        start = end

    return task_blocks


cdef list allocate_threads(int total_threads, int num_tasks):
    """
    Allocate threads for tasks based on the total number of threads and tasks.

    Parameters:
        total_threads (int): Total number of threads available.
        num_tasks (int): Number of tasks to distribute threads across.

    Returns:
        list[int]: A list of thread counts for each task.
    """
    base_thread = total_threads // num_tasks
    extra_threads = total_threads % num_tasks
    return [base_thread + 1 if i < extra_threads else base_thread for i in range(num_tasks)]


###########################
### Main Pipeline Logic ###
###########################
cpdef object run_in_parallel(object cmd_args):
    """
    Execute the main pipeline in parallel.

    Parameters:
        cmd_args (object): Command-line arguments object.

    Returns:
        dict: All cluster data.
    """
    cdef dict cluster_data_by_tid = {}, background_info

    # 1. Get Background Info
    background_info = get_background_info(cmd_args.genome_bam_fn, cmd_args.num_thread)
    logger.info("Background divergence: %f", background_info["average_divergence"])
    logger.info("Background depth: %f", background_info["average_depth"])
    logger.info("Background read length: %f", background_info["median_read_length"])

    # 2. Define LTR size
    define_ltr_size(cmd_args)

    # 3. Build TE reference
    build_te_reference(cmd_args)

    # 4. Determine max_workers
    cdef int num_chromosomes = background_info["num_chromosomes"]
    cdef int max_workers = min(cmd_args.num_thread, num_chromosomes)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 5. Build Clusters
        cluster_data_by_tid = build_clusters(executor, cmd_args, background_info, max_workers)

        # 6. Local Assembly
        high_quality_clusters = get_high_quality_clusters(cluster_data_by_tid)
        assemble_clusters(executor, cmd_args, high_quality_clusters, cluster_data_by_tid, max_workers)

        # 7. Output Sequences
        output_sequences(executor, cmd_args, cluster_data_by_tid, high_quality_clusters, num_chromosomes, max_workers)

        # 8. Annotate Clusters
        annotate_clusters(executor, cmd_args, high_quality_clusters, max_workers)

    # 9. Merge Output
    merge_output()
    return cluster_data_by_tid


cdef dict get_background_info(str genome_bam_fn, int num_thread):
    """
    Extract background information from the genome BAM file.

    Parameters:
        genome_bam_fn (str): Path to the genome BAM file.
        num_thread (int): Number of threads to use.

    Returns:
        dict: A dictionary containing background information such as
              number of chromosomes, divergence, depth, and read length.
    """
    cdef BamFile genome_bam = BamFile(genome_bam_fn, "rb", num_thread)
    cdef int i, tid = 0, max_chromosome_length = 0

    # Find the chromosome with the maximum length
    for i in range(genome_bam.header.n_targets):
        chrom_length = sam_hdr_tid2len(genome_bam.header, i)
        if max_chromosome_length < chrom_length:
            max_chromosome_length = chrom_length
            tid = i

    cdef Iterator iterator = Iterator(genome_bam, tid)
    cdef int size = 0, capacity = 200000, threshold = int(0.9 * capacity)
    cdef object read_length_array = np.zeros(capacity, dtype=np.int32)
    cdef object divergence_array = np.zeros(capacity, dtype=np.float32)
    cdef float[::1] divergence_view = divergence_array
    cdef int[::1] read_length_view = read_length_array
    cdef int return_value, alignment_length
    cdef int64_t total_alignment_length = 0
    cdef float divergence
    cdef dict background_info = {}

    while True:
        return_value = iterator.next_record_by_tid()
        if return_value < 0:
            # Compute background information
            if size == 0:
                raise ValueError("No valid alignments found in the BAM file.")
            background_info["num_chromosomes"] = genome_bam.header.n_targets
            background_info["average_divergence"] = np.mean(divergence_array[:size])
            background_info["average_depth"] = float(total_alignment_length) / max_chromosome_length
            background_info["median_read_length"] = np.median(read_length_array[:size])
            genome_bam.close()
            del iterator
            del genome_bam
            return background_info

        if bam_is_invalid(iterator.bam_record):
            continue

        if size >= threshold:
            # Dynamically resize arrays
            capacity = int(capacity * 1.5)
            threshold = int(0.9 * capacity)
            divergence_array.resize((capacity,), refcheck=False)
            read_length_array.resize((capacity,), refcheck=False)
            divergence_view = divergence_array
            read_length_view = read_length_array

        # Extract alignment length and divergence
        get_mapping_length_and_divergence(&alignment_length, &divergence, iterator.bam_record)
        total_alignment_length += alignment_length
        divergence_view[size] = divergence
        read_length_view[size] = iterator.bam_record.core.l_qseq
        size += 1


cdef void define_ltr_size(cmd_args):
    """
    Define LTR size based on input arguments.
    """
    cdef bytes te_fn = cmd_args.te_fn.encode()
    cdef bytes class_fn = cmd_args.class_fn.encode()
    define_ltr(te_fn, class_fn)


cdef void build_te_reference(object cmd_args):
    """
    Build a temporary TE reference file.

    Parameters:
        cmd_args (object): Command-line arguments object containing paths and parameters.
    """
    try:
        # Create a temporary FASTA file
        with open("tmp_build/tmp.fa", "w") as output_fasta:
            output_fasta.write(">0\nAAAAAAAAAAAAAA\n")

        # Define input and output file paths
        query_filename = "tmp_build/tmp.fa"
        output_filename = "tmp_build/tmp.bam"

        # Construct the command for minimap2 and samtools
        command = (
            "minimap2 -t 1 -aY {} {} | samtools view -bhS -o {} -".format(
                cmd_args.te_fn, query_filename, output_filename
            )
        )

        # Execute the command
        subprocess.run(command, stderr=subprocess.DEVNULL, shell=True, executable="/bin/bash", check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error occurred while building TE reference: {e}")


def build_clusters(executor, cmd_args, background_info, max_workers):
    """
    Build clusters in parallel.

    Parameters:
        executor (ProcessPoolExecutor): Executor for parallel processing.
        cmd_args (object): Command-line arguments object.
        background_info (dict): Background information.
        max_workers (int): Maximum number of workers.

    Returns:
        dict: Cluster data.
    """
    cdef dict cluster_data_by_tid = {}
    thread_allocation = allocate_threads(cmd_args.num_thread - max_workers, background_info["num_chromosomes"])

    subprocess_tuple = set()
    for tid, extra_thread in enumerate(thread_allocation):
        subprocess_tuple.add(executor.submit(
            build_cluster, background_info["average_divergence"], background_info["average_depth"],
            background_info["median_read_length"], cmd_args, tid, extra_thread
        ))

    for subprocess in as_completed(subprocess_tuple):
        cluster_data = subprocess.result()
        cluster_data_by_tid.update(cluster_data)

    return cluster_data_by_tid


def assemble_clusters(executor, cmd_args, high_quality_clusters, cluster_data_by_tid, max_workers):
    """
    Perform local assembly for high-quality clusters.

    Parameters:
        executor (ProcessPoolExecutor): Executor for parallel processing.
        cmd_args (object): Command-line arguments object.
        high_quality_clusters (object): High-quality clusters.
        cluster_data_by_tid (dict): All cluster data.
        max_workers (int): Maximum number of workers.
    """
    task_blocks = divide_tasks(high_quality_clusters.shape[0], max_workers)
    thread_allocation = allocate_threads(cmd_args.num_thread - max_workers, len(task_blocks))

    subprocess_tuple = set()
    for block, extra_thread in zip(task_blocks, thread_allocation):
        subprocess_tuple.add(executor.submit(
            assemble_cluster, high_quality_clusters, cluster_data_by_tid, block, cmd_args, extra_thread
        ))

    for subprocess in as_completed(subprocess_tuple):
        subprocess.result()


def output_sequences(executor, cmd_args, cluster_data_by_tid, high_quality_clusters, num_chromosomes, max_workers):
    """
    Output sequences for clusters.

    Parameters:
        executor (ProcessPoolExecutor): Executor for parallel processing.
        cmd_args (object): Command-line arguments object.
        cluster_data_by_tid (dict): All cluster data.
        high_quality_clusters (object): High-quality clusters.
        num_chromosomes (int): Number of chromosomes.
        max_workers (int): Maximum number of workers.
    """
    # Output sequences for clusters without assembly
    thread_allocation = allocate_threads(cmd_args.num_thread - max_workers, num_chromosomes)

    subprocess_tuple = set()
    for tid, extra_thread in enumerate(thread_allocation):
        subprocess_tuple.add(executor.submit(
            output_lowfreq_clusters_seq, cluster_data_by_tid[tid][0], cluster_data_by_tid[tid][1], cmd_args, tid, extra_thread
        ))

    for subprocess in as_completed(subprocess_tuple):
        subprocess.result()

    # Output reference flank sequences for high-quality clusters
    task_blocks = divide_tasks(high_quality_clusters.shape[0], max_workers)
    subprocess_tuple = set([executor.submit(
        output_reference_flank, high_quality_clusters, cluster_data_by_tid, block, cmd_args
    ) for block in task_blocks])

    for subprocess in as_completed(subprocess_tuple):
        subprocess.result()


def annotate_clusters(executor, cmd_args, high_quality_clusters, max_workers):
    """
    Output results for clusters and high-quality clusters.

    Parameters:
        executor (ProcessPoolExecutor): Executor for parallel processing.
        cmd_args (object): Command-line arguments object.
        high_quality_clusters (object): High-quality clusters.
        max_workers (int): Maximum number of workers.
    """
    task_blocks = divide_tasks(high_quality_clusters.shape[0], max_workers)
    subprocess_tuple = set([executor.submit(
        annotate_cluster, high_quality_clusters, block, cmd_args
    ) for block in task_blocks])

    for subprocess in as_completed(subprocess_tuple):
        subprocess.result()