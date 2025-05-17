#!/usr/bin/python
import os
import argparse
import shutil
import logging
from .ParallelTaskExecutor import run_in_parallel

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

###########################
### Helper Functions ###
###########################
def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="LOCATE: A tool for transposon insertion detection.")
    parser.add_argument('-b', '--bam', dest='genome_bam_fn', type=str, required=True,
                        help='Genomic alignment in BAM format, aligned by minimap2 -Y')
    parser.add_argument('-r', '--repeat', dest='repeat_fn', type=str, default='',
                        help='Repeat annotation in BED format, annotated by RepeatMasker')
    parser.add_argument('-g', '--gap', dest='gap_fn', type=str, default='',
                        help='Gap annotation in BED format')
    parser.add_argument('-C', '--class', dest='class_fn', type=str, required=True,
                        help='Tab-delimited TE class file, the order should be consistent with the TE consensus fasta')
    parser.add_argument('-T', '--te_fn', dest='te_fn', type=str, required=True,
                        help='Transposon consensus sequences in FASTA format')
    parser.add_argument('-R', '--ref_fa', dest='ref_fn', type=str, required=True,
                        help='Reference genome in FASTA format')
    parser.add_argument('-H', '--high', dest='high_freq_model', type=str, required=True,
                        help='Path to autogluon model for high frequency insertions')
    parser.add_argument('-L', '--low', dest='low_freq_model', type=str, required=True,
                        help='Path to autogluon model for low frequency insertions')
    parser.add_argument('-B', '--blacklist', dest='blacklist_fn', type=str, default='',
                        help='Blacklist in BED format')
    parser.add_argument('-e', '--min_edge', dest='min_edge', type=int, default=0,
                        help='Min read depth of a valid edge for wtdbg2, automatically estimated if not provided')
    parser.add_argument('-n', '--node_len', dest='node_len', type=int, default=256,
                        help='Node length for wtdbg2, times of 256bp')
    parser.add_argument('-o', '--outpath', dest='out_path', type=str, default='./',
                        help='Output directory')
    parser.add_argument('-t', '--num_thread', dest='num_thread', type=int, default=1,
                        help='Max number of extra threads to use in each sub-process')
    parser.add_argument('-l', '--min_seg_len', dest='min_seg_len', type=int, default=100,
                        help='Min segment length, reads with clip-/insert-segment < min_seg_len will be ignored')
    parser.add_argument('-d', '--max_dist', dest='max_dist', type=int, default=50,
                        help='Reads (breakpoints) within maxDist will be merged as a cluster')
    parser.add_argument('-O', '--overhang', dest='overhang', type=int, default=200,
                        help='Min overhang length, reads with genomic-mapping-length < overhang will be ignored')
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate the parsed arguments.

    Parameters:
        args (argparse.Namespace): Parsed arguments.

    Raises:
        FileNotFoundError: If required files or directories are missing.
        ValueError: If argument values are invalid.
    """
    required_files = [
        args.genome_bam_fn, args.class_fn, args.te_fn, args.ref_fn,
        args.high_freq_model, args.low_freq_model
    ]
    optional_files = [args.repeat_fn, args.gap_fn, args.blacklist_fn]

    for file in required_files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"Required file not found: {file}")

    for file in optional_files:
        if file and not os.path.exists(file):
            raise FileNotFoundError(f"Optional file not found: {file}")

    if args.min_edge < 0:
        raise ValueError(f"min_edge must be >= 0: {args.min_edge}")
    if args.node_len <= 0 or args.node_len % 256 != 0:
        raise ValueError(f"node_len must be a positive multiple of 256: {args.node_len}")
    if args.min_seg_len <= 0:
        raise ValueError(f"min_seg_len must be > 0: {args.min_seg_len}")
    if args.max_dist <= 0:
        raise ValueError(f"max_dist must be > 0: {args.max_dist}")
    if args.overhang <= 0:
        raise ValueError(f"overhang must be > 0: {args.overhang}")


def prepare_directories(out_path: str) -> None:
    """
    Prepare output and temporary directories.

    Parameters:
        out_path (str): Output directory path.

    Raises:
        OSError: If directory creation fails.
    """
    try:
        os.makedirs(out_path, exist_ok=True)
        os.chdir(out_path)
        logger.info(f"Changed working directory to: {os.getcwd()}")

        for tmp_dir in ["tmp_build", "tmp_assm", "tmp_anno"]:
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)
            logger.info(f"Created temporary directory: {tmp_dir}")
    except Exception as e:
        raise OSError(f"Failed to prepare directories: {e}")


def summarize_results(tid_to_result: dict) -> None:
    """
    Summarize the results of the pipeline.

    Parameters:
        tid_to_result (dict): Results from the pipeline.
    """
    num_clusters = sum(result[0].shape[0] for result in tid_to_result.values())
    num_segments = sum(result[1].shape[0] for result in tid_to_result.values())

    logger.info(f"Total clusters constructed: {num_clusters}")
    logger.info(f"Total segments constructed: {num_segments}")


###########################
### Main Function ###
###########################
def main():
    """
    Main entry point for the script.
    """
    # 1. Parse command-line arguments
    args = parse_args()

    # 2. Validate arguments
    validate_args(args)

    # 3. Convert paths to absolute paths
    args.genome_bam_fn = os.path.abspath(args.genome_bam_fn)
    args.class_fn = os.path.abspath(args.class_fn)
    args.te_fn = os.path.abspath(args.te_fn)
    args.ref_fn = os.path.abspath(args.ref_fn)
    args.high_freq_model = os.path.abspath(args.high_freq_model) + "/"
    args.low_freq_model = os.path.abspath(args.low_freq_model) + "/"
    if args.repeat_fn:
        args.repeat_fn = os.path.abspath(args.repeat_fn)
    if args.gap_fn:
        args.gap_fn = os.path.abspath(args.gap_fn)
    if args.blacklist_fn:
        args.blacklist_fn = os.path.abspath(args.blacklist_fn)

    # 4. Prepare directories
    prepare_directories(args.out_path)

    # 5. Run the pipeline
    tid_to_result = run_in_parallel(args)

    # 6. Summarize results
    summarize_results(tid_to_result)


if __name__ == "__main__":
    main()