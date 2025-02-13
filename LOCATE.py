#!/usr/bin/python
import os
import argparse
import shutil
from LOCATE.ParallelTaskExecutor import runInParallel

def parseArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', '--bam', dest='genomeBamFn', type=str, required=True,
                                    help='Genomic alignment in BAM format, aligned by minimap2 -Y')
    parser.add_argument('-r', '--repeat', dest='repeatFn', type=str, default='',
                                    help='Repeat annotation in BED format, annotated by RepeatMasker')
    parser.add_argument('-g', '--gap', dest='gapFn', type=str, default='',
                                    help='Gap annotation in BED format')
    parser.add_argument('-C', '--class', dest='classFn', type=str, required=True,
                                    help='Tab-delimited TE class file, the order should be consistent with the TE consensus fasta')
    parser.add_argument('-T', '--teFn', dest='teFn', type=str, required=True,
                                    help='Transposon consensus sequences in FASTA format')
    parser.add_argument('-R', '--refFa', dest='refFn', type=str, required=True,
                                    help='Reference genome in FASTA format')
    parser.add_argument('-H', '--high', dest='highFreqModel', type=str, required=True,
                                    help='Path to autogluon model for high frequency insertions')
    parser.add_argument('-L', '--low', dest='lowFreqModel', type=str, required=True,
                                    help='Path to autogluon model for low frequency insertions')
    parser.add_argument('-B', '--blacklist', dest='blacklistFn', type=str, default='',
                                    help='Blacklist in BED format')
    parser.add_argument('-e', '--minEdge', dest='minEdge', type=int, default=0,
                                    help='Min read depth of a valid edge for wtdbg2, automatically estimated if not provided')
    parser.add_argument('-n', '--nodeLen', dest='nodeLen', type=int, default=256,
                                    help='Node length for wtdbg2, times of 256bp')
    parser.add_argument('-o', '--outpath', dest='outPath', type=str, default='./',
                                    help='Output directory')
    parser.add_argument('-p', '--numProcess', dest='numProcess', type=int, default=1,
                                    help='Max number of worker processes')
    parser.add_argument('-t', '--numThread', dest='numThread', type=int, default=1,
                                    help='Max number of extra threads to use in each sub-process')
    parser.add_argument('-l', '--minSegLen', dest='minSegLen', type=int, default=100,
                                    help='Min segment length, reads with clip-/insert-segment < minSegLen will be ignored')
    parser.add_argument('-d', '--maxDist', dest='maxDistance', type=int, default=50,
                                    help='Reads (breakpoints) within maxDist will be merged as a cluster')
    parser.add_argument('-O', '--overhang', dest='minOverhang', type=int, default=200,
                                    help='Min overhang length, reads with genomic-mapping-length < overhang will be ignored')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # 1. Parse command args
    cmdArgs = parseArgs()
    
    # 2. Check args
    if not os.path.exists(cmdArgs.genomeBamFn):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.genomeBamFn}")
    if not os.path.exists(cmdArgs.repeatFn) and cmdArgs.repeatFn:
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.repeatFn}")
    if not os.path.exists(cmdArgs.gapFn) and cmdArgs.gapFn:
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.gapFn}")
    if not os.path.exists(cmdArgs.classFn):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.classFn}")
    if not os.path.exists(cmdArgs.teFn):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.teFn}")
    if not os.path.exists(cmdArgs.refFn):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.refFn}")
    if not os.path.exists(cmdArgs.highFreqModel):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.highFreqModel}")
    if not os.path.exists(cmdArgs.lowFreqModel):
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.lowFreqModel}")
    if not os.path.exists(cmdArgs.blacklistFn) and cmdArgs.blacklistFn:
        raise FileNotFoundError(f"No such file or directory: {cmdArgs.blacklistFn}")
    if cmdArgs.minEdge < 0:
        raise ValueError(f"minEdge less than 0: {cmdArgs.minEdge}")
    if cmdArgs.nodeLen <= 0 or (cmdArgs.nodeLen % 256) != 0:
        raise ValueError(f"nodeLen should be times of 256: {cmdArgs.nodeLen}")
    if cmdArgs.minSegLen <= 0:
        raise ValueError(f"minSegLen less or equal to 0: {cmdArgs.minSegLen}")
    if cmdArgs.maxDistance <= 0:
        raise ValueError(f"maxDistance less or equal to 0: {cmdArgs.maxDistance}")
    if cmdArgs.minOverhang <= 0:
        raise ValueError(f"minOverhang less or equal to 0: {cmdArgs.minOverhang}")
    
    # 3. Convert to absolute path
    cmdArgs.genomeBamFn = os.path.abspath(cmdArgs.genomeBamFn)
    if cmdArgs.repeatFn:
        cmdArgs.repeatFn = os.path.abspath(cmdArgs.repeatFn)
    if cmdArgs.gapFn:
        cmdArgs.gapFn = os.path.abspath(cmdArgs.gapFn)
    cmdArgs.classFn = os.path.abspath(cmdArgs.classFn)
    cmdArgs.refFn = os.path.abspath(cmdArgs.refFn)
    cmdArgs.teFn = os.path.abspath(cmdArgs.teFn)
    cmdArgs.refFn = os.path.abspath(cmdArgs.refFn)
    cmdArgs.highFreqModel = os.path.abspath(cmdArgs.highFreqModel) + "/"
    cmdArgs.lowFreqModel = os.path.abspath(cmdArgs.lowFreqModel) + "/"
    if cmdArgs.blacklistFn:
        cmdArgs.blacklistFn = os.path.abspath(cmdArgs.blacklistFn)
    
    # 4. Change working directory
    try:
        if os.path.exists(cmdArgs.outPath) == False:
            os.makedirs(cmdArgs.outPath)
        os.chdir(cmdArgs.outPath)
        print("Current working directory:", os.getcwd())
    except:
        raise OSError(f"Failed to change working dierectory: {cmdArgs.outPath}")
    
    # 5. Make temporary directories
    try:
        if os.path.exists("tmp_build"):
            shutil.rmtree("tmp_build")
            os.makedirs("tmp_build")
        else:
            os.makedirs("tmp_build")

        if os.path.exists("tmp_assm"):
            shutil.rmtree("tmp_assm")
            os.makedirs("tmp_assm")
        else:
            os.makedirs("tmp_assm")

        if os.path.exists("tmp_anno"):
            shutil.rmtree("tmp_anno")
            os.makedirs("tmp_anno")
        else:
            os.makedirs("tmp_anno")
    except:
        raise OSError("Failed to temporary directories")
    
    # 6. Run all jobs
    tidToResult = runInParallel(cmdArgs)
    
    # 7. Simple overview
    print("当前版本构建的 clusters 总数:")
    numCluster = 0
    for i in range(len(tidToResult)):
        numCluster += tidToResult[i][0].shape[0]

    print(numCluster, "\n")

    print("当前版本构建的 segments 总数:")
    numSeg = 0
    for i in range(len(tidToResult)):
        numSeg += tidToResult[i][1].shape[0]
    
    print(numSeg)