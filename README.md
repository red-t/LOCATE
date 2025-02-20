# LOCATE

LOCATE (Long-read to Characterize All Transposable Elements) is a mapping-based method using long-read whole genome sequencing data (ONT / PacBio) to detect and assemble transposon insertions.

## 1. Installation by conda / mamba

You can install LOCATE using conda or mamba. If there's no `mamba` in your environment, we recommend that you start with the [Miniforge distribution](https://github.com/conda-forge/miniforge).

```shell
mamba create -n locate
mamba activate locate
mamba install locate
```

## 2. Installation from source

Or, you can install from source code.

### 2.1 Preparing dependencies

```shell
mamba create -n locate python=3.10.13 cython=3.0.6 scikit-learn=1.3.2 autogluon=1.0.0 samtools=1.21 minimap2=2.28 wtdbg=2.5
mamba activate locate
```
### 2.2 Clone the repository

```shell
git@github.com:red-t/LOCATE.git
cd LOCATE
```

### 2.3 Build LOCATE

```shell
make
make clean
```

## 3. Quick Start
To call transposon insertions from long read alignments (PacBio / ONT), you can use:

```shell
locate -b sorted.bam -r sorted_rmsk.bed -r sorted_gap.bed -C transposon.class -T transposon.fa -R reference.fa -H highfreq_model_dir -L lowfreq_model_dir -o output_dir
```

**Note:**
- Currently, LOCATE compatible `annotations/models` for `GRCh38/Dm6` can be found [here]().
- Currently, LOCATE requires alignment mapped by `minimap2 -Y` option, which use soft clipping for supplementary alignments, for example:

```shell
minimap2 -aYx $PRESET $REF $QUERY | samtools view -bhS - | samtools sort -o sorted.bam -
samtools index sorted.bam
```

## 4. Output

The tab-delimited file `${output_dir}/result.txt` stores the result of LOCATE.

```shell
Column  Value               Description

1       chrom               chromosome
2       refStart            insertion start on reference sequence (0-based, included)
3       refEnd              insertion end on reference sequence (0-based, not-included)
4       class               transposon class of the insertion, separated by ","
5       predProb            predicted probability to be TP by ML models, useless
6       orientation         orientation of the inserted transposon fragment
7       id                  insertion id
8       annoReg             annotated regions on the insertion sequence, follow the pattern: "{orientation}:{start}-{end}"
9       annoInfo            annotation of each region, follow the pattern: "{Source}:{start}-{end}"
10      numSupport          number of support reads
11      numLeftClip         number of left-clipping support reads
12      numSpan             number of spanning support reads
13      numRightClip        number of right-clipping support reads
14      isAssembled         a flag indicating whether the sequence is "assembled (1)" or "one of the support reads (0)"
15      tsdSeq              annotated TSD sequence, from reference genome ("ref:refStart-refEnd". "." if no annotated tsd)
16      insSeq              annotated insertion sequence, from the assembled result
17      upSeq               upstream sequence on the left of insSeq, from the assembled result (has the same orientation as reference)
18      downSeq             downstream sequence on the right of insSeq, from the assembled result (has the same orientation as reference)
19      flag                bitwise flag (refer to "LOCATE/src/cluster_utils.h#L16-L41")
```
