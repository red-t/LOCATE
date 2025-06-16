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

## 3. Download annotations
LOCATE compatible `annotations/models` can be downloaded from [here](https://users.wenglab.org/boxu/LOCATE/data.html).

## 4. Quick Start

### 4.1 Most common way

The most common way to call transposon insertions from long read alignments (PacBio / ONT), you can use:

```shell
# For example, using GRCh38 as reference
locate -b sorted.bam -r GRCh38.rmsk.bed -g GRCh38.gap.bed -C GRCh38.transposon.class -T GRCh38.transposon.fa -R GRCh38_no_alt.fa -H GRCh38_HighFreq -L GRCh38_LowFreq -o output_path
```

**Note:**
- Currently, LOCATE requires alignment mapped by `minimap2 -Y` option, which use soft clipping for supplementary alignments, for example:

```shell
minimap2 -aYx $PRESET $REF $QUERY | samtools view -bhS - | samtools sort -o sorted.bam -
samtools index sorted.bam
```

### 4.2 No pretrained model

Currently, LOCATE provides pretrained models for GRCh38 and Dm6.
If no models are available for the genome assembly or species you are working with, one alternative is to use the existing models.

For example, the GRCh38 model can be used if your sequencing data was generated from an individual library, while the Dm6 model may be appropriate for data generated from a pooled library.

## 5. Output

The tab-delimited file `output_path/result.tsv` stores the result of LOCATE.

```shell
Column  Value               Description

1       chrom               chromosome
2       start               insertion start site on reference sequence (0-based, included)
3       end                 insertion end site on reference sequence (0-based, not-included)
4       family              transposon family of the insertion, separated by ","
5       frequency           insertion frequency.
6       strand              orientation of the inserted transposon fragment
7       genotype            genotype determined by the frequency (0/0, 0/1, 1/1)
8       passed              whether this insertion pass the post-filtering (True/False)
9       query_region        annotated regions on the insertion sequence, follow the pattern: "{+/-}:{start}-{end}"
10      target_region       target regions of each query region, follow the pattern: "{source}:{start}-{end}"
11      total_support       total support reads of the insertion
12      tsd_seq             annotated TSD sequence (corresponding to "chrom:start-end". "." if no annotated tsd)
16      insertion_seq       annotated insertion sequence, from the assembled sequence
17      upstream_seq        upstream sequence of the insertion sequence, from the assembled sequence (has the same orientation as the reference)
18      downstream_seq      downstream sequence of the insertion sequence, from the assembled sequence (has the same orientation as the reference)
19      extra_info          extra information
```
