# LOCATE

LOCATE (Long-read to Characterize All Transposable Elements) is a mapping-based method using long-read whole genome sequencing data (ONT / PacBio) to detect, assemble, and genotype transposon insertions.

## 1. Installation

### 1.1 Installation by conda / mamba

You can install LOCATE using conda or mamba. If there's no `mamba` in your environment, we recommend that you start with the [Miniforge distribution](https://github.com/conda-forge/miniforge).

```shell
mamba install huzr::locate
```

### 1.2 Installation by container

Or, you can install through Docker or Singularity.

```shell
# Docker pull from dockerhub
docker pull huzr14830/locate
docker run huzr14830/locate --help

# Singularity pull from dockerhub
singularity pull locate.sif docker://huzr14830/locate
singularity run --cleanenv locate.sif --help

# Build locally
git clone git@github.com:red-t/LOCATE.git
cd LOCATE
docker build -t locate .
docker run locate --help
```

### 1.3 Installation from source

Or, you can install from source code.

```shell
git clone https://github.com/red-t/LOCATE.git
cd LOCATE
mamba env create -f environment.yml
mamba activate locate
make && make clean
locate --help
```

Note: There's a small test data set inside `tests/data`.

## 2. Download annotations

LOCATE compatible `annotations` and `models` can be downloaded from [here](https://users.wenglab.org/boxu/LOCATE/data.html).

## 3. Quick Start

The most common way to call transposon insertions from long read alignments (PacBio / ONT):

```shell
# For example, using GRCh38 as reference
locate -b sorted.bam -r GRCh38.rmsk.bed -g GRCh38.gap.bed \
       -C GRCh38.transposon.class -T GRCh38.transposon.fa \
       -R GRCh38_no_alt.fa -H GRCh38_HighFreq -L GRCh38_LowFreq \
       -G bayesian -F both -o output_path
```

**Note:**
- Currently, LOCATE requires alignment mapped by `minimap2 -Y` option, which uses soft clipping for supplementary alignments:

```shell
minimap2 -aYx $PRESET $REF $QUERY | samtools view -bhS - | samtools sort -o sorted.bam -
samtools index sorted.bam
```

## 4. Command-line options

| Flag | Long | Required | Default | Description |
|------|------|----------|---------|-------------|
| `-b` | `--bam` | Yes | | Genomic alignment BAM (aligned with `minimap2 -Y`) |
| `-C` | `--class` | Yes | | TE class file, tab-delimited (order matches TE FASTA) |
| `-T` | `--te_fn` | Yes | | TE consensus sequences FASTA |
| `-R` | `--ref_fa` | Yes | | Reference genome FASTA |
| `-H` | `--high` | Yes | | AutoGluon model for high-frequency insertions |
| `-L` | `--low` | Yes | | AutoGluon model for low-frequency insertions |
| `-r` | `--repeat` | No | "" | Repeat annotation BED (RepeatMasker output) |
| `-g` | `--gap` | No | "" | Gap annotation BED |
| `-B` | `--blacklist` | No | "" | Blacklist BED |
| `-o` | `--outpath` | No | ./ | Output directory |
| `-t` | `--num_thread` | No | 1 | Max number of threads |
| `-G` | `--genotyper` | No | bayesian | Genotyping method: `bayesian` or `threshold` |
| `-F` | `--output-format` | No | both | Output format: `tsv`, `vcf`, or `both` |
| `-e` | `--min_edge` | No | 0 | Min read depth for wtdbg2 edges (auto if 0) |
| `-n` | `--node_len` | No | 256 | wtdbg2 node length (multiple of 256) |
| `-l` | `--min_seg_len` | No | 100 | Min segment length to consider |
| `-d` | `--max_dist` | No | 50 | Max distance to merge breakpoints into a cluster |
| `-O` | `--overhang` | No | 200 | Min overhang length |
| `-v` | `--verbose` | No | INFO | Enable debug logging |

## 5. Pipeline

LOCATE processes the input BAM in 8 stages:

1. **Background information extraction** — computes average divergence, depth, and median read length from the BAM
2. **LTR size definition** — splits LTR TE consensus in half for paired mapping
3. **TE reference building** — builds a temporary indexed TE FASTA reference
4. **Cluster building** — parses CIGAR strings, extracts segments, merges nearby breakpoints, computes features, and filters clusters using AutoGluon ML models
5. **Local assembly** — assembles insertion sequences with wtdbg2, polishes with minimap2 + samtools consensus, and recalibrates homopolymer regions
6. **Sequence output** — outputs assembled sequences for high-frequency and pseudo-assemblies for low-frequency clusters
7. **Annotation** — maps assembled sequences to TE consensus library, annotates TE fragments, polyA/T tails, TSDs, and computes insertion frequency
8. **Output merging & genotyping** — merges all data and determines genotypes using a Beta-Binomial Bayesian model (or threshold-based method)

## 6. Output

By default, LOCATE produces both a tab-delimited file (`result.tsv`) and a VCF 4.3 file (`result.vcf`) in the output directory. Use the `-F` flag to select a single format.

### result.tsv

| Column | Value | Description |
|--------|-------|-------------|
| 1 | chrom | Chromosome |
| 2 | start | Insertion start site on reference sequence (0-based, included) |
| 3 | end | Insertion end site on reference sequence (0-based, not-included) |
| 4 | family | Transposon family of the insertion, separated by "," |
| 5 | frequency | Insertion allele frequency |
| 6 | strand | Orientation of the inserted transposon fragment |
| 7 | genotype | Genotype determined by the frequency (0/0, 0/1, 1/1) |
| 8 | genotype_quality | genotype quality |
| 9 | passed | Whether this insertion passes the post-filtering (True/False) |
| 10 | query_region | Annotated regions on the insertion sequence, format: "{+/-}:{start}-{end}" |
| 11 | target_region | Target regions of each query region, format: "{source}:{start}-{end}" |
| 12 | total_support | Total support reads of the insertion |
| 13 | tsd_seq | Annotated TSD sequence. "." if no annotated TSD |
| 14 | insertion_seq | Annotated insertion sequence, from the assembled sequence |
| 15 | upstream_seq | Upstream sequence of the insertion (same orientation as reference) |
| 16 | downstream_seq | Downstream sequence of the insertion (same orientation as reference) |
| 17 | extra_info | Extra information |

### result.vcf

Each insertion is represented as a VCF 4.3 record with `SVTYPE=INS`. The `INFO` field includes:

| INFO tag | Type | Description |
|----------|------|-------------|
| END | Integer | End position of the variant |
| SVTYPE | String | Structural variant type (INS) |
| SVLEN | Integer | Insertion length |
| FAMILY | String | TE family name(s) |
| AF | Float | Allele frequency |
| STRAND | String | Insertion strand orientation (+/-) |
| TSD | String | Target site duplication sequence |
| TE_CLASS | String | TE class (DNA/LTR/LINE/SINE/Retroposon/unknown) |
| TRUNCATION | String | Truncation status |
| RECONSTRUCTED_ENDS | String | Reconstructed ends status |
| HAS_POLYA | Integer | Has polyA tail (0/1) |
| HAS_TSD | Integer | Has target site duplication (0/1) |
| ASSEMBLED | Integer | Insertion was assembled (0/1) |
| SINGLETON | Integer | Singleton insertion (0/1) |
| SELF2SELF | Integer | Self-to-self insertion (0/1) |
| SOLO_LTR | Integer | Solo LTR (0/1) |
| SUPPORT | Integer | Total supporting reads |
| LEFT_CLIP | Integer | Left-clipped reads |
| SPANNING | Integer | Spanning/mid-insert reads |
| RIGHT_CLIP | Integer | Right-clipped reads |
| QV | Float | ML model probability |

The `FORMAT` field contains `GT` (genotype) and `GQ` (genotype quality). The `FILTER` field is `PASS` or `FAIL` based on post-filtering.
