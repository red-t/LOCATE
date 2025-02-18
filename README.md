# LOCATE

**(TO DO)** A brief project description describing the purpose and function of the project.

## 1. Dependencies

**(TO DO)** Some dependencies update frequently, especially machine learning related packages, we need to find the stable version(s) of some "key pakages", or host a conda/docker container.

### 1.1 Miniforge

If there's no `mamba` in your environment, we recommend that you start with the [Miniforge distribution](https://github.com/conda-forge/miniforge).

Don't forget to add `bioconda` channel after installation:

```shell
mamba config --add channels bioconda
```

### 1.2 scikit-learn (1.3.2)

```shell
mamba create -n locate python=3.10.13
mamba activate locate

# Must be 1.3.2 (should be consist with the version used for model training)
mamba install -c conda-forge scikit-learn=1.3.2
```

### 1.3 [AutoGluon (1.0.0)](https://auto.gluon.ai/stable/install.html)

```shell
mamba install -c conda-forge autogluon=1.0.0
```

### 1.4 Samtools

```shell
# 1.17 is not available now
mamba install -c bioconda samtools=1.17

# 1.21 is available
mamba install -c bioconda samtools=1.21
mamba install -c bioconda samtools
```

### 1.5 Minimap2

```shell
# 2.26 is not available now
mamba install -c bioconda minimap2=2.26

# 2.28-r1209, 2.1.1 are available
mamba install -c bioconda minimap2=2.1.1
mamba install -c bioconda minimap2
```

### 1.6 Cython (3.0.6)

```shell
pip install cython==3.0.6
```

## 2. LOCATE Installation

**(TO DO)** Only installing from source code is available now, maybe we should host a conda/docker container, and find a way to share our models.

### 2.1 Clone LOCATE repository

```shell
git clone git@github.com:red-t/LOCATE.git
```

### 2.2 Compile

```shell
cd LOCATE
make
make clean
```

### 2.3 Download models

**(TO DO)** Pre-trained models for LOCATE can be downloaded from [link]().

Now all the models are placed in the following paths:

```shell
# High frequency model for GRCh38
/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/For_ML/germ_v4/Dedup/GRCh38_G_1V1
# Low frequency model for GRCh38
/data/tusers/zhongrenhu/for_SMS/dna/simulation/GRCh38/For_ML/soma_v4/Dedup/GRCh38_S_1V30

# High frequency model for dm6
/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm6/For_ML/germ_v4/Dedup/dm6_G_1V1
# Low frequency model for dm6
/data/tusers/zhongrenhu/for_SMS/dna/simulation/dm6/For_ML/soma_v4/Dedup/dm6_S_1V30
```

## 3. Usage

```shell
python LOCATE.py -b BAM -r REPEAT.bed -g GAP.bed -B BlackList.bed -T TE.fa --germ GERM_MODEL --soma SOMA_MODEL -p NUM_PROCESS -t NUM_THREAD

usage: LOCATE.py [-h] -b GENOMEBAMFN [-r REPEATFN] [-g GAPFN] -C CLASSFN -T TEFN -R REFFN -H HIGHFREQMODEL -L LOWFREQMODEL [-B BLACKLISTFN] [-e MINEDGE] [-n NODELEN] [-o OUTPATH] [-p NUMPROCESS] [-t NUMTHREAD] [-l MINSEGLEN] [-d MAXDISTANCE] [-O MINOVERHANG]

options:
  -h, --help            show this help message and exit
  -b GENOMEBAMFN, --bam GENOMEBAMFN
                        Genomic alignment in BAM format, mapped by: minimap2 -aYx map-hifi/map-pb/map-ont
  -r REPEATFN, --repeat REPEATFN
                        Repeat annotation in BED format, annotated by RepeatMasker
  -g GAPFN, --gap GAPFN
                        Gap annotation in BED format
  -C CLASSFN, --class CLASSFN
                        Tab-delimited TE class file, the order should be consistent with the TE consensus fasta
  -T TEFN, --teFn TEFN  Transposon consensus sequences in FASTA format
  -R REFFN, --refFa REFFN
                        Reference genome in FASTA format
  -H HIGHFREQMODEL, --high HIGHFREQMODEL
                        Path to the directory of high frequency model
  -L LOWFREQMODEL, --low LOWFREQMODEL
                        Path to the directory of low frequency model
  -B BLACKLISTFN, --blacklist BLACKLISTFN
                        Blacklist in BED format
  -e MINEDGE, --minEdge MINEDGE
                        Min read depth of a valid edge for wtdbg2, automatically estimated if not provided
  -n NODELEN, --nodeLen NODELEN
                        Node length for wtdbg2, times of 256bp
  -o OUTPATH, --outpath OUTPATH
                        Output directory
  -p NUMPROCESS, --numProcess NUMPROCESS
                        Max number of worker processes
  -t NUMTHREAD, --numThread NUMTHREAD
                        Max number of extra threads to use in each sub-process
  -l MINSEGLEN, --minSegLen MINSEGLEN
                        Min segment length, reads with clip-/insert-segment < minSegLen will be ignored
  -d MAXDISTANCE, --maxDist MAXDISTANCE
                        Reads (breakpoints) within maxDist will be merged as a cluster
  -O MINOVERHANG, --overhang MINOVERHANG
                        Min overhang length, reads with genomic-mapping-length < overhang will be ignored
```


**Test Case For TianXiong:**

```shell
# After installation, I successfully tested 2 cases with simulated data

# dm6
bash /home/zhongrenhu/test/run_dm6_Germ.sh
cd workdir_dm6 && bash /home/zhongrenhu/test/results/summary_2.sh tmp_anno/result.txt > summary.txt

# GRCh38
bash /home/zhongrenhu/test/run_GRCh38_Germ.sh
cd workdir_GRCh38 && bash /home/zhongrenhu/test/results/summary_2.sh tmp_anno/result.txt > summary.txt
```

## 4. Output

Tab-delimited file `{OUTPUT_DIR}/tmp_anno/result.txt` stores the result of LOCATE.

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
