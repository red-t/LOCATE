#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <string.h>
#include "htslib/faidx.h"
#include "cluster_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/
#define is_lowqual_clt(clt) (((clt)->isInBlacklist) || ((clt)->probability <= 0.5) || ((clt)->numSeg >= 20))
#define is_lowfreq_clt(clt) ((clt)->cltType != 0)

/// @brief select a segment to output from cluster
int get_ouput_segidx(Cluster *clt, Segment *seg_arr, Args args);

/// @brief get extended region of the segment
void setTrimRegion(Segment *segment, int *start, int *end, int flank_size);


/*************************
 *** Flank Sequence IO ***
 *************************/

/// @brief Data container for reference flank region
typedef struct FlankRegion
{
    int start1;
    int start2;
    int end1;
    int end2;
} FlankRegion;

/// @brief output flank-seq and local-seq for all clusters
void extract_ref_flankseq(char *ref_fn, Cluster *clt_arr, int start_idx, int end_idx);

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *clt, FlankRegion *region);

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *clt, faidx_t *refFa, FlankRegion region);

/// @brief output +-500bp local-seq around cluster position for tsd annotation
void outputLocal(Cluster *clt, faidx_t *refFa, FlankRegion region);


/*****************************
 *** Insertion Sequence IO ***
 *****************************/
#define bamIsSup(bam) (((bam)->core.flag & BAM_FSUPPLEMENTARY) != 0)
#define isLeftFlank(bam) (strcmp(bam_get_qname((bam)), "0") == 0)
#define isClipInFlank(cigar, length) (bam_cigar_op((cigar)) == BAM_CSOFT_CLIP && bam_cigar_oplen((cigar)) > (length))

/// @brief Data container for insertion region on assembled contig
typedef struct InsRegion
{
    int rightMost;
    int leftMost;
    int tid1;
    int tid2;
    int len1;
    uint32_t cigar1;
    uint32_t cigar2;
    uint32_t flag;
} InsRegion;

/// @brief output insertion-seq and tsd-containing-seq from contig
void extract_ins_seq(Cluster *clt);

/// @brief define insertion-seq region by Flank-To-Assm alignments
void setInsRegion(Cluster *clt, InsRegion *region);

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region);

/// @brief output insertion-seq in FASTA format
void outputInsSeq(faidx_t *assmFa, Cluster *clt);

/// @brief get insertion-seq
char *getInsSeq(faidx_t *assmFa, Cluster *clt);

/// @brief output flank-seqs of insertion-seq from assembled-contig for re-deining insertion region
void outputAssmFlank(faidx_t *assmFa, Cluster *clt);

/// @brief refine and ouput insertion-seq for single-flank-mapped cases
void re_extract_ins_seq(Cluster *clt);

/// @brief re-set insertion-seq region
void reSetInsRegion(Cluster *clt, faidx_t *assmFa);

/*******************
 *** Cluster I/O ***
 *******************/

/// @brief Output formated cluster records
void output_clusters(Cluster *clt_arr, int start_idx, int end_idx, const char *ref_fn, const char *te_fn);

/// @brief Fetch TSD sequence from reference genome
char *fetchTsdSeq(faidx_t *refFa, Cluster *clt);

/// @brief Fetch flank sequence from temporary file
void fetchSeqs(Cluster *clt, char **insSeq, char **leftSeq, char **rightSeq);

/// @brief get adjusted insertion-seq
char *getAdjustSeq(faidx_t *assmFa, Cluster *clt);

#endif // IO_UTILS_H