#ifndef ANNO_UTILS_H
#define ANNO_UTILS_H

#include "htslib/faidx.h"
#include "io_utils.h"

/******************
 *** Structures ***
 ******************/

/// @brief Data container for annotation record
typedef struct Annotation
{
    int         idx;
    int         cltTid;
    int         cltIdx;
    int         query_start;
    int         query_end;
    uint8_t     strand;
    int         tid;
    int         ref_start;
    int         ref_end;
    uint32_t    flag;
    int         extra;
} __attribute__((packed)) Annotation;

/// @brief Data container for polyA
typedef struct PolyA
{
    int idx;
    int cltTid;
    int cltIdx;
    int leftAnnoStart;
    int leftIdx;
    int rightAnnoEnd;
    int rightIdx;
    int isA;
    int seqLen;
    int existPolyT;
} PolyA;


/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/
#define isRevAnno(anno) ((anno).strand == 1)

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fill_anno_arr(Cluster *clt, Annotation *anno_arr, uint32_t *class_arr, int idx);

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *clt, Annotation *anno_arr, int num_anno, PolyA *polyA);

/// @brief Find and record single polyA/polyT
int setPolyA(char *flankSeq, Annotation *anno_arr, Cluster *clt, int num_anno, PolyA *polyA);

/// @brief Add polyA candidate to anno_arr
void addPolyA(Annotation *anno_arr, int num_anno, Cluster *clt, PolyA *polyA, int position, int length);

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *clt, PolyA *polyA, Annotation *anno_arr, int num_anno);

/// @brief Adjust annotation position
void adjustAnno(Annotation *anno_arr, int num_anno, int leftDelta);

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annotate_tsd(Cluster *clt, Annotation *anno_arr, int num_anno);

/// @brief Find TSD and refine breakpoint
int setTsd(Cluster *clt, int localStart, int leftEnd, int rightStart);

/// @brief Set ins-seq structure based on annotations
void set_ins_structure(Cluster *clt, Annotation *anno_arr, int num_anno, uint32_t *class_arr, int *size_arr, int *ltr_arr);

/// @brief Compare function for sorting annotations
int compare(const void *a, const void *b);

/// @brief Check whether the ins-seq contains large gap
void checkGap(Annotation *anno_arr, int num_anno, Cluster *clt, int *size_arr);

/// @brief Get index of the left-most TE annotation
int getLeftIdx(Annotation *anno_arr, int num_anno);

/// @brief Get index of the right-most TE annotation
int getRightIdx(Annotation *anno_arr, int num_anno);

/// @brief Check whether the ins-seq contains valid polyA
void checkPolyA(Annotation *anno_arr, int num_anno, Cluster *clt);

/// @brief Check whether the ins-seq contains complete ends
void checkEnd(Annotation *anno_arr, int num_anno, Cluster *clt);

/// @brief Check which TE class the insertion belongs to
void checkTEClass(Annotation *anno_arr, int num_anno, Cluster *clt, uint32_t *class_arr);

/// @brief Check whether left-/right- assm-flank-seq contains valid polyT/A
void checkFlankPolyA(Annotation *anno_arr, int num_anno, Cluster *clt);

/// @brief Search polyT/polyA in left-/right- assm-flank-seq sequence
int searchFlankPolyA(char *flankSeq, int isA, int seqLen);

/// @brief Check whether the insertion is SOLO LTR
void checkSoloLtr(Annotation *anno_arr, int num_anno, Cluster *clt, int *size_arr, int *ltr_arr);


/**********************
 *** Annotation I/O ***
 **********************/

/// @brief Output formated annotation records
void output_annotations(Annotation *anno_arr, int num_anno, int start_idx, const char *te_fn);

/// @brief Change single annotation record into specified format
void formatSingleAnno(Annotation anno, char *queryTmp, char *refTmp, faidx_t *teFa, int *teTable, int *strandFlag);

/// @brief Write annotation for a single cluster
void writeSingleCltAnno(int strandFlag, int numTe, int *teTable, faidx_t *teFa, char *queryStr, char *refStr, FILE *fp, Annotation anno);

/// @brief Get cluster strand
char getCltStrand(int strandFlag);

/// @brief Generate a string which represents cluster TE-class
char *getCltClass(int numTe, int *teTable, faidx_t *teFa);

#endif // ANNO_UTILS_H