#ifndef SEG_UTILS_H
#define SEG_UTILS_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "htslib/sam.h"
#include "AIList.h"

/******************
 *** Structures ***
 ******************/

/*! @brief Data container for segment extracted from CIGAR.
 @field  map_qual            alignment mapping quality
 @field  query_start         segment start on query sequence (0-based, include)
 @field  query_end           segment end on query sequence (0-based, not-include)
 @field  ref_position        segment position on reference (0-based, include)
 @field  seg_type            bitwise flag representing segment type
                                1: left clip
                                2: internal insert
                                4: right clip
 @field  aln_type            bitwise flag representing read/alignment type, combination of seg_type(s)
 @field  file_offset         file offset of the alignment
 @field  aln_refstart        alignment reference start (0-based)
 @field  aln_refend          alignment reference end (0-based)
 @field  order              segment order in the same alignment
 @field  numSeg             number of segments from the same alignment
 @field  overhang           length of the shorter anchor part around segment breakpoint
 @field  match_len           matched length of the alignment (M/=/X)
 @field  read_len            read length of the alignment
 @field  aln_location_type    bitwise flag representing alignment location
                                1:  at least one end inside normal region
                                2:  both ends inside repeat/gap region
                                4:  both ends at repeat/gap boundary
                                8:  one end at repeat/gap boundary, the other inside normal region
                                16: one end at repeat/gap boundary, the other inside repeat/gap
 @field  num_te_aln     number of alignments mapped to TE
 @field  sum_query_maplen     summary query mapped-length of the TE alignments (M/I/=/X, no overlap)
 @field  sum_aln_score        summary per base alignment score of the TE alignments
 @field  sum_divergence      summary per base divergence of the TE alignments
 @field  start_idx           start index in TE alignments arr (0-based, include)
 @field  end_idx             end index in TE alignments arr (0-based, not-include)
 */
typedef struct Segment
{
    uint8_t     map_qual;
    int         query_start;
    int         query_end;
    int         ref_position;
    uint8_t     seg_type;
    uint8_t     aln_type;
    int64_t     file_offset;
    int         aln_refstart;
    int         aln_refend;
    uint8_t     order;
    uint8_t     numSeg;
    int         overhang;
    int         match_len;
    int         read_len;
    uint8_t     aln_location_type;
    uint8_t     num_te_aln;
    int         sum_query_maplen;
    float       sum_aln_score;
    float       sum_divergence;
    int         start_idx;
    int         end_idx;
} __attribute__((packed)) Segment;

/// @brief Data container used in constrcuting Segment
typedef struct SegValues
{
    uint8_t     map_qual;
    uint8_t     numSeg;
    uint8_t     seg_type;
    uint8_t     aln_type;
    int         queryPosition;
    int         ref_position;
    int         aln_refstart;
    int         match_len;
    int         read_len;
    int64_t     file_offset;
} SegValues;

/*! @brief Data container for TE alignment.
 @field  seg_idx         index of corresponding segment in segments arr
 @field  aln_score       alignment score
 @field  query_start     query start (original direction of segment)
 @field  query_end       query end (original direction of segment)
 @field  map_len         mapped length (M/I/D/=/X)
 @field  divergence     per-base divergence
 */
typedef struct TeAlignment
{
    int     seg_idx;
    int     aln_score;
    int     query_start;
    int     query_end;
    int     map_len;
    float   divergence;
} __attribute__((packed)) TeAlignment;


/***************************
 *** Initialize Segments ***
 ***************************/
#define LEFT_CLIP   1
#define MID_INSERT  2
#define RIGHT_CLIP  4
#define DUAL_CLIP   5

/// @brief Extract all segments from CIGAR, record in seg_arr
int fill_seg_arr(bam1_t *bam, Segment *seg_arr, int64_t file_offset, int min_seg_len);

/// @brief Init all segments from CIGAR
SegValues initSegmentsFromCigar(bam1_t *bam, Segment *seg_arr, int64_t file_offset, int min_seg_len);

/// @brief Init single segment from CIGAR
void initSegment(Segment *segment, SegValues segValues, int cigarLen);

/// @brief Set same values for segments extracted from the same alignment
void setSameSegValues(Segment *seg_arr, SegValues segValues);


/***********************
 *** Update Segments ***
 ***********************/
#define isLeftClip(segment) (((segment)->seg_type & LEFT_CLIP) != 0)
#define isMidInsert(segment) (((segment)->seg_type & MID_INSERT) != 0)
#define isRightClip(segment) (((segment)->seg_type & RIGHT_CLIP) != 0)
#define isFirstSegment(segment) ((segment)->order == 0)
#define isSingleSegment(segment) ((segment)->numSeg == 1)

/// @brief Update segment's overhang and location type
void update_segment(Segment *seg_arr, AiList *repeat_ailist, AiList *gap_ailist);

/// @brief Get length of the shorter ref-anchor part as overhang
int getOverhang(int overhang, int match_len);

/// @brief Set location type of the alignment
void setAlnLocationType(AiList *repeat_ailist, AiList *gap_ailist, Segment *segment);

/// @brief Get location type of a single point
uint8_t getPointLocationType(AiList *repeat_ailist, AiList *gap_ailist, int point);

/// @brief Get location type of a alignment
uint8_t getAlnLocationType(uint8_t startLocationType, uint8_t endLocationType);


/***************************************
 *** Update Segments By TeAlignments ***
 ***************************************/
#define isFirstTeAlign(segment) ((segment)->num_te_aln == 0)
#define isOverlap(te_alignment, prevTeAlignment) ((te_alignment)->query_start < (prevTeAlignment)->query_end)
#define isCover(te_alignment, prevTeAlignment) ((te_alignment)->query_end <= (prevTeAlignment)->query_end)

/// @brief Update segment values using Seg-To-TE alignments
void update_seg_by_te_arr(Segment *seg_arr, TeAlignment *te_arr, int teIdx);

/// @brief Compute query-map-len of a Seg-To-TE alignment
int getQueryMapLen(TeAlignment *te_alignment);

/// @brief Compute query-map-len of the non-verlap part between two alignments
int getOverlapQueryMapLen(TeAlignment *te_alignment, TeAlignment *prevTeAlignment);

/// @brief Update segment values using single Seg-To-TE alignment
void updateSegByTeAlignment(Segment *segment, TeAlignment *te_alignment, int teIdx, int queryMapLen);


/*******************************
 *** Initialize TeAlignments ***
 *******************************/
#define alingnedBits 0x3c03f
#define bam_is_invalid(bam) (((bam)->core.flag & BAM_FSECONDARY) != 0 || ((bam)->core.flag & BAM_FUNMAP) != 0)
#define isCigarAligned(cigar) ((alingnedBits >> (bam_cigar_op((cigar))<<1) & 3) & 1)

/// @brief Parsing and record a Seg-To-TE alignment
void fill_te_arr(bam1_t *bam, TeAlignment *te_arr);

/// @brief Get query map region on the original segment sequence
void getQueryPosition(int *queryStartPtr, int *queryEndPtr, bam1_t *bam);

/// @brief Check if first cigar is clip
int firstCigarIsClip(uint32_t *cigarArr);

/// @brief Check if final cigar is clip
int lastCigarIsClip(uint32_t *cigarArr, int numCigar);

/// @brief Compute map length and divergence of the alignment
void get_mapping_length_and_divergence(int *mapLenPtr, float *divergencePtr, bam1_t *bam);

/// @brief Compute divergence of the alignment
float getDivergence(bam1_t *bam, int map_len);

/// @brief Record a Seg-To-TE alignment
void initTeAlignment(TeAlignment *te_alignment, bam1_t *bam, int query_start, int query_end, int map_len, float divergence);


/********************
 *** Trim Segment ***
 ********************/
#define isOdd(number) ((number) & 1)

/// @brief Trim source_record from source_start to source_end, record in dest_record
int trim_segment(bam1_t *source_record, bam1_t *dest_record, int seg_idx, int source_start, int source_end);

/// @brief Reallocate memory for bam record
int samReallocBamData(bam1_t *bam, size_t desired);

static inline int reallocBamData(bam1_t *bam, size_t desired)
{
    if (desired <= bam->m_data) return 0;
    return samReallocBamData(bam, desired);
}

static inline uint32_t bamGetMemPolicy(bam1_t *bam)
{ return bam->mempolicy; }

static inline void bamSetMemPolicy(bam1_t *bam, uint32_t policy)
{ bam->mempolicy = policy; }

/// @brief Set values of dest_record
void setDestValues(bam1_t *dest_record, int destNameLen, int numNulls, int destDataLen);

/// @brief Set name of dest_record
uint8_t *setDestName(bam1_t *dest_record, char *destName, int destNameLen, int numNulls);

/// @brief Copy sequence from source_record to dest_record
void copySequence(bam1_t *source_record, bam1_t *dest_record, uint8_t *destDataPtr, int source_start, int destSeqLen);


/*********************
 *** Alinged Pairs ***
 *********************/

/// @brief Get aligned pairs of the read(alignment)
int get_aligned_pairs(bam1_t *read, int *query_arr, int *ref_arr);

#endif // SEG_UTILS_H