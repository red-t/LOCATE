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
 @field  mapQual            alignment mapping quality
 @field  queryStart         segment start on query sequence (0-based, include)
 @field  queryEnd           segment end on query sequence (0-based, not-include)
 @field  refPosition        segment position on reference (0-based, include)
 @field  segType            bitwise flag representing segment type
                                1: left clip
                                2: internal insert
                                4: right clip
 @field  alnType            bitwise flag representing read/alignment type, combination of segType(s)
 @field  fileOffset         file offset of the alignment
 @field  alnRefStart        alignment reference start (0-based)
 @field  alnRefEnd          alignment reference end (0-based)
 @field  order              segment order in the same alignment
 @field  numSeg             number of segments from the same alignment
 @field  overhang           length of the shorter anchor part around segment breakpoint
 @field  matchLen           matched length of the alignment (M/=/X)
 @field  readLen            read length of the alignment
 @field  alnLocationType    bitwise flag representing alignment location
                                1:  at least one end inside normal region
                                2:  both ends inside repeat/gap region
                                4:  both ends at repeat/gap boundary
                                8:  one end at repeat/gap boundary, the other inside normal region
                                16: one end at repeat/gap boundary, the other inside repeat/gap
 @field  numTeAlignment     number of alignments mapped to TE
 @field  sumQueryMapLen     summary query mapped-length of the TE alignments (M/I/=/X, no overlap)
 @field  sumAlnScore        summary per base alignment score of the TE alignments
 @field  sumDivergence      summary per base divergence of the TE alignments
 @field  startIdx           start index in TE alignments arr (0-based, include)
 @field  endIdx             end index in TE alignments arr (0-based, not-include)
 */
typedef struct Segment
{
    uint8_t     mapQual;
    int         queryStart;
    int         queryEnd;
    int         refPosition;
    uint8_t     segType;
    uint8_t     alnType;
    int64_t     fileOffset;
    int         alnRefStart;
    int         alnRefEnd;
    uint8_t     order;
    uint8_t     numSeg;
    int         overhang;
    int         matchLen;
    int         readLen;
    uint8_t     alnLocationType;
    uint8_t     numTeAlignment;
    int         sumQueryMapLen;
    float       sumAlnScore;
    float       sumDivergence;
    int         startIdx;
    int         endIdx;
} __attribute__((packed)) Segment;

/// @brief Data container used in constrcuting Segment
typedef struct SegValues
{
    uint8_t     mapQual;
    uint8_t     numSeg;
    uint8_t     segType;
    uint8_t     alnType;
    int         queryPosition;
    int         refPosition;
    int         alnRefStart;
    int         matchLen;
    int         readLen;
    int64_t     fileOffset;
} SegValues;

/*! @brief Data container for TE alignment.
 @field  segIdx         index of corresponding segment in segments arr
 @field  AlnScore       alignment score
 @field  queryStart     query start (original direction of segment)
 @field  queryEnd       query end (original direction of segment)
 @field  mapLen         mapped length (M/I/D/=/X)
 @field  divergence     per-base divergence
 */
typedef struct TeAlignment
{
    int     segIdx;
    int     AlnScore;
    int     queryStart;
    int     queryEnd;
    int     mapLen;
    float   divergence;
} __attribute__((packed)) TeAlignment;


/***************************
 *** Initialize Segments ***
 ***************************/
#define LEFT_CLIP   1
#define MID_INSERT  2
#define RIGHT_CLIP  4
#define DUAL_CLIP   5

/// @brief Extract all segments from CIGAR, record in segArr
int fillSegArr(bam1_t *bam, Segment *segArr, int64_t fileOffset, int minSegLen);

/// @brief Init all segments from CIGAR
SegValues initSegmentsFromCigar(bam1_t *bam, Segment *segArr, int64_t fileOffset, int minSegLen);

/// @brief Init single segment from CIGAR
void initSegment(Segment *segment, SegValues segValues, int cigarLen);

/// @brief Set same values for segments extracted from the same alignment
void setSameSegValues(Segment *segArr, SegValues segValues);


/***********************
 *** Update Segments ***
 ***********************/
#define isLeftClip(segment) (((segment)->segType & LEFT_CLIP) != 0)
#define isMidInsert(segment) (((segment)->segType & MID_INSERT) != 0)
#define isRightClip(segment) (((segment)->segType & RIGHT_CLIP) != 0)
#define isFirstSegment(segment) ((segment)->order == 0)
#define isSingleSegment(segment) ((segment)->numSeg == 1)

/// @brief Update segment's overhang and location type
void updateSegment(Segment *segArr, AiList *repeatAiList, AiList *gapAiList);

/// @brief Get length of the shorter ref-anchor part as overhang
int getOverhang(int overhang, int matchLen);

/// @brief Set location type of the alignment
void setAlnLocationType(AiList *repeatAiList, AiList *gapAiList, Segment *segment);

/// @brief Get location type of a single point
uint8_t getPointLocationType(AiList *repeatAiList, AiList *gapAiList, int point);

/// @brief Get location type of a alignment
uint8_t getAlnLocationType(uint8_t startLocationType, uint8_t endLocationType);


/***************************************
 *** Update Segments By TeAlignments ***
 ***************************************/
#define isFirstTeAlign(segment) ((segment)->numTeAlignment == 0)
#define isOverlap(teAlignment, prevTeAlignment) ((teAlignment)->queryStart < (prevTeAlignment)->queryEnd)
#define isCover(teAlignment, prevTeAlignment) ((teAlignment)->queryEnd <= (prevTeAlignment)->queryEnd)

/// @brief Update segment values using Seg-To-TE alignments
void updateSegByTeArr(Segment *segArr, TeAlignment *teArr, int teIdx);

/// @brief Compute query-map-len of a Seg-To-TE alignment
int getQueryMapLen(TeAlignment *teAlignment);

/// @brief Compute query-map-len of the non-verlap part between two alignments
int getOverlapQueryMapLen(TeAlignment *teAlignment, TeAlignment *prevTeAlignment);

/// @brief Update segment values using single Seg-To-TE alignment
void updateSegByTeAlignment(Segment *segment, TeAlignment *teAlignment, int teIdx, int queryMapLen);


/*******************************
 *** Initialize TeAlignments ***
 *******************************/
#define alingnedBits 0x3c03f
#define bamIsInvalid(bam) (((bam)->core.flag & BAM_FSECONDARY) != 0 || ((bam)->core.flag & BAM_FUNMAP) != 0)
#define isCigarAligned(cigar) ((alingnedBits >> (bam_cigar_op((cigar))<<1) & 3) & 1)

/// @brief Parsing and record a Seg-To-TE alignment
void fillTeArr(bam1_t *bam, TeAlignment *teArr);

/// @brief Get query map region on the original segment sequence
void getQueryPosition(int *queryStartPtr, int *queryEndPtr, bam1_t *bam);

/// @brief Check if first cigar is clip
int firstCigarIsClip(uint32_t *cigarArr);

/// @brief Check if final cigar is clip
int lastCigarIsClip(uint32_t *cigarArr, int numCigar);

/// @brief Compute map length and divergence of the alignment
void getMapLenAndDiv(int *mapLenPtr, float *divergencePtr, bam1_t *bam);

/// @brief Compute divergence of the alignment
float getDivergence(bam1_t *bam, int mapLen);

/// @brief Record a Seg-To-TE alignment
void initTeAlignment(TeAlignment *teAlignment, bam1_t *bam, int queryStart, int queryEnd, int mapLen, float divergence);


/********************
 *** Trim Segment ***
 ********************/
#define isOdd(number) ((number) & 1)

/// @brief Trim sourceRecord from sourceStart to sourceEnd, record in destRecord
int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIdx, int sourceStart, int sourceEnd);

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

/// @brief Set values of destRecord
void setDestValues(bam1_t *destRecord, int destNameLen, int numNulls, int destDataLen);

/// @brief Set name of destRecord
uint8_t *setDestName(bam1_t *destRecord, char *destName, int destNameLen, int numNulls);

/// @brief Copy sequence from sourceRecord to destRecord
void copySequence(bam1_t *sourceRecord, bam1_t *destRecord, uint8_t *destDataPtr, int sourceStart, int destSeqLen);


/*********************
 *** Alinged Pairs ***
 *********************/

/// @brief Get aligned pairs of the read(alignment)
int getAlignedPairs(bam1_t *read, int *queryArr, int *refArr);

#endif // SEG_UTILS_H