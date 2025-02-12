#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H

#include <stdint.h>
#include <math.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "seg_utils.h"
#include "AIList.h"

/******************************
 *** Cluster related macros ***
 ******************************/
#define CLT_PASS                1
#define CLT_ASSEMBLED           2
#define CLT_LEFT_FLANK_MAP      4
#define CLT_RIGHT_FLANK_MAP     8
#define CLT_DIFF_FLANK_MAP      16
#define CLT_SAME_FLANK_MAP      32
#define CLT_TE_MAP              64
#define CLT_POLYA               128
#define CLT_TSD                 256
#define CLT_5P_FULL             512
#define CLT_3P_FULL             1024
#define CLT_5P_UNKNOWN          2048
#define CLT_3P_UNKNOWN          4096
#define CLT_SINGLE_TE           8192
#define CLT_SELF_TO_SELF        16384
#define CLT_LINE                32768
#define CLT_SINE                65536
#define CLT_RETROPOSON          131072
#define CLT_LTR                 262144
#define CLT_DNA                 524288
#define CLT_AT_RICH             1048576
#define CLT_LEFT_NEAR_END       2097152
#define CLT_RIGHT_NEAR_END      4194304
#define CLT_SOLO_LTR            8388608
#define CLT_SECONDARY           16777216
#define CLT_POLYA_ONLY          33554432

#define isTEMapped(flag) (((flag) & CLT_TE_MAP) != 0)
#define isFlankMapped(flag) (((flag) & 60) != 0)
#define isBothFlankMapped(flag) (((flag) & 48) != 0)
#define isLeftFlankMapped(flag) (((flag) & CLT_LEFT_FLANK_MAP) != 0)
#define isRightFlankMapped(flag) (((flag) & CLT_RIGHT_FLANK_MAP) != 0)
#define hasPolyA(flag) (((flag) & CLT_POLYA) != 0)
#define hasTSD(flag) (((flag) & CLT_TSD) != 0)
#define hasFull5P(flag) (((flag) & CLT_5P_FULL) != 0)
#define hasFull3P(flag) (((flag) & CLT_3P_FULL) != 0)
#define hasUnknown5P(flag) (((flag) & CLT_5P_UNKNOWN) != 0)
#define hasUnknown3P(flag) (((flag) & CLT_3P_UNKNOWN) != 0)
#define hasSingleTE(flag) (((flag) & CLT_SINGLE_TE) != 0)
#define isSelfToSelf(flag) (((flag) & CLT_SELF_TO_SELF) != 0)
#define isLINE(flag) (((flag) & CLT_LINE) != 0)
#define isSINE(flag) (((flag) & CLT_SINE) != 0)
#define isRETROPOSON(flag) (((flag) & CLT_RETROPOSON) != 0)
#define isLTR(flag) (((flag) & CLT_LTR) != 0)
#define isDNA(flag) (((flag) & CLT_DNA) != 0)
#define isRetroTE(flag) (((flag) & (CLT_LINE | CLT_SINE | CLT_RETROPOSON)) != 0)
#define isATRich(flag) (((flag) & CLT_AT_RICH) != 0)
#define isRightNearEnd(flag) (((flag) & CLT_RIGHT_NEAR_END) != 0)
#define isLeftNearEnd(flag) (((flag) & CLT_LEFT_NEAR_END) != 0)
#define isSoloLtr(flag) (((flag) & CLT_SOLO_LTR) != 0)
#define isHighFreq(clt) ((clt)->cltType == 0)


/******************
 *** Structures ***
 ******************/

/*! @brief Data container for cluster merged from segments.
 @field  tid                target id of corresponding chromosome
 @field  idx                cluster index in cluster arr (0-based)
 @field  refStart           cluster start on reference genome (0-based)
 @field  refEnd             cluster end on reference genome (0-based)
 @field  startIdx           start index in segment arr (0-based, include)
 @field  endIdx             end index in segment arr (0-based, not-include)
 @field  numSeg             number of segments in the cluster (normalized by bg depth)
 @field  cltType            cluster type
                                0: germline (multiple support reads)
                                1: somatic (1 support read & 1 alignment)
                                2: somatic (1 support read & 2 alignments)
 @field  locationType       bitwise flag representing cluster location
                                1: inside normal region
                                2: at repeat/gap boundary
                                4: inside repeat/gap
 @field  numSegType         number of different segment types
 @field  entropy            entropy based on fraction of different type segments
 @field  balanceRatio       balance ratio based on number of left- & right-clip segments
 @field  lowMapQualFrac     fraction of segments with low mapQual (<5)
 @field  dualClipFrac       fraction of "dual-clip" alignments
 @field  alnFrac1           fraction of segments with alnLocationType=1
 @field  alnFrac2           fraction of segments with alnLocationType=2
 @field  alnFrac4           fraction of segments with alnLocationType=4
 @field  alnFrac8           fraction of segments with alnLocationType=8
 @field  alnFrac16          fraction of segments with alnLocationType=16
 @field  meanMapQual        mean mapQual of cluster
 @field  meanAlnScore       mean per-base alignment score (based on teAlignments)
 @field  meanQueryMapFrac   mean query mapped fraction (based on teAlignments)
 @field  meanDivergence     mean per-base divergence ((#mismatches + #I + #D) / (#mismatches + #I + #D + #matches))
 @field  bgDiv              background divergence (for normalization)
 @field  bgDepth            background depth (for normalization)
 @field  bgReadLen          background read length
 @field  teAlignedFrac      fraction of TE-aligned segments
 @field  isInBlacklist      whether cluster intersects with blacklist
 @field  probability        the probability of the cluster to be a positive insertion
 @field  flag               bitwise flag representing cluster features
 @field  numSegRaw          number of segments in the cluster (no normalization)
 @field  numLeft            number of left-clipped segments in the cluster (no normalization)
 @field  numMiddle          number of mid-inserted segments in the cluster (no normalization)
 @field  numRight           number of right-clipped segments in the cluster (no normalization)
 @field  tid1               tid of the assembled contig which left-genome-flank-seq maps to
 @field  leftMost           left most position of the insertion seq on the assembled contig (map end of left-genome-flank-seq, included)
 @field  tid2               tid of the assembled contig which right-genome-flank-seq maps to
 @field  rightMost          right most position of the insertion seq on the assembled contig (map start of right-genome-flank-seq, not included)
 @field  insLen             length of insertion sequence
 @field  repTid             repeat tid, TE tid of the repeat which the cluster overlaps with
 */
typedef struct Cluster
{
    int         tid;
    int         idx;
    int         refStart;
    int         refEnd;
    int         startIdx;
    int         endIdx;
    float       numSeg;
    uint8_t     cltType;
    uint8_t     locationType;
    uint8_t     numSegType;
    float       entropy;
    float       balanceRatio;
    float       lowMapQualFrac;
    float       dualClipFrac;
    float       alnFrac1;
    float       alnFrac2;
    float       alnFrac4;
    float       alnFrac8;
    float       alnFrac16;
    float       meanMapQual;
    float       meanAlnScore;
    float       meanQueryMapFrac;
    float       meanDivergence;
    float       bgDiv;
    float       bgDepth;
    float       bgReadLen;
    float       teAlignedFrac;
    uint8_t     isInBlacklist;
    float       probability;
    uint32_t    flag;
    int         numSegRaw;
    int         numLeft;
    int         numMiddle;
    int         numRight;
    int         tid1;
    int         leftMost;
    int         tid2;
    int         rightMost;
    int         insLen;
    int         repTid;
} __attribute__((packed)) Cluster;

/// @brief Data container for arguments
typedef struct Args
{
    int         numThread;
    int         tid;
    int         minSegLen;
    int         maxDistance;
    int         minOverhang;
    float       bgDiv;
    float       bgDepth;
    float       bgReadLen;
    htsFile     *genomeBam;
    bam1_t      *firstBamRecord;
    bam1_t      *secondBamRecord;
    AiList      *repeatAiList;
    AiList      *gapAiList;
    AiList      *blackAiList;
} Args;


/**********************
 *** Update Cluster ***
 **********************/

/// @brief Update cluster values by segments features and background info
void updateCluster(Cluster *cltArr, Segment *segArr, Args args);


/************************
 *** Define Candidate ***
 ************************/
#define overhangIsShort(segment, minOverhang) ((segment)->overhang < (minOverhang))
#define nameIsSame(record1, record2) (strcmp(bam_get_qname((record1)), bam_get_qname((record2))) == 0)
#define isValidCandidate(clt) ((clt)->teAlignedFrac >= 0.8)

/// @brief Compute TE-aligned-fraction of a cluster
void setTeAlignedFrac(Cluster *clt, Segment *segArr, Args args);

/// @brief Set cluster's cltType, which represent the cluster is germ/soma
void setCltType(Cluster *clt, Segment *segArr, Args args);


/**********************************
 *** Update Cluster By Segments ***
 **********************************/
#define isDualClip(segment) (((segment)->alnType & DUAL_CLIP) == 5)
#define noTeAlignment(segment) ((segment)->numTeAlignment == 0)

/// @brief Update cluster values by all segments
void updateBySegArr(Cluster *clt, Segment *segArr, Args args);

/// @brief Update cluster values by single segment
void countValuesFromSeg(Cluster *clt, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight);

/// @brief Count number of different type segments
void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment);

/// @brief Count the number of segments with different alnLocationType
void countAlnFracs(Cluster *clt, Segment *segment);

/// @brief Compute entropy based on the number of different type segments
float getEntropy(int numLeft, int numMiddle, int numRight, int numSeg);

/// @brief Set entropy of the cluster
void setEntropy(Cluster *clt, int numLeft, int numMiddle, int numRight);

/// @brief Compute BalanceRatio of the cluster
void setBalanceRatio(Cluster *clt, int numLeft, int numRight);

/// @brief Compute the numer of segment type of the cluster
void setNumSegType(Cluster *clt);


/*********************************
 *** Update By Background Info ***
 *********************************/

/// @brief Set location type of a cluster
void setCltLocationType(Cluster *clt, Args args);

/// @brief Divide cluster values by numSeg
void divideByNumSeg(Cluster *clt);

/// @brief Divide cluster values by background info
void divideByBgInfo(Cluster *clt, Args args);

/// @brief Set background info of a cluster
void setBackbgInfo(Cluster *clt, Args args);


/*****************
 *** Filtering ***
 *****************/

/// @brief Check if the cluster inersect with blacklist
void intersectBlackList(Cluster *clt, Args args);

#endif // CLUSTER_UTILS_H