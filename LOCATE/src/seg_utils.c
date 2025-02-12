#include "seg_utils.h"

/***************************
 *** Initialize Segments ***
 ***************************/

/// @brief Extract all segments from CIGAR, record in segArr
int fillSegArr(bam1_t *bam, Segment *segArr, int64_t fileOffset, int minSegLen)
{
    SegValues sameSegValues = initSegmentsFromCigar(bam, segArr, fileOffset, minSegLen);
    if (sameSegValues.numSeg > 0)
        setSameSegValues(segArr, sameSegValues);
    return (int)sameSegValues.numSeg;
}

/// @brief Init all segments from CIGAR
SegValues initSegmentsFromCigar(bam1_t *bam, Segment *segArr, int64_t fileOffset, int minSegLen)
{
    int numCigar = bam->core.n_cigar;
    int lastCigarIdx = numCigar - 1;
    uint32_t *cigarArr = bam_get_cigar(bam);

    SegValues segValues;
    segValues.numSeg = 0;
    segValues.segType = 0;
    segValues.alnType = 0;
    segValues.matchLen = 0;
    segValues.queryPosition = 0;
    segValues.refPosition = bam->core.pos;

    for (int i = 0; i < numCigar; i++)
    {
        int cigarLen = bam_cigar_oplen(cigarArr[i]);
        switch (bam_cigar_op(cigarArr[i]))
        {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                segValues.queryPosition += cigarLen;
                segValues.refPosition += cigarLen;
                segValues.matchLen += cigarLen;
                break;

            case BAM_CSOFT_CLIP:
            case BAM_CINS:
                if (cigarLen < minSegLen) {
                    segValues.queryPosition += cigarLen;
                    break;
                }

                if (i == 0)
                    segValues.segType = LEFT_CLIP;
                else if (i == lastCigarIdx)
                    segValues.segType = RIGHT_CLIP;
                else
                    segValues.segType = MID_INSERT;
                
                initSegment(&segArr[segValues.numSeg], segValues, cigarLen);
                segValues.numSeg += 1;
                segValues.alnType |= segValues.segType;
                segValues.queryPosition += cigarLen;
                break;

            case BAM_CDEL:
            case BAM_CREF_SKIP:
                segValues.refPosition += cigarLen;
                break;
            default:
                break;
        }
    }

    segValues.fileOffset = fileOffset;
    segValues.mapQual = bam->core.qual;
    segValues.alnRefStart = bam->core.pos;
    segValues.readLen = bam->core.l_qseq;
    return segValues;
}

/// @brief Init single segment from CIGAR
void initSegment(Segment *segment, SegValues segValues, int cigarLen)
{
    segment->queryStart = segValues.queryPosition;
    segment->queryEnd = segValues.queryPosition + cigarLen;
    segment->refPosition = segValues.refPosition;
    segment->overhang = segValues.matchLen;
    segment->segType = segValues.segType;
    segment->order = segValues.numSeg;
}

/// @brief Set same values for segments extracted from the same alignment
void setSameSegValues(Segment *segArr, SegValues segValues)
{
    for (uint8_t i = 0; i < segValues.numSeg; i++)
    {
        Segment *segment = &segArr[i];
        segment->mapQual = segValues.mapQual;
        segment->alnType = segValues.alnType;
        segment->numSeg = segValues.numSeg;
        segment->alnRefStart = segValues.alnRefStart;
        segment->alnRefEnd = segValues.refPosition;
        segment->matchLen = segValues.matchLen;
        segment->readLen = segValues.readLen;
        segment->fileOffset = segValues.fileOffset;
    }
}


/***********************
 *** Update Segments ***
 ***********************/

/// @brief Update segment's overhang and location type
void updateSegment(Segment *segArr, AiList *repeatAiList, AiList *gapAiList)
{
    Segment *segment = &segArr[0];
    if (isLeftClip(segment))
        segment->overhang = segment->matchLen;
    else if (isMidInsert(segment))
        segment->overhang = getOverhang(segment->overhang, segment->matchLen);

    if (!isFirstSegment(segment)) return;
    setAlnLocationType(repeatAiList, gapAiList, segment);
    
    if (isSingleSegment(segment)) return;
    for (uint8_t i = 1; i < segment->numSeg; i++)
        segArr[i].alnLocationType = segment->alnLocationType;
}

/// @brief Get length of the shorter ref-anchor part as overhang
int getOverhang(int overhang, int matchLen)
{ return (matchLen - overhang) < overhang ? (matchLen - overhang) : overhang; }

/// @brief Set location type of the alignment
void setAlnLocationType(AiList *repeatAiList, AiList *gapAiList, Segment *segment)
{
    uint8_t startLocationType = getPointLocationType(repeatAiList, gapAiList, segment->alnRefStart);
    uint8_t endLocationType = getPointLocationType(repeatAiList, gapAiList, segment->alnRefEnd);
    segment->alnLocationType = getAlnLocationType(startLocationType, endLocationType);
}

/// @brief Get location type of a single point
uint8_t getPointLocationType(AiList *repeatAiList, AiList *gapAiList, int point)
{
    int numOverlap = 0;
    int minDistanceToOverlap = 0x7fffffff;
    ailistQueryPoint(repeatAiList, point, 50, &numOverlap, &minDistanceToOverlap);
    ailistQueryPoint(gapAiList, point, 50, &numOverlap, &minDistanceToOverlap);

    if (numOverlap == 0) return 1;
    if (minDistanceToOverlap < 50) return 2;
    return 4;
}

/// @brief Get location type of a alignment
uint8_t getAlnLocationType(uint8_t startLocationType, uint8_t endLocationType)
{
    switch (startLocationType | endLocationType)
    {
        case 1:
        case 5:
            return 1;   // at least one side inside normal region
        case 4:
            return 2;   // both ends inside repeat/gap region
        case 2:
            return 4;   // both ends at repeat/gap boundary
        case 3:
            return 8;   // one side at repeat/gap boundary, the other inside normal region
        case 6:
            return 16;  // one side at repeat/gap boundary, the other inside repeat/gap
        default:
            return 0;
    }
}


/***************************************
 *** Update Segments By TeAlignments ***
 ***************************************/

/// @brief Update segment values using all Seg-To-TE alignments
void updateSegByTeArr(Segment *segArr, TeAlignment *teArr, int teIdx)
{
    TeAlignment *teAlignment = &teArr[teIdx];
    int segIdx = teAlignment->segIdx;
    Segment *segment = &segArr[segIdx];

    if (isFirstTeAlign(segment)) {
        segment->startIdx = teIdx;
        int queryMapLen = getQueryMapLen(teAlignment);
        updateSegByTeAlignment(segment, teAlignment, teIdx, queryMapLen);
        return;
    }

    int prevIdx = teIdx - 1;
    TeAlignment *prevTeAlignment = &teArr[prevIdx];
    
    if (!isOverlap(teAlignment, prevTeAlignment)) {
        int queryMapLen = getQueryMapLen(teAlignment);
        updateSegByTeAlignment(segment, teAlignment, teIdx, queryMapLen);
        return;
    }

    if (isCover(teAlignment, prevTeAlignment)) {
        teAlignment->queryEnd = prevTeAlignment->queryEnd;
        return;
    }

    int queryMapLen = getOverlapQueryMapLen(teAlignment, prevTeAlignment);
    updateSegByTeAlignment(segment, teAlignment, teIdx, queryMapLen);
}

/// @brief Compute query-map-len of a Seg-To-TE alignment
int getQueryMapLen(TeAlignment *teAlignment)
{ return teAlignment->queryEnd - teAlignment->queryStart; }

/// @brief Compute query-map-len of the non-verlap part between two alignments
int getOverlapQueryMapLen(TeAlignment *teAlignment, TeAlignment *prevTeAlignment)
{ return teAlignment->queryEnd - prevTeAlignment->queryEnd; }

/// @brief Update segment values using single Seg-To-TE alignment
void updateSegByTeAlignment(Segment *segment, TeAlignment *teAlignment, int teIdx, int queryMapLen)
{
    segment->endIdx = teIdx + 1;
    segment->numTeAlignment += 1;
    segment->sumQueryMapLen += queryMapLen;
    segment->sumDivergence += teAlignment->divergence;
    segment->sumAlnScore += (float)teAlignment->AlnScore / teAlignment->mapLen;
}


/*******************************
 *** Initialize TeAlignments ***
 *******************************/

/// @brief Parsing and record a Seg-To-TE alignment
void fillTeArr(bam1_t *bam, TeAlignment *teArr)
{
    int queryStart, queryEnd;
    getQueryPosition(&queryStart, &queryEnd, bam);

    int mapLen;
    float divergence;
    getMapLenAndDiv(&mapLen, &divergence, bam);

    initTeAlignment(&teArr[0], bam, queryStart, queryEnd, mapLen, divergence);
}

/// @brief Get query map region on the original segment sequence
void getQueryPosition(int *queryStart, int *queryEnd, bam1_t *bam)
{
    int numCigar = bam->core.n_cigar;
    uint32_t *cigarArr = bam_get_cigar(bam);
    *queryStart = 0;
    *queryEnd = bam->core.l_qseq;

    if (!bam_is_rev(bam)) {
        if (firstCigarIsClip(cigarArr))
            *queryStart = bam_cigar_oplen(cigarArr[0]);
        if (lastCigarIsClip(cigarArr, numCigar))
            *queryEnd = bam->core.l_qseq - bam_cigar_oplen(cigarArr[numCigar-1]);
        return;
    }

    if (lastCigarIsClip(cigarArr, numCigar))
        *queryStart = bam_cigar_oplen(cigarArr[numCigar-1]);
    if (firstCigarIsClip(cigarArr))
        *queryEnd = bam->core.l_qseq - bam_cigar_oplen(cigarArr[0]);
}

/// @brief Check if first cigar is clip
int firstCigarIsClip(uint32_t *cigarArr)
{ return bam_cigar_op(cigarArr[0]) == BAM_CSOFT_CLIP; }

/// @brief Check if final cigar is clip
int lastCigarIsClip(uint32_t *cigarArr, int numCigar)
{ return bam_cigar_op(cigarArr[numCigar - 1]) == BAM_CSOFT_CLIP; }

/// @brief Compute map length and divergence of the alignment
void getMapLenAndDiv(int *mapLen, float *divergence, bam1_t *bam)
{
    uint32_t *cigarArr = bam_get_cigar(bam);
    int numCigar = bam->core.n_cigar;
    int i, len;

    for (i=len=0; i < numCigar; i++)
        if (isCigarAligned(cigarArr[i]))
            len += bam_cigar_oplen(cigarArr[i]);

    *mapLen = len;
    *divergence = getDivergence(bam, len);
}

/// @brief Compute divergence of the alignment
float getDivergence(bam1_t *bam, int mapLen)
{ return (float)(bam_aux2i(bam_aux_get(bam, "NM")) - bam_aux2i(bam_aux_get(bam, "nn"))) / mapLen; }

/// @brief Record a Seg-To-TE alignment
void initTeAlignment(TeAlignment *teAlignment, bam1_t *bam, int queryStart, int queryEnd, int mapLen, float divergence)
{
    teAlignment->segIdx = atoi(bam_get_qname(bam));
    teAlignment->AlnScore = bam_aux2i(bam_aux_get(bam, "AS"));
    teAlignment->queryStart = queryStart;
    teAlignment->queryEnd = queryEnd;
    teAlignment->mapLen = mapLen;
    teAlignment->divergence = divergence;
}


/********************
 *** Trim Segment ***
 ********************/

/// @brief Trim sourceRecord from sourceStart to sourceEnd, record in destRecord
int trimSegment(bam1_t *sourceRecord, bam1_t *destRecord, int segIdx, int sourceStart, int sourceEnd)
{
    // use segIdx as destName
    char destName[12];
    sprintf(destName, "%d", segIdx);
    
    int destNameLen = strlen(destName);
    int numNulls = 4 - destNameLen % 4;
    int destSeqLen = sourceEnd - sourceStart;
    int destDataLen = destNameLen + numNulls + (destSeqLen + 1)/2 + 1;

    // re-allocate the data buffer as needed.
    if (reallocBamData(destRecord, destDataLen) < 0) return -1;

    setDestValues(destRecord, destNameLen, numNulls, destDataLen);
    uint8_t *destDataPtr = setDestName(destRecord, destName, destNameLen, numNulls);
    // if sourceRecord is reverse, the sequence is reverse-complementary
    copySequence(sourceRecord, destRecord, destDataPtr, sourceStart, destSeqLen);

    return (int)destDataLen;
}

/// @brief Reallocate memory for bam record
int samReallocBamData(bam1_t *bam, size_t desired)
{
    // similar to sam_realloc_bam_data in htslib sam.c
    uint32_t newDataMem;
    uint8_t *newData;
    newDataMem = desired;
    kroundup32(newDataMem);
    if (newDataMem < desired) { errno = ENOMEM; return -1; }
    if ((bamGetMemPolicy(bam) & BAM_USER_OWNS_DATA) == 0) {
        newData = realloc(bam->data, newDataMem);
    } else {
        if ((newData = malloc(newDataMem)) != NULL) {
            if (bam->l_data > 0)
                memcpy(newData, bam->data,
                       bam->l_data < (int)bam->m_data ? bam->l_data : (int)bam->m_data);
            bamSetMemPolicy(bam, bamGetMemPolicy(bam) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!newData) return -1;
    bam->data = newData;
    bam->m_data = newDataMem;
    return 0;
}

/// @brief Set values of destRecord
void setDestValues(bam1_t *destRecord, int destNameLen, int numNulls, int destDataLen)
{
    destRecord->core.l_qname = (uint16_t)(destNameLen + numNulls);
    destRecord->core.l_extranul = (uint8_t)(numNulls - 1);
    destRecord->l_data = (int)destDataLen;
    destRecord->core.flag = BAM_FUNMAP;
}

/// @brief Set name of destRecord
uint8_t *setDestName(bam1_t *destRecord, char *destName, int destNameLen, int numNulls)
{
    uint8_t *destDataPtr = destRecord->data;
    
    memcpy(destDataPtr, destName, destNameLen);
    for (int i = 0; i < numNulls; i++)
        destDataPtr[destNameLen + i] = '\0';
    
    destDataPtr += destNameLen + numNulls;
    return destDataPtr;
}

/// @brief Copy sequence from sourceRecord to destRecord
void copySequence(bam1_t *sourceRecord, bam1_t *destRecord, uint8_t *destDataPtr, int sourceStart, int destSeqLen)
{
    uint8_t *sourceSeqPtr = bam_get_seq(sourceRecord);
    sourceSeqPtr += (sourceStart >> 1);

    // 1. when sourceStart is even, segment actually start from 1-st base of sourceSeqPtr[0].
    //    else, start from 2-nd base of sourceSeqPtr[0].
    // 2. when destSeqLen is even, memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1) will copy destSeqLen bases.
    //    else, will copy (destSeqLen+1) bases.
    if (!isOdd(sourceStart)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
        destRecord->core.l_qseq = (int)destSeqLen;
        return;
    }

    if (!isOdd(destSeqLen)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+2) >> 1);
        destRecord->core.l_qseq = (int)destSeqLen + 1;
        return;
    }

    memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
    destRecord->core.l_qseq = (int)destSeqLen + 1;
}


/*********************
 *** Alinged Pairs ***
 *********************/

/// @brief Get aligned pairs of the read(alignment)
int getAlignedPairs(bam1_t *read, int *queryArr, int *refArr)
{
    uint32_t *cigarArr = bam_get_cigar(read);
    int numCigar = read->core.n_cigar;
    int refPos = read->core.pos;
    int queryPos = 0;
    int size = 0;

    for (int i = 0; i < numCigar; i++)
    {
        int cigarLen = bam_cigar_oplen(cigarArr[i]);
        int end = refPos + cigarLen;

        switch (bam_cigar_op(cigarArr[i]))
        {
        case BAM_CSOFT_CLIP:
            queryPos += cigarLen;
            break;

        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            for (int j = refPos; j < end; j++)
            {
                queryArr[size] = queryPos++;
                refArr[size] = j;
                size++;
            }
            refPos += cigarLen;
            break;

        case BAM_CINS:
            for (int j = refPos; j < end; j++)
            {
                queryArr[size] = queryPos++;
                refArr[size] = -1;
                size++;
            }
            break;

        case BAM_CDEL:
        case BAM_CREF_SKIP:
            for (int j = refPos; j < end; j++)
            {
                queryArr[size] = -1;
                refArr[size] = j;
                size++;
            }
            refPos += cigarLen;
            break;

        default:
            break;
        }
    }

    return size;
}