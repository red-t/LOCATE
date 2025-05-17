#include "seg_utils.h"

/***************************
 *** Initialize Segments ***
 ***************************/

/// @brief Extract all segments from CIGAR, record in seg_arr
int fill_seg_arr(bam1_t *bam, Segment *seg_arr, int64_t file_offset, int min_seg_len)
{
    SegValues sameSegValues = initSegmentsFromCigar(bam, seg_arr, file_offset, min_seg_len);
    if (sameSegValues.numSeg > 0)
        setSameSegValues(seg_arr, sameSegValues);
    return (int)sameSegValues.numSeg;
}

/// @brief Init all segments from CIGAR
SegValues initSegmentsFromCigar(bam1_t *bam, Segment *seg_arr, int64_t file_offset, int min_seg_len)
{
    int numCigar = bam->core.n_cigar;
    int lastCigarIdx = numCigar - 1;
    uint32_t *cigarArr = bam_get_cigar(bam);

    SegValues segValues;
    segValues.numSeg = 0;
    segValues.seg_type = 0;
    segValues.aln_type = 0;
    segValues.match_len = 0;
    segValues.queryPosition = 0;
    segValues.ref_position = bam->core.pos;

    for (int i = 0; i < numCigar; i++)
    {
        int cigarLen = bam_cigar_oplen(cigarArr[i]);
        switch (bam_cigar_op(cigarArr[i]))
        {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                segValues.queryPosition += cigarLen;
                segValues.ref_position += cigarLen;
                segValues.match_len += cigarLen;
                break;

            case BAM_CSOFT_CLIP:
            case BAM_CINS:
                if (cigarLen < min_seg_len) {
                    segValues.queryPosition += cigarLen;
                    break;
                }

                if (i == 0)
                    segValues.seg_type = LEFT_CLIP;
                else if (i == lastCigarIdx)
                    segValues.seg_type = RIGHT_CLIP;
                else
                    segValues.seg_type = MID_INSERT;
                
                initSegment(&seg_arr[segValues.numSeg], segValues, cigarLen);
                segValues.numSeg += 1;
                segValues.aln_type |= segValues.seg_type;
                segValues.queryPosition += cigarLen;
                break;

            case BAM_CDEL:
            case BAM_CREF_SKIP:
                segValues.ref_position += cigarLen;
                break;
            default:
                break;
        }
    }

    segValues.file_offset = file_offset;
    segValues.map_qual = bam->core.qual;
    segValues.aln_refstart = bam->core.pos;
    segValues.read_len = bam->core.l_qseq;
    return segValues;
}

/// @brief Init single segment from CIGAR
void initSegment(Segment *segment, SegValues segValues, int cigarLen)
{
    segment->query_start = segValues.queryPosition;
    segment->query_end = segValues.queryPosition + cigarLen;
    segment->ref_position = segValues.ref_position;
    segment->overhang = segValues.match_len;
    segment->seg_type = segValues.seg_type;
    segment->order = segValues.numSeg;
}

/// @brief Set same values for segments extracted from the same alignment
void setSameSegValues(Segment *seg_arr, SegValues segValues)
{
    for (uint8_t i = 0; i < segValues.numSeg; i++)
    {
        Segment *segment = &seg_arr[i];
        segment->map_qual = segValues.map_qual;
        segment->aln_type = segValues.aln_type;
        segment->numSeg = segValues.numSeg;
        segment->aln_refstart = segValues.aln_refstart;
        segment->aln_refend = segValues.ref_position;
        segment->match_len = segValues.match_len;
        segment->read_len = segValues.read_len;
        segment->file_offset = segValues.file_offset;
    }
}


/***********************
 *** Update Segments ***
 ***********************/

/// @brief Update segment's overhang and location type
void update_segment(Segment *seg_arr, AiList *repeat_ailist, AiList *gap_ailist)
{
    Segment *segment = &seg_arr[0];
    if (isLeftClip(segment))
        segment->overhang = segment->match_len;
    else if (isMidInsert(segment))
        segment->overhang = getOverhang(segment->overhang, segment->match_len);

    if (!isFirstSegment(segment)) return;
    setAlnLocationType(repeat_ailist, gap_ailist, segment);
    
    if (isSingleSegment(segment)) return;
    for (uint8_t i = 1; i < segment->numSeg; i++)
        seg_arr[i].aln_location_type = segment->aln_location_type;
}

/// @brief Get length of the shorter ref-anchor part as overhang
int getOverhang(int overhang, int match_len)
{ return (match_len - overhang) < overhang ? (match_len - overhang) : overhang; }

/// @brief Set location type of the alignment
void setAlnLocationType(AiList *repeat_ailist, AiList *gap_ailist, Segment *segment)
{
    uint8_t startLocationType = getPointLocationType(repeat_ailist, gap_ailist, segment->aln_refstart);
    uint8_t endLocationType = getPointLocationType(repeat_ailist, gap_ailist, segment->aln_refend);
    segment->aln_location_type = getAlnLocationType(startLocationType, endLocationType);
}

/// @brief Get location type of a single point
uint8_t getPointLocationType(AiList *repeat_ailist, AiList *gap_ailist, int point)
{
    int numOverlap = 0;
    int minDistanceToOverlap = 0x7fffffff;
    ailistQueryPoint(repeat_ailist, point, 50, &numOverlap, &minDistanceToOverlap);
    ailistQueryPoint(gap_ailist, point, 50, &numOverlap, &minDistanceToOverlap);

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
void update_seg_by_te_arr(Segment *seg_arr, TeAlignment *te_arr, int teIdx)
{
    TeAlignment *te_alignment = &te_arr[teIdx];
    int seg_idx = te_alignment->seg_idx;
    Segment *segment = &seg_arr[seg_idx];

    if (isFirstTeAlign(segment)) {
        segment->start_idx = teIdx;
        int queryMapLen = getQueryMapLen(te_alignment);
        updateSegByTeAlignment(segment, te_alignment, teIdx, queryMapLen);
        return;
    }

    int prevIdx = teIdx - 1;
    TeAlignment *prevTeAlignment = &te_arr[prevIdx];
    
    if (!isOverlap(te_alignment, prevTeAlignment)) {
        int queryMapLen = getQueryMapLen(te_alignment);
        updateSegByTeAlignment(segment, te_alignment, teIdx, queryMapLen);
        return;
    }

    if (isCover(te_alignment, prevTeAlignment)) {
        te_alignment->query_end = prevTeAlignment->query_end;
        return;
    }

    int queryMapLen = getOverlapQueryMapLen(te_alignment, prevTeAlignment);
    updateSegByTeAlignment(segment, te_alignment, teIdx, queryMapLen);
}

/// @brief Compute query-map-len of a Seg-To-TE alignment
int getQueryMapLen(TeAlignment *te_alignment)
{ return te_alignment->query_end - te_alignment->query_start; }

/// @brief Compute query-map-len of the non-verlap part between two alignments
int getOverlapQueryMapLen(TeAlignment *te_alignment, TeAlignment *prevTeAlignment)
{ return te_alignment->query_end - prevTeAlignment->query_end; }

/// @brief Update segment values using single Seg-To-TE alignment
void updateSegByTeAlignment(Segment *segment, TeAlignment *te_alignment, int teIdx, int queryMapLen)
{
    segment->end_idx = teIdx + 1;
    segment->num_te_aln += 1;
    segment->sum_query_maplen += queryMapLen;
    segment->sum_divergence += te_alignment->divergence;
    segment->sum_aln_score += (float)te_alignment->aln_score / te_alignment->map_len;
}


/*******************************
 *** Initialize TeAlignments ***
 *******************************/

/// @brief Parsing and record a Seg-To-TE alignment
void fill_te_arr(bam1_t *bam, TeAlignment *te_arr)
{
    int query_start, query_end;
    getQueryPosition(&query_start, &query_end, bam);

    int map_len;
    float divergence;
    get_mapping_length_and_divergence(&map_len, &divergence, bam);

    initTeAlignment(&te_arr[0], bam, query_start, query_end, map_len, divergence);
}

/// @brief Get query map region on the original segment sequence
void getQueryPosition(int *query_start, int *query_end, bam1_t *bam)
{
    int numCigar = bam->core.n_cigar;
    uint32_t *cigarArr = bam_get_cigar(bam);
    *query_start = 0;
    *query_end = bam->core.l_qseq;

    if (!bam_is_rev(bam)) {
        if (firstCigarIsClip(cigarArr))
            *query_start = bam_cigar_oplen(cigarArr[0]);
        if (lastCigarIsClip(cigarArr, numCigar))
            *query_end = bam->core.l_qseq - bam_cigar_oplen(cigarArr[numCigar-1]);
        return;
    }

    if (lastCigarIsClip(cigarArr, numCigar))
        *query_start = bam_cigar_oplen(cigarArr[numCigar-1]);
    if (firstCigarIsClip(cigarArr))
        *query_end = bam->core.l_qseq - bam_cigar_oplen(cigarArr[0]);
}

/// @brief Check if first cigar is clip
int firstCigarIsClip(uint32_t *cigarArr)
{ return bam_cigar_op(cigarArr[0]) == BAM_CSOFT_CLIP; }

/// @brief Check if final cigar is clip
int lastCigarIsClip(uint32_t *cigarArr, int numCigar)
{ return bam_cigar_op(cigarArr[numCigar - 1]) == BAM_CSOFT_CLIP; }

/// @brief Compute map length and divergence of the alignment
void get_mapping_length_and_divergence(int *map_len, float *divergence, bam1_t *bam)
{
    uint32_t *cigarArr = bam_get_cigar(bam);
    int numCigar = bam->core.n_cigar;
    int i, len;

    for (i=len=0; i < numCigar; i++)
        if (isCigarAligned(cigarArr[i]))
            len += bam_cigar_oplen(cigarArr[i]);

    *map_len = len;
    *divergence = getDivergence(bam, len);
}

/// @brief Compute divergence of the alignment
float getDivergence(bam1_t *bam, int map_len)
{ return (float)(bam_aux2i(bam_aux_get(bam, "NM")) - bam_aux2i(bam_aux_get(bam, "nn"))) / map_len; }

/// @brief Record a Seg-To-TE alignment
void initTeAlignment(TeAlignment *te_alignment, bam1_t *bam, int query_start, int query_end, int map_len, float divergence)
{
    te_alignment->seg_idx = atoi(bam_get_qname(bam));
    te_alignment->aln_score = bam_aux2i(bam_aux_get(bam, "AS"));
    te_alignment->query_start = query_start;
    te_alignment->query_end = query_end;
    te_alignment->map_len = map_len;
    te_alignment->divergence = divergence;
}


/********************
 *** Trim Segment ***
 ********************/

/// @brief Trim source_record from source_start to source_end, record in dest_record
int trim_segment(bam1_t *source_record, bam1_t *dest_record, int seg_idx, int source_start, int source_end)
{
    // use seg_idx as destName
    char destName[12];
    sprintf(destName, "%d", seg_idx);
    
    int destNameLen = strlen(destName);
    int numNulls = 4 - destNameLen % 4;
    int destSeqLen = source_end - source_start;
    int destDataLen = destNameLen + numNulls + (destSeqLen + 1)/2 + 1;

    // re-allocate the data buffer as needed.
    if (reallocBamData(dest_record, destDataLen) < 0) return -1;

    setDestValues(dest_record, destNameLen, numNulls, destDataLen);
    uint8_t *destDataPtr = setDestName(dest_record, destName, destNameLen, numNulls);
    // if source_record is reverse, the sequence is reverse-complementary
    copySequence(source_record, dest_record, destDataPtr, source_start, destSeqLen);

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

/// @brief Set values of dest_record
void setDestValues(bam1_t *dest_record, int destNameLen, int numNulls, int destDataLen)
{
    dest_record->core.l_qname = (uint16_t)(destNameLen + numNulls);
    dest_record->core.l_extranul = (uint8_t)(numNulls - 1);
    dest_record->l_data = (int)destDataLen;
    dest_record->core.flag = BAM_FUNMAP;
}

/// @brief Set name of dest_record
uint8_t *setDestName(bam1_t *dest_record, char *destName, int destNameLen, int numNulls)
{
    uint8_t *destDataPtr = dest_record->data;
    
    memcpy(destDataPtr, destName, destNameLen);
    for (int i = 0; i < numNulls; i++)
        destDataPtr[destNameLen + i] = '\0';
    
    destDataPtr += destNameLen + numNulls;
    return destDataPtr;
}

/// @brief Copy sequence from source_record to dest_record
void copySequence(bam1_t *source_record, bam1_t *dest_record, uint8_t *destDataPtr, int source_start, int destSeqLen)
{
    uint8_t *sourceSeqPtr = bam_get_seq(source_record);
    sourceSeqPtr += (source_start >> 1);

    // 1. when source_start is even, segment actually start from 1-st base of sourceSeqPtr[0].
    //    else, start from 2-nd base of sourceSeqPtr[0].
    // 2. when destSeqLen is even, memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1) will copy destSeqLen bases.
    //    else, will copy (destSeqLen+1) bases.
    if (!isOdd(source_start)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
        dest_record->core.l_qseq = (int)destSeqLen;
        return;
    }

    if (!isOdd(destSeqLen)) {
        memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+2) >> 1);
        dest_record->core.l_qseq = (int)destSeqLen + 1;
        return;
    }

    memcpy(destDataPtr, sourceSeqPtr, (destSeqLen+1) >> 1);
    dest_record->core.l_qseq = (int)destSeqLen + 1;
}


/*********************
 *** Alinged Pairs ***
 *********************/

/// @brief Get aligned pairs of the read(alignment)
int get_aligned_pairs(bam1_t *read, int *query_arr, int *ref_arr)
{
    uint32_t *cigarArr = bam_get_cigar(read);
    int numCigar = read->core.n_cigar;
    int ref_pos = read->core.pos;
    int query_pos = 0;
    int size = 0;

    for (int i = 0; i < numCigar; i++)
    {
        int cigarLen = bam_cigar_oplen(cigarArr[i]);
        int end = ref_pos + cigarLen;

        switch (bam_cigar_op(cigarArr[i]))
        {
        case BAM_CSOFT_CLIP:
            query_pos += cigarLen;
            break;

        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            for (int j = ref_pos; j < end; j++)
            {
                query_arr[size] = query_pos++;
                ref_arr[size] = j;
                size++;
            }
            ref_pos += cigarLen;
            break;

        case BAM_CINS:
            for (int j = ref_pos; j < end; j++)
            {
                query_arr[size] = query_pos++;
                ref_arr[size] = -1;
                size++;
            }
            break;

        case BAM_CDEL:
        case BAM_CREF_SKIP:
            for (int j = ref_pos; j < end; j++)
            {
                query_arr[size] = -1;
                ref_arr[size] = j;
                size++;
            }
            ref_pos += cigarLen;
            break;

        default:
            break;
        }
    }

    return size;
}