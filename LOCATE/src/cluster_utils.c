#include "cluster_utils.h"

/**********************
 *** Update Cluster ***
 **********************/

/// @brief Update cluster values by segments features and background info
void update_cluster(Cluster *clt_arr, Segment *seg_arr, Args args)
{
    Cluster *clt = &clt_arr[0];
    setTeAlignedFrac(clt, seg_arr, args);
    if (clt->numSeg <= 2)
        setCltType(clt, seg_arr, args);
        
    if (!isValidCandidate(clt))
        return;

    updateBySegArr(clt, seg_arr, args);
    setCltLocationType(clt, args);
    divideByNumSeg(clt);
    divideByBgInfo(clt, args);
    setBackbgInfo(clt, args);
}


/************************
 *** Define Candidate ***
 ************************/

/// @brief Compute TE-aligned-fraction of a cluster
void setTeAlignedFrac(Cluster *clt, Segment *seg_arr, Args args)
{
    float numTeAlignedSeg = 0;
    for (int i = clt->start_idx; i < clt->end_idx; i++) {
        if (seg_arr[i].overhang < args.overhang)
            continue;
        if (seg_arr[i].num_te_aln > 0)
            numTeAlignedSeg += 1;
        clt->numSeg += 1;
    }

    clt->teAlignedFrac = numTeAlignedSeg / clt->numSeg;
}

/// @brief Set cluster's cltType, which represent the cluster is germ/soma
void setCltType(Cluster *clt, Segment *seg_arr, Args args)
{
    if (clt->numSeg == 1) {
        clt->cltType = 1;
        return;
    }

    int isFirst = 1;
    for (int i = clt->start_idx; i < clt->end_idx; i++) {
        Segment *segment = &seg_arr[i];
        if (overhang_too_short(segment, args.overhang))
            continue;
        if (isFirst) {
            bgzf_seek(args.genome_bam->fp.bgzf, segment->file_offset, SEEK_SET);
            bam_read1(args.genome_bam->fp.bgzf, args.first_bam_record);
            isFirst = 0;
            continue;
        }
        bgzf_seek(args.genome_bam->fp.bgzf, segment->file_offset, SEEK_SET);
        bam_read1(args.genome_bam->fp.bgzf, args.second_bam_record);
    }
    
    if (nameIsSame(args.first_bam_record, args.second_bam_record))
        clt->cltType = 2;
}


/**********************************
 *** Update Cluster By Segments ***
 **********************************/

/// @brief Update cluster values by all segments
void updateBySegArr(Cluster *clt, Segment *seg_arr, Args args)
{
    int numLeft = 0, numMiddle = 0, numRight = 0;

    for (int i = clt->start_idx; i < clt->end_idx; i++)
    {
        Segment *segment = &seg_arr[i];
        if (overhang_too_short(segment, args.overhang))
            continue;
        countValuesFromSeg(clt, args, segment, &numLeft, &numMiddle, &numRight);
    }

    setEntropy(clt, numLeft, numMiddle, numRight);
    setBalanceRatio(clt, numLeft, numRight);
    setNumSegType(clt);
}

/// @brief Update cluster values by single segment
void countValuesFromSeg(Cluster *clt, Args args, Segment *segment, int *numLeft, int *numMiddle, int *numRight)
{
    countDifferentSeg(numLeft, numMiddle, numRight, segment);
    clt->numSegType |= segment->seg_type;
    clt->meanMapQual += segment->map_qual;
    countAlnFracs(clt, segment);

    if (segment->map_qual < 5)
        clt->lowMapQualFrac += 1;
    if (isDualClip(segment))
        clt->dualClipFrac += 1;
    if (noTeAlignment(segment)) {
        clt->meanDivergence += args.bg_div;
        return;
    }

    clt->meanAlnScore += segment->sum_aln_score / segment->num_te_aln;
    clt->meanQueryMapFrac += (float)segment->sum_query_maplen / (segment->query_end - segment->query_start);
    clt->meanDivergence += segment->sum_divergence / segment->num_te_aln;
}

/// @brief Count the number of segments with different type
void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment)
{
    switch (segment->seg_type)
    {
    case LEFT_CLIP:
        *numLeft += 1; break;
    case RIGHT_CLIP:
        *numRight += 1; break;
    case MID_INSERT:
        *numMiddle += 1; break;
    default:
        break;
    }
}

/// @brief Count the number of segments with different aln_location_type
void countAlnFracs(Cluster *clt, Segment *segment)
{
    switch (segment->aln_location_type)
    {
    case 1:
        clt->alnFrac1 += 1; break;
    case 2:
        clt->alnFrac2 += 1; break;
    case 4:
        clt->alnFrac4 += 1; break;
    case 8:
        clt->alnFrac8 += 1; break;
    case 16:
        clt->alnFrac16 += 1; break;
    default:
        break;
    }
}

/// @brief Compute entropy based on the number of different type segments
float getEntropy(int numLeft, int numMiddle, int numRight, int numSeg)
{
    float entropy = 0;
    if (numLeft > 0)
        entropy -= ((float)numLeft / numSeg) * log2((float)numLeft / numSeg);
    if (numMiddle > 0)
        entropy -= ((float)numMiddle / numSeg) * log2((float)numMiddle / numSeg);
    if (numRight > 0)
        entropy -= ((float)numRight / numSeg) * log2((float)numRight / numSeg);
    return entropy;
}

/// @brief Set entropy of the cluster
void setEntropy(Cluster *clt, int numLeft, int numMiddle, int numRight)
{
    clt->entropy = getEntropy(numLeft, numMiddle, numRight, clt->numSeg);
    clt->numSegRaw = clt->numSeg;
    clt->numLeft = numLeft;
    clt->numMiddle = numMiddle;
    clt->numRight = numRight;
}

/// @brief Compute BalanceRatio of the cluster
void setBalanceRatio(Cluster *clt, int numLeft, int numRight)
{ clt->balanceRatio = (MIN(numLeft, numRight) + 0.01) / (MAX(numLeft, numRight) + 0.01); }

/// @brief Compute the numer of segment type of the cluster
void setNumSegType(Cluster *clt)
{ clt->numSegType = (clt->numSegType & 1) + ((clt->numSegType & 2) >> 1) + ((clt->numSegType & 4) >> 2); }


/*********************************
 *** Update By Background Info ***
 *********************************/

/// @brief Set location type of a cluster
void setCltLocationType(Cluster *clt, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = INT_MAX, repTid = -1;
    repTid = ailistQueryInterval(args.gap_ailist, clt->ref_start, clt->ref_end, 50, &numOverlap, &minDistanceToOverlap);
    repTid = ailistQueryInterval(args.repeat_ailist, clt->ref_start, clt->ref_end, 50, &numOverlap, &minDistanceToOverlap);
    clt->repTid = repTid;

    if (numOverlap == 0) {
        clt->locationType = 1;
        return;
    }

    if (minDistanceToOverlap < 50) {
        clt->locationType = 2;
        return;
    }

    clt->locationType = 4;
}

/// @brief Divide cluster values by numSeg
void divideByNumSeg(Cluster *clt)
{
    if (clt->alnFrac1 > 0)
        clt->alnFrac1 = clt->alnFrac1 / clt->numSeg;
    if (clt->alnFrac2 > 0)
        clt->alnFrac2 = clt->alnFrac2 / clt->numSeg;
    if (clt->alnFrac4 > 0)
        clt->alnFrac4 = clt->alnFrac4 / clt->numSeg;
    if (clt->alnFrac8 > 0)
        clt->alnFrac8 = clt->alnFrac8 / clt->numSeg;
    if (clt->alnFrac16 > 0)
        clt->alnFrac16 = clt->alnFrac16 / clt->numSeg;

    clt->dualClipFrac = clt->dualClipFrac / clt->numSeg;
    clt->lowMapQualFrac = clt->lowMapQualFrac / clt->numSeg;
    clt->meanQueryMapFrac = clt->meanQueryMapFrac / clt->numSeg;
    clt->meanDivergence = clt->meanDivergence / clt->numSeg;
    clt->meanMapQual = clt->meanMapQual / clt->numSeg;
    clt->meanAlnScore = clt->meanAlnScore / clt->numSeg;
}

/// @brief Divide cluster values by background info
void divideByBgInfo(Cluster *clt, Args args)
{
    clt->meanDivergence = clt->meanDivergence / args.bg_div;
    clt->numSeg = clt->numSeg / args.bg_depth;
}

/// @brief Set background info of a cluster
void setBackbgInfo(Cluster *clt, Args args)
{
    clt->bgDiv = args.bg_div;
    clt->bgDepth = args.bg_depth;
    clt->bgReadLen = args.bg_read_len;
}


/*****************
 *** Filtering ***
 *****************/

/// @brief Check if the cluster inersect with blacklist
void intersect_black_list(Cluster *clt, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.blacklist_ailist, clt->ref_start, clt->ref_end, 2, &numOverlap, &minDistanceToOverlap);
    if (numOverlap > 0)
        clt->isInBlacklist = 1;
}
