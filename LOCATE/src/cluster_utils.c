#include "cluster_utils.h"

/**********************
 *** Update Cluster ***
 **********************/

/// @brief Update cluster values by segments features and background info
void updateCluster(Cluster *cltArr, Segment *segArr, Args args)
{
    Cluster *clt = &cltArr[0];
    setTeAlignedFrac(clt, segArr, args);
    if (clt->numSeg <= 2)
        setCltType(clt, segArr, args);
        
    if (!isValidCandidate(clt))
        return;

    updateBySegArr(clt, segArr, args);
    setCltLocationType(clt, args);
    divideByNumSeg(clt);
    divideByBgInfo(clt, args);
    setBackbgInfo(clt, args);
}


/************************
 *** Define Candidate ***
 ************************/

/// @brief Compute TE-aligned-fraction of a cluster
void setTeAlignedFrac(Cluster *clt, Segment *segArr, Args args)
{
    float numTeAlignedSeg = 0;
    for (int i = clt->startIdx; i < clt->endIdx; i++) {
        if (segArr[i].overhang < args.minOverhang)
            continue;
        if (segArr[i].numTeAlignment > 0)
            numTeAlignedSeg += 1;
        clt->numSeg += 1;
    }

    clt->teAlignedFrac = numTeAlignedSeg / clt->numSeg;
}

/// @brief Set cluster's cltType, which represent the cluster is germ/soma
void setCltType(Cluster *clt, Segment *segArr, Args args)
{
    if (clt->numSeg == 1) {
        clt->cltType = 1;
        return;
    }

    int isFirst = 1;
    for (int i = clt->startIdx; i < clt->endIdx; i++) {
        Segment *segment = &segArr[i];
        if (overhangIsShort(segment, args.minOverhang))
            continue;
        if (isFirst) {
            bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
            bam_read1(args.genomeBam->fp.bgzf, args.firstBamRecord);
            isFirst = 0;
            continue;
        }
        bgzf_seek(args.genomeBam->fp.bgzf, segment->fileOffset, SEEK_SET);
        bam_read1(args.genomeBam->fp.bgzf, args.secondBamRecord);
    }
    
    if (nameIsSame(args.firstBamRecord, args.secondBamRecord))
        clt->cltType = 2;
}


/**********************************
 *** Update Cluster By Segments ***
 **********************************/

/// @brief Update cluster values by all segments
void updateBySegArr(Cluster *clt, Segment *segArr, Args args)
{
    int numLeft = 0, numMiddle = 0, numRight = 0;

    for (int i = clt->startIdx; i < clt->endIdx; i++)
    {
        Segment *segment = &segArr[i];
        if (overhangIsShort(segment, args.minOverhang))
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
    clt->numSegType |= segment->segType;
    clt->meanMapQual += segment->mapQual;
    countAlnFracs(clt, segment);

    if (segment->mapQual < 5)
        clt->lowMapQualFrac += 1;
    if (isDualClip(segment))
        clt->dualClipFrac += 1;
    if (noTeAlignment(segment)) {
        clt->meanDivergence += args.bgDiv;
        return;
    }

    clt->meanAlnScore += segment->sumAlnScore / segment->numTeAlignment;
    clt->meanQueryMapFrac += (float)segment->sumQueryMapLen / (segment->queryEnd - segment->queryStart);
    clt->meanDivergence += segment->sumDivergence / segment->numTeAlignment;
}

/// @brief Count the number of segments with different type
void countDifferentSeg(int *numLeft, int *numMiddle, int *numRight, Segment *segment)
{
    switch (segment->segType)
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

/// @brief Count the number of segments with different alnLocationType
void countAlnFracs(Cluster *clt, Segment *segment)
{
    switch (segment->alnLocationType)
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
    repTid = ailistQueryInterval(args.gapAiList, clt->refStart, clt->refEnd, 50, &numOverlap, &minDistanceToOverlap);
    repTid = ailistQueryInterval(args.repeatAiList, clt->refStart, clt->refEnd, 50, &numOverlap, &minDistanceToOverlap);
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
    clt->meanDivergence = clt->meanDivergence / args.bgDiv;
    clt->numSeg = clt->numSeg / args.bgDepth;
}

/// @brief Set background info of a cluster
void setBackbgInfo(Cluster *clt, Args args)
{
    clt->bgDiv = args.bgDiv;
    clt->bgDepth = args.bgDepth;
    clt->bgReadLen = args.bgReadLen;
}


/*****************
 *** Filtering ***
 *****************/

/// @brief Check if the cluster inersect with blacklist
void intersectBlackList(Cluster *clt, Args args)
{
    int numOverlap = 0, minDistanceToOverlap = 0x7fffffff;
    ailistQueryInterval(args.blackAiList, clt->refStart, clt->refEnd, 2, &numOverlap, &minDistanceToOverlap);
    if (numOverlap > 0)
        clt->isInBlacklist = 1;
}
