#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "anno_utils.h"

/***********************************
 *** Annotate Insertion sequence ***
 ***********************************/

/// @brief Initiate PolyA
PolyA initPolyA(int idx)
{
    PolyA polyA;
    polyA.leftAnnoStart = INT_MAX;
    polyA.rightAnnoEnd = 0;
    polyA.idx = idx;
    return polyA;
}

/// @brief Record single TE annotation
void initAnno(bam1_t *bam, sam_hdr_t *header, Cluster *clt, Annotation *anno, int idx, uint32_t *classArr)
{
    int numCigar = bam->core.n_cigar;
    uint32_t *cigarArr = bam_get_cigar(bam);

    int queryStart = 0;
    if (bam_cigar_op(cigarArr[0]) == BAM_CSOFT_CLIP)
        queryStart = bam_cigar_oplen(cigarArr[0]);

    int queryEnd = bam->core.l_qseq;
    if (bam_cigar_op(cigarArr[numCigar - 1]) == BAM_CSOFT_CLIP)
        queryEnd -= bam_cigar_oplen(cigarArr[numCigar - 1]);

    int i, refLen;
    for (i = refLen = 0; i < numCigar; i++)
        if (bam_cigar_type(bam_cigar_op(cigarArr[i])) & 2)
            refLen += bam_cigar_oplen(cigarArr[i]);

    anno->idx = idx;
    anno->cltTid = clt->tid;
    anno->cltIdx = clt->idx;
    anno->queryStart = queryStart;
    anno->queryEnd = queryEnd;
    anno->strand = bam_is_rev(bam);
    anno->tid = bam->core.tid;
    anno->refStart = bam->core.pos;
    anno->refEnd = bam->core.pos + refLen;
    anno->flag |= classArr[anno->tid];

    if (bam_is_rev(bam)) {
        anno->queryStart = bam->core.l_qseq - queryEnd;
        anno->queryEnd = bam->core.l_qseq - queryStart;
    }

    int teLen = sam_hdr_tid2len(header, anno->tid);
    int truncSize = (((float)teLen * 0.08) < 100) ? (teLen * 0.08) : 100;
    if (anno->refStart < truncSize)
        anno->flag |= CLT_5P_FULL;
    if ((teLen - anno->refEnd) < truncSize)
        anno->flag |= CLT_3P_FULL;
    
    if (anno->tid == clt->repTid)
        clt->flag |= CLT_SELF_TO_SELF;
}

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fillAnnoArr(Cluster *clt, Annotation *annoArr, uint32_t *classArr, int idx)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_InsToTE.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int numAnno = 0;
    PolyA polyA = initPolyA(idx);
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam))
            continue;

        initAnno(bam, header, clt, &annoArr[numAnno], idx, classArr);
        if (annoArr[numAnno].queryStart < polyA.leftAnnoStart) {
            polyA.leftAnnoStart = annoArr[numAnno].queryStart;
            polyA.leftIdx = numAnno;
        }
        if (annoArr[numAnno].queryEnd > polyA.rightAnnoEnd) {
            polyA.rightAnnoEnd = annoArr[numAnno].queryEnd;
            polyA.rightIdx = numAnno;
        }
        numAnno++;
    }

    if (numAnno == 0)
        goto END;

    if (numAnno == 1)
        clt->flag |= CLT_SINGLE_TE;
    clt->flag |= CLT_TE_MAP;
    numAnno = annoPolyA(clt, annoArr, numAnno, &polyA);
    outputTsdSeq(clt, &polyA, annoArr, numAnno);
    goto END;

    END:
    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    return numAnno;
}

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *clt, Annotation *annoArr, int numAnno, PolyA *polyA)
{
    char insFn[100] = {'\0'};
    sprintf(insFn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);
    faidx_t *insFa = fai_load((const char *)insFn);
    const char *insID = faidx_iseq(insFa, 0);
    char *flankSeq = NULL;
    hts_pos_t seqLen;

    if (polyA->leftAnnoStart >= 5) {
        flankSeq = faidx_fetch_seq64(insFa, insID, 0, polyA->leftAnnoStart-1, &seqLen);
        polyA->isA = 0;
        polyA->seqLen = seqLen;
        numAnno = setPolyA(flankSeq, annoArr, clt, numAnno, polyA);
    }

    clt->insLen = faidx_seq_len64(insFa, insID);
    if ((clt->insLen - polyA->rightAnnoEnd) >= 5) {
        flankSeq = faidx_fetch_seq64(insFa, insID, polyA->rightAnnoEnd, clt->insLen-1, &seqLen);
        polyA->isA = 1;
        polyA->seqLen = seqLen;
        numAnno = setPolyA(flankSeq, annoArr, clt, numAnno, polyA);
    }

    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    if (flankSeq != NULL) {free(flankSeq); flankSeq=NULL;}
    return numAnno;
}

/// @brief Find and record single polyA/polyT
int setPolyA(char *flankSeq, Annotation *annoArr, Cluster *clt, int numAnno, PolyA *polyA)
{
    // 1) Right-most TE is reverse OR 2) Is not retro-TE --> polyA is invalid
    if (polyA->isA && (isRevAnno(annoArr[polyA->rightIdx]) || !isRetroTE(annoArr[polyA->rightIdx].flag)))
        return numAnno;
    // 1) Left-most TE is forward OR 2) Is not retro-TE --> polyT is invalid
    if (!polyA->isA && (!isRevAnno(annoArr[polyA->leftIdx]) || !isRetroTE(annoArr[polyA->leftIdx].flag)))
        return numAnno;

    int thisLen = 0, maxLen = 0;
    int thisSum = 0, maxSum = 0;
    int numOther = 0, position = 0;

    char targetChar = (polyA->isA) ? 'A' : 'T';
    char altChar = (polyA->isA) ? 'T' : 'A';
    int start = (polyA->isA) ? 0 : polyA->seqLen-1;
    int end = (polyA->isA) ? polyA->seqLen : -1;
    int step = (polyA->isA) ? 1 : -1;
    int numOrigin = numAnno;

    for (int i = start; i != end; i += step)
    {
        thisLen++;
        if (toupper(flankSeq[i]) == targetChar) {
            thisSum++;
        } else {
            thisSum -= toupper(flankSeq[i]) == altChar ? 1 : 3;
            numOther++;
        }

        // Update maxSum
        if (thisSum > maxSum) {
            maxLen = thisLen;
            maxSum = thisSum;
            position = i;
        }

        // Search from the new start
        if (thisSum < 0 || numOther > 10) {
            double fdr = (polyA->seqLen - maxLen + 1) * pow(0.25, maxLen);
            if (maxLen < 5 || fdr > 0.01)
                goto RESET;

            addPolyA(annoArr, numAnno++, clt, polyA, position, maxLen);

            RESET:
            thisLen = maxLen = 0;
            thisSum = maxSum = 0;
            numOther = 0;
        }
    }

    // The final polyA candiadte
    double fdr = (polyA->seqLen - maxLen + 1) * pow(0.25, maxLen);
    if (maxLen >= 5 && fdr <= 0.01)
        addPolyA(annoArr, numAnno++, clt, polyA, position, maxLen);

    if (numAnno - numOrigin == 0)
        return numAnno;

    polyA->rightAnnoEnd = polyA->isA ? annoArr[numAnno-1].queryEnd : polyA->rightAnnoEnd;
    polyA->leftAnnoStart = polyA->isA ? polyA->leftAnnoStart : annoArr[numAnno-1].queryStart;
    return numAnno;
}

/// @brief Add polyA candidate to annoArr
void addPolyA(Annotation *annoArr, int numAnno, Cluster *clt, PolyA *polyA, int position, int length)
{
    if (polyA->isA) {
        annoArr[numAnno].queryStart = position + polyA->rightAnnoEnd + 1 - length;
        annoArr[numAnno].queryEnd = position + polyA->rightAnnoEnd + 1;
        annoArr[numAnno].strand = 0;
        annoArr[numAnno].tid = -1;
    }
    else {
        annoArr[numAnno].queryStart = position;
        annoArr[numAnno].queryEnd = position + length;
        annoArr[numAnno].strand = 1;
        annoArr[numAnno].tid = -2;
    }
    annoArr[numAnno].idx = polyA->idx;
    annoArr[numAnno].cltTid = clt->tid;
    annoArr[numAnno].cltIdx = clt->idx;
}

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *clt, PolyA *polyA, Annotation *annoArr, int numAnno)
{
    if (!isBothFlankMapped(clt->flag))
        return;

    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);
    hts_pos_t seqLen;

    int leftStart = clt->leftMost - 100;
    int leftEnd = clt->leftMost + polyA->leftAnnoStart;
    int rightStart = clt->rightMost - (clt->insLen - polyA->rightAnnoEnd);
    int rightEnd = clt->rightMost + 100;
    char *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), leftStart, leftEnd-1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), rightStart, rightEnd-1, &seqLen);

    if ((leftEnd > faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1))) || (rightStart < 0))
        goto END;

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_%d_tsd.fa", clt->tid, clt->idx);
    FILE *fp = fopen(outFn, "w");
    fprintf(fp, ">0_%i\n%s\n", polyA->leftAnnoStart, leftSeq);
    fprintf(fp, ">1_%i\n%s\n", (clt->insLen - polyA->rightAnnoEnd), rightSeq);
    fclose(fp);

    adjustAnno(annoArr, numAnno, -polyA->leftAnnoStart);
    clt->leftMost = leftEnd;
    clt->rightMost = rightStart;
    clt->insLen = clt->rightMost - clt->leftMost;
    goto END;

    END:
    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}

/// @brief Adjust annotation position
void adjustAnno(Annotation *annoArr, int numAnno, int leftDelta)
{
    for (int i = 0; i < numAnno; i++)
    {
        annoArr[i].queryStart += leftDelta;
        annoArr[i].queryEnd += leftDelta;
    }
}

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annoTsd(Cluster *clt, Annotation *annoArr, int numAnno)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_%d_TsdToLocal.bam", clt->tid, clt->idx);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int leftEnd = -1, rightStart = -1;
    int leftDelta = 0, rightDelta = 0;
    int orgLeftDelta = 0, orgRightDelta = 0;
    int mid = (clt->refStart + clt->refEnd) / 2;
    int localStart = atoi(sam_hdr_tid2name(header, 0));

    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;

        int readId, delta;
        sscanf((char *)bam->data, "%d_%d", &readId, &delta);
        if (readId == 0)
            orgLeftDelta = delta;
        else
            orgRightDelta = delta;

        if (bamIsInvalid(bam) || bam_is_rev(bam))
            continue;
        
        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArr = bam_get_cigar(bam);
        
        if (readId == 0) {
            if (abs(localStart + bam_endpos(bam) - mid) > 50)
                continue;
            if (isClipInFlank(cigarArr[0], 20))
                continue;
            if (isClipInFlank(cigarArr[numCigar-1], 0))
                leftDelta = bam_cigar_oplen(cigarArr[numCigar-1]);
            leftEnd = bam_endpos(bam);
        } else {
            if (abs(localStart + bam->core.pos - mid) > 50)
                continue;
            if (isClipInFlank(cigarArr[numCigar-1], 20))
                continue;
            if (isClipInFlank(cigarArr[0], 0))
                rightDelta = bam_cigar_oplen(cigarArr[0]);
            rightStart = bam->core.pos;
        }
    }

    int retValue = setTsd(clt, localStart, leftEnd, rightStart);
    leftDelta = (retValue == 1 || retValue == 3) ? leftDelta : orgLeftDelta;
    rightDelta = (retValue == 2 || retValue == 3) ? rightDelta : orgRightDelta;
    adjustAnno(annoArr, numAnno, leftDelta);
    clt->leftMost -= leftDelta;
    clt->rightMost += rightDelta;
    clt->insLen += leftDelta + rightDelta;

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
}

/// @brief Find TSD and refine breakpoint
int setTsd(Cluster *clt, int localStart, int leftEnd, int rightStart)
{
    if (leftEnd < 0 && rightStart < 0)
        return 0;

    int tmpStart, tmpEnd, retValue;
    uint32_t hasTsd = 0;
    if (leftEnd < 0 || rightStart < 0) {
        tmpEnd = localStart;
        tmpEnd += (leftEnd < 0) ? rightStart : leftEnd;
        tmpStart = tmpEnd - 1;
        retValue = (leftEnd < 0) ? 2 : 1;
    }

    if (rightStart >= 0 && leftEnd >= 0) {
        tmpStart = tmpEnd = localStart;
        tmpStart += (rightStart < leftEnd) ? rightStart : leftEnd;
        tmpEnd += (rightStart < leftEnd) ? leftEnd : rightStart;
        hasTsd = ((rightStart < leftEnd) && ((leftEnd - rightStart) < 50)) ? CLT_TSD : 0;
        retValue = 3;
    }

    clt->flag |= hasTsd;
    clt->refStart = tmpStart;
    clt->refEnd = tmpEnd;
    return retValue;
}

/// @brief Set ins-seq structure based on annotations
void setInsStruc(Cluster *clt, Annotation *annoArr, int numAnno, uint32_t *classArr, int *sizeArr, int *ltrArr)
{
    if (numAnno == 0)
        return;

    qsort(annoArr, numAnno, sizeof(Annotation), compare);
    checkGap(annoArr, numAnno, clt, sizeArr);
    checkPolyA(annoArr, numAnno, clt);
    checkEnd(annoArr, numAnno, clt);
    checkTEClass(annoArr, numAnno, clt, classArr);
    checkFlankPolyA(annoArr, numAnno, clt);
    checkSoloLtr(annoArr, numAnno, clt, sizeArr, ltrArr);

    if (isRightFlankMapped(clt->flag))
        adjustAnno(annoArr, numAnno, 50);
}

/// @brief Compare function for sorting annotations
int compare(const void *a, const void *b)
{
    Annotation *pa = (Annotation *)a;
    Annotation *pb = (Annotation *)b;
    if (pa->queryStart != pb->queryStart)
        return pa->queryStart - pb->queryStart;
    else
        return pa->queryEnd - pb->queryEnd;
}

/// @brief Check whether left-/right-most annotation is close to insSeq end
void checkGap(Annotation *annoArr, int numAnno, Cluster *clt, int *sizeArr)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(annoArr, numAnno);
    int rightIdx = getRightIdx(annoArr, numAnno);

    int leftCutOff = (int)MAX(0.2 * sizeArr[annoArr[leftIdx].tid], 100);
    if (annoArr[0].queryStart < leftCutOff)
        clt->flag |= CLT_LEFT_NEAR_END;

    int rightCutOff = (int)MAX(0.2 * sizeArr[annoArr[rightIdx].tid], 100);
    if ((clt->insLen - annoArr[numAnno-1].queryEnd) < rightCutOff)
        clt->flag |= CLT_RIGHT_NEAR_END;
}

/// @brief Get index of the left-most TE annotation
int getLeftIdx(Annotation *annoArr, int numAnno)
{
    for (int i = 0; i < numAnno; i++) {
        if (annoArr[i].tid != -2)
            return i;
    }
    return 0;
}

/// @brief Get index of the right-most TE annotation
int getRightIdx(Annotation *annoArr, int numAnno)
{
    for (int i = numAnno-1; i >= 0; i--) {
        if (annoArr[i].tid != -1)
            return i;
    }
    return numAnno - 1;
}

/// @brief Check whether the ins-seq contains valid polyA
void checkPolyA(Annotation *annoArr, int numAnno, Cluster *clt)
{
    if (numAnno == 1)
        return;

    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(annoArr, numAnno);
    int rightIdx = getRightIdx(annoArr, numAnno);

    // Define valid polyT
    int hasPolyT = (annoArr[0].tid == -2) ? 1 : 0;
    if (hasPolyT && hasFull3P(annoArr[leftIdx].flag)) {
        int gap1 = annoArr[leftIdx].queryStart - annoArr[leftIdx-1].queryEnd;
        int gap2 = annoArr[0].queryStart;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }

    // Define valid polyA
    int hasPolyA = (annoArr[numAnno-1].tid == -1) ? 1 : 0;
    if (hasPolyA && hasFull3P(annoArr[rightIdx].flag)) {
        int gap1 = annoArr[rightIdx+1].queryStart - annoArr[rightIdx].queryEnd;
        int gap2 = clt->insLen - annoArr[numAnno-1].queryEnd;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }
}

/// @brief Check whether the ins-seq contains complete ends
void checkEnd(Annotation *annoArr, int numAnno, Cluster *clt)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(annoArr, numAnno);
    int rightIdx = getRightIdx(annoArr, numAnno);

    if (!isRevAnno(annoArr[leftIdx]) && hasFull5P(annoArr[leftIdx].flag))
        clt->flag |= CLT_5P_FULL;
    if (isRevAnno(annoArr[rightIdx]) && hasFull5P(annoArr[rightIdx].flag))
        clt->flag |= CLT_5P_FULL;

    if (isRevAnno(annoArr[leftIdx]) && hasFull3P(annoArr[leftIdx].flag))
        clt->flag |= CLT_3P_FULL;
    if (!isRevAnno(annoArr[rightIdx]) && hasFull3P(annoArr[rightIdx].flag))
        clt->flag |= CLT_3P_FULL;

    if ((clt->flag & CLT_5P_FULL) == 0) {
        clt->flag |= (!isRevAnno(annoArr[leftIdx]) && isRightFlankMapped(clt->flag)) ? CLT_5P_UNKNOWN : 0;
        clt->flag |= (isRevAnno(annoArr[rightIdx]) && isLeftFlankMapped(clt->flag)) ? CLT_5P_UNKNOWN : 0;
    }

    if ((clt->flag & CLT_3P_FULL) == 0) {
        clt->flag |= (isRevAnno(annoArr[leftIdx]) && isRightFlankMapped(clt->flag)) ? CLT_3P_UNKNOWN : 0;
        clt->flag |= (!isRevAnno(annoArr[rightIdx]) && isLeftFlankMapped(clt->flag)) ? CLT_3P_UNKNOWN : 0;
    }
}

/// @brief Check which TE class the insertion belongs to
void checkTEClass(Annotation *annoArr, int numAnno, Cluster *clt, uint32_t *classArr)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(annoArr, numAnno);
    int rightIdx = getRightIdx(annoArr, numAnno);

    int leftLen = annoArr[leftIdx].refEnd - annoArr[leftIdx].refStart;
    int rightLen = annoArr[rightIdx].refEnd - annoArr[rightIdx].refStart;
    int teTid = (leftLen > rightLen) ? annoArr[leftIdx].tid : annoArr[rightIdx].tid;

    clt->flag |= classArr[teTid];
}

/// @brief Check whether left-/right- assm-flank-seq contains valid polyT/A
void checkFlankPolyA(Annotation *annoArr, int numAnno, Cluster *clt)
{
    if (hasPolyA(clt->flag))
        return;

    char assmFn[100] = {'\0'};
    sprintf(assmFn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assmFn);
    hts_pos_t leftLen = 0, rightLen = 0;
    char *leftSeq = NULL, *rightSeq = NULL;

    if (isLeftFlankMapped(clt->flag) || isBothFlankMapped(clt->flag))
        leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 20), (clt->leftMost - 1), &leftLen);
    
    if (isRightFlankMapped(clt->flag) || isBothFlankMapped(clt->flag))
        rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 19), &rightLen);

    int existPolyT = searchFlankPolyA(leftSeq, 0, leftLen);
    int existPolyA = searchFlankPolyA(rightSeq, 1, rightLen);
    int leftIdx = getLeftIdx(annoArr, numAnno);
    int rightIdx = getRightIdx(annoArr, numAnno);

    // 1) Left-most TE is reverse AND 2) Has full-3-Prime AND 3) Is retro-TE --> T-Rich is valid
    if (existPolyT && isRevAnno(annoArr[leftIdx]) && hasFull3P(annoArr[leftIdx].flag) && isRetroTE(annoArr[leftIdx].flag))
        clt->flag |= CLT_AT_RICH;

    // 1) Right-most TE is forward AND 2) Has full-3-Prime AND 3) Is retro-TE --> A-Rich is valid
    if (existPolyA && !isRevAnno(annoArr[rightIdx]) && hasFull3P(annoArr[rightIdx].flag) && isRetroTE(annoArr[rightIdx].flag))
        clt->flag |= CLT_AT_RICH;

    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
}

/// @brief Search polyT/polyA in left-/right- assm-flank-seq sequence
int searchFlankPolyA(char *flankSeq, int isA, int seqLen)
{
    if (seqLen < 20)
        return 0;

    int thisLen = 0, maxLen = 0;
    int thisSum = 0, maxSum = 0;
    int numOther = 0, position = 0;

    char targetChar = isA ? 'A' : 'T';
    char altChar = isA ? 'T' : 'A';
    int start = isA ? 0 : seqLen-1;
    int end = isA ? seqLen : -1;
    int step = isA ? 1 : -1;

    for (int i = start; i != end; i += step)
    {
        thisLen++;
        if (toupper(flankSeq[i]) == targetChar) {
            thisSum++;
        } else {
            thisSum -= toupper(flankSeq[i]) == altChar ? 1 : 3;
            numOther++;
        }
        
        if (thisSum > maxSum) {
            maxLen = thisLen;
            maxSum = thisSum;
            position = i;
        }

        if (thisSum < 0 || numOther > 10)
            thisSum = thisLen = numOther = 0;
    }

    int gapToStart = isA ? position : (seqLen - position - maxLen);
    double fdr = (seqLen - maxLen + 1) * pow(0.25, maxLen);

    if (gapToStart > 5 || maxLen < 5 || fdr > 0.01)
        return 0;

    return 1;
}

/// @brief Check whether the insertion is SOLO LTR
void checkSoloLtr(Annotation *annoArr, int numAnno, Cluster *clt, int *sizeArr, int *ltrArr)
{
    if (!isLTR(clt->flag) || !hasSingleTE(clt->flag))
        return;

    int idx = getLeftIdx(annoArr, numAnno);
    Annotation anno = annoArr[idx];
    int ltrLen = ltrArr[anno.tid], teLen = sizeArr[anno.tid];

    // failed to define LTR length
    if (ltrLen == 0)
        return;

    // boundary close to left-LTR
    if (anno.refStart <= 20 && abs(anno.refEnd - ltrLen) <= 20)
        clt->flag |= CLT_SOLO_LTR;

    // boundary close to right-LTR
    if (abs(anno.refStart - (teLen - ltrLen)) <= 20 && abs(anno.refEnd - teLen) <= 20)
        clt->flag |= CLT_SOLO_LTR;
}

/**********************
 *** Annotation I/O ***
 **********************/

/// @brief Output formated annotation records
void outputAnno(Annotation *annoArr, int numAnno, int startIdx, const char *teFn)
{   
    faidx_t *teFa = fai_load(teFn);
    int numTe = faidx_nseq(teFa), strandFlag = 0, prevIdx = annoArr[0].idx;
    int *teTable = malloc(numTe * sizeof(int));
    memset(teTable, 0, numTe * sizeof(int));

    char queryTmp[100] = {'\0'}, refTmp[100] = {'\0'};
    char queryStr[1000] = {'\0'}, refStr[1000] = {'\0'};
    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/%d_annoFormated.txt", startIdx);
    FILE *fp = fopen(outFn, "w");

    for (int i = 0; i < numAnno; i++)
    {
        if (annoArr[i].idx != prevIdx) {
            writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, annoArr[i-1]);
            prevIdx = annoArr[i].idx;
            strandFlag = 0;
            memset(teTable, 0, numTe * sizeof(int));
            memset(queryStr, '\0', 1000);
            memset(refStr, '\0', 1000);
        }
        formatSingleAnno(annoArr[i], queryTmp, refTmp, teFa, teTable, &strandFlag);
        strcat(queryStr, queryTmp);
        strcat(refStr, refTmp);
    }
    // Output final anno record
    writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, annoArr[numAnno-1]);
    fclose(fp);

    if (teFa != NULL) {fai_destroy(teFa); teFa = NULL;}
    if (teTable != NULL) {free(teTable); teTable = NULL;}
}

/// @brief Change single annotation record into specified format
void formatSingleAnno(Annotation anno, char *queryTmp, char *refTmp, faidx_t *teFa, int *teTable, int *strandFlag)
{
    char strand = isRevAnno(anno) ? '-' : '+';
    sprintf(queryTmp, "%c:%d-%d,", strand, anno.queryStart, anno.queryEnd);

    if (anno.tid == -1)
        sprintf(refTmp, "PolyA:%d-%d,", anno.refStart, anno.refEnd);
    else if (anno.tid == -2)
        sprintf(refTmp, "PolyT:%d-%d,", anno.refStart, anno.refEnd);
    else
    {
        sprintf(refTmp, "%s:%d-%d,", faidx_iseq(teFa, anno.tid), anno.refStart, anno.refEnd);
        teTable[anno.tid] = 1;
        *strandFlag += isRevAnno(anno) ? -1 : 1;
    }
}

/// @brief Write annotation for a single cluster
void writeSingleCltAnno(int strandFlag, int numTe, int *teTable, faidx_t *teFa, char *queryStr, char *refStr, FILE *fp, Annotation anno)
{
    char strand = getCltStrand(strandFlag);
    char *cltClass = getCltClass(numTe, teTable, teFa);
    queryStr[strlen(queryStr) - 1] = '\0';
    refStr[strlen(refStr) - 1] = '\0';
    fprintf(fp, "%d-%d\t%c\t%s\t%s\t%s\n", anno.cltTid, anno.cltIdx, strand, cltClass, queryStr, refStr);
}

/// @brief Get cluster strand
char getCltStrand(int strandFlag)
{
    char strand = '.';
    if (strandFlag > 0)
        strand = '+';
    if (strandFlag < 0)
        strand = '-';

    return strand;
}

/// @brief Generate a string which represents cluster TE-class
char *getCltClass(int numTe, int *teTable, faidx_t *teFa)
{
    static char cltClass[500] = {'\0'};
    memset(cltClass, '\0', 500);
    for (int tid = 0; tid < numTe; tid++)
    {
        if (teTable[tid] == 0)
            continue;
        strcat(cltClass, faidx_iseq(teFa, tid));
        cltClass[strlen(cltClass)] = ',';
    }

    cltClass[strlen(cltClass) - 1] = '\0';
    return cltClass;
}