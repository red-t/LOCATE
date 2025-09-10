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
    polyA.existPolyT = 0;
    return polyA;
}

/// @brief Record single TE annotation
void initAnno(bam1_t *bam, sam_hdr_t *header, Cluster *clt, Annotation *anno, int idx, uint32_t *class_arr)
{
    int numCigar = bam->core.n_cigar;
    uint32_t *cigarArr = bam_get_cigar(bam);

    int query_start = 0;
    if (bam_cigar_op(cigarArr[0]) == BAM_CSOFT_CLIP)
        query_start = bam_cigar_oplen(cigarArr[0]);

    int query_end = bam->core.l_qseq;
    if (bam_cigar_op(cigarArr[numCigar - 1]) == BAM_CSOFT_CLIP)
        query_end -= bam_cigar_oplen(cigarArr[numCigar - 1]);

    int i, ref_len;
    for (i = ref_len = 0; i < numCigar; i++)
        if (bam_cigar_type(bam_cigar_op(cigarArr[i])) & 2)
            ref_len += bam_cigar_oplen(cigarArr[i]);

    anno->idx = idx;
    anno->cltTid = clt->tid;
    anno->cltIdx = clt->idx;
    anno->query_start = query_start;
    anno->query_end = query_end;
    anno->strand = bam_is_rev(bam);
    anno->tid = bam->core.tid;
    anno->ref_start = bam->core.pos;
    anno->ref_end = bam->core.pos + ref_len;
    anno->flag |= class_arr[anno->tid];

    if (bam_is_rev(bam)) {
        anno->query_start = bam->core.l_qseq - query_end;
        anno->query_end = bam->core.l_qseq - query_start;
    }

    int teLen = sam_hdr_tid2len(header, anno->tid);
    int truncSize = (((float)teLen * 0.08) < 100) ? (teLen * 0.08) : 100;
    if (anno->ref_start < truncSize)
        anno->flag |= CLT_5P_FULL;
    if ((teLen - anno->ref_end) < truncSize)
        anno->flag |= CLT_3P_FULL;
    
    if (anno->tid == clt->repTid)
        clt->flag |= CLT_SELF_TO_SELF;
}

/// @brief Find and record all TE annotations and polyA/polyT by parsing Ins-To-TE alignments
int fill_anno_arr(Cluster *clt, Annotation *anno_arr, uint32_t *class_arr, int idx)
{
    char input_fn[100] = {'\0'};
    sprintf(input_fn, "tmp_anno/%d_%d_InsToTE.bam", clt->tid, clt->idx);
    htsFile *input_bam = sam_open(input_fn, "rb");
    sam_hdr_t *header = sam_hdr_read(input_bam);
    bam1_t *bam = bam_init1();

    int num_anno = 0;
    PolyA polyA = initPolyA(idx);
    while (1)
    {
        int return_value = bam_read1(input_bam->fp.bgzf, bam);
        if (return_value < 0)
            break;
        if (bam_is_invalid(bam))
            continue;

        initAnno(bam, header, clt, &anno_arr[num_anno], idx, class_arr);
        if (anno_arr[num_anno].query_start < polyA.leftAnnoStart) {
            polyA.leftAnnoStart = anno_arr[num_anno].query_start;
            polyA.leftIdx = num_anno;
        }
        if (anno_arr[num_anno].query_end > polyA.rightAnnoEnd) {
            polyA.rightAnnoEnd = anno_arr[num_anno].query_end;
            polyA.rightIdx = num_anno;
        }
        num_anno++;
    }

    if (num_anno == 0)
        goto END;

    if (num_anno == 1)
        clt->flag |= CLT_SINGLE_TE;
    clt->flag |= CLT_TE_MAP;
    num_anno = annoPolyA(clt, anno_arr, num_anno, &polyA);
    outputTsdSeq(clt, &polyA, anno_arr, num_anno);
    goto END;

    END:
    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (input_bam != NULL) {sam_close(input_bam); input_bam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    return num_anno;
}

/// @brief Find and record all polyA/polyT
int annoPolyA(Cluster *clt, Annotation *anno_arr, int num_anno, PolyA *polyA)
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
        num_anno = setPolyA(flankSeq, anno_arr, clt, num_anno, polyA);
    }

    clt->insLen = faidx_seq_len64(insFa, insID);
    if ((clt->insLen - polyA->rightAnnoEnd) >= 5) {
        flankSeq = faidx_fetch_seq64(insFa, insID, polyA->rightAnnoEnd, clt->insLen-1, &seqLen);
        polyA->isA = 1;
        polyA->seqLen = seqLen;
        num_anno = setPolyA(flankSeq, anno_arr, clt, num_anno, polyA);
    }

    if (insFa != NULL) {fai_destroy(insFa); insFa=NULL;}
    if (flankSeq != NULL) {free(flankSeq); flankSeq=NULL;}
    return num_anno;
}

/// @brief Find and record single polyA/polyT
int setPolyA(char *flankSeq, Annotation *anno_arr, Cluster *clt, int num_anno, PolyA *polyA)
{
    // 1) Right-most TE is reverse OR 2) Is not retro-TE --> polyA is invalid
    if (polyA->isA && (isRevAnno(anno_arr[polyA->rightIdx]) || !isRetroTE(anno_arr[polyA->rightIdx].flag)))
        return num_anno;
    // 1) Left-most TE is forward OR 2) Is not retro-TE --> polyT is invalid
    if (!polyA->isA && (!isRevAnno(anno_arr[polyA->leftIdx]) || !isRetroTE(anno_arr[polyA->leftIdx].flag)))
        return num_anno;

    int thisLen = 0, max_len = 0;
    int thisSum = 0, maxSum = 0;
    int numOther = 0, position = 0;

    char targetChar = (polyA->isA) ? 'A' : 'T';
    char altChar = (polyA->isA) ? 'T' : 'A';
    int start = (polyA->isA) ? 0 : polyA->seqLen-1;
    int end = (polyA->isA) ? polyA->seqLen : -1;
    int step = (polyA->isA) ? 1 : -1;
    int numOrigin = num_anno;

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
            max_len = thisLen;
            maxSum = thisSum;
            position = i;
        }

        // Search from the new start
        if (thisSum < 0 || numOther > 10) {
            double fdr = (polyA->seqLen - max_len + 1) * pow(0.25, max_len);
            if (max_len < 5 || fdr > 0.01)
                goto RESET;

            addPolyA(anno_arr, num_anno++, clt, polyA, position, max_len);

            RESET:
            thisLen = max_len = 0;
            thisSum = maxSum = 0;
            numOther = 0;
        }
    }

    // The final polyA candiadte
    double fdr = (polyA->seqLen - max_len + 1) * pow(0.25, max_len);
    if (max_len >= 5 && fdr <= 0.01)
        addPolyA(anno_arr, num_anno++, clt, polyA, position, max_len);

    if (num_anno - numOrigin == 0)
        return num_anno;

    if (!polyA->isA) {
        polyA->existPolyT = 1;
        polyA->leftAnnoStart = anno_arr[num_anno - 1].query_start;
    }

    return num_anno;
}

/// @brief Add polyA candidate to anno_arr
void addPolyA(Annotation *anno_arr, int num_anno, Cluster *clt, PolyA *polyA, int position, int length)
{
    if (polyA->isA) {
        anno_arr[num_anno].query_start = position + polyA->rightAnnoEnd + 1 - length;
        anno_arr[num_anno].query_end = position + polyA->rightAnnoEnd + 1;
        anno_arr[num_anno].strand = 0;
        anno_arr[num_anno].tid = -1;
    }
    else {
        anno_arr[num_anno].query_start = position;
        anno_arr[num_anno].query_end = position + length;
        anno_arr[num_anno].strand = 1;
        anno_arr[num_anno].tid = -2;
        anno_arr[num_anno].extra = 0;
    }
    anno_arr[num_anno].idx = polyA->idx;
    anno_arr[num_anno].cltTid = clt->tid;
    anno_arr[num_anno].cltIdx = clt->idx;
}

/// @brief output tsd-containing-seq for tsd annotation
void outputTsdSeq(Cluster *clt, PolyA *polyA, Annotation *anno_arr, int num_anno)
{
    if (!isBothFlankMapped(clt->flag))
        return;

    char assembly_fn[100] = {'\0'};
    sprintf(assembly_fn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assembly_fn);
    hts_pos_t seqLen;

    // Adjust polyA->leftAnnoStart, when polyT is found（due to minimap2 bias）
    if (polyA->existPolyT)
        polyA->leftAnnoStart += 5;

    //                       polyA->leftAnnoStart
    //                                |          polyA->rightAnnoEnd
    //                                |            |
    // Ins-seq                      =================
    //
    //                         clt->leftMost     clt->rightMost
    //                              |                 |
    // Assembly     --------------- ================= ---------------
    //                    |           |            |             |
    //                 leftStart   leftEnd     rightStart     rightEnd
    int leftStart = clt->leftMost - 100;
    int leftEnd = clt->leftMost + polyA->leftAnnoStart;
    int rightStart = clt->rightMost - (clt->insLen - polyA->rightAnnoEnd);
    int rightEnd = clt->rightMost + 100;
    
    char *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), leftStart, leftEnd-1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), rightStart, rightEnd-1, &seqLen);

    if ((leftEnd > faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1))) || (rightStart < 0))
        goto END;

    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_%d_tsd.fa", clt->tid, clt->idx);
    FILE *fp = fopen(output_fn, "w");
    fprintf(fp, ">0_%i\n%s\n", polyA->leftAnnoStart, leftSeq);
    fprintf(fp, ">1_%i\n%s\n", (clt->insLen - polyA->rightAnnoEnd), rightSeq);
    fclose(fp);

    adjustAnno(anno_arr, num_anno, -polyA->leftAnnoStart);
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
void adjustAnno(Annotation *anno_arr, int num_anno, int leftDelta)
{
    for (int i = 0; i < num_anno; i++)
    {
        anno_arr[i].query_start += leftDelta;
        anno_arr[i].query_end += leftDelta;
    }
}

/// @brief Annotate TSD and refine breakpoint by parsing Tsd-To-Local alignments
void annotate_tsd(Cluster *clt, Annotation *anno_arr, int num_anno)
{
    char input_fn[100] = {'\0'};
    sprintf(input_fn, "tmp_anno/%d_%d_TsdToLocal.bam", clt->tid, clt->idx);
    htsFile *input_bam = sam_open(input_fn, "rb");
    sam_hdr_t *header = sam_hdr_read(input_bam);
    bam1_t *bam = bam_init1();

    int leftEnd = -1, rightStart = -1;
    int leftDelta = 0, rightDelta = 0;
    int orgLeftDelta = 0, orgRightDelta = 0;
    int mid = (clt->ref_start + clt->ref_end) / 2;
    int localStart = atoi(sam_hdr_tid2name(header, 0));

    while (1)
    {
        int return_value = bam_read1(input_bam->fp.bgzf, bam);
        if (return_value < 0)
            break;

        int readId, delta;
        sscanf((char *)bam->data, "%d_%d", &readId, &delta);
        if (readId == 0)
            orgLeftDelta = delta;
        else
            orgRightDelta = delta;

        if (bam_is_invalid(bam) || bam_is_rev(bam))
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

    int return_value = setTsd(clt, localStart, leftEnd, rightStart);
    leftDelta = (return_value == 1 || return_value == 3) ? leftDelta : orgLeftDelta;
    rightDelta = (return_value == 2 || return_value == 3) ? rightDelta : orgRightDelta;
    adjustAnno(anno_arr, num_anno, leftDelta);
    clt->leftMost -= leftDelta;
    clt->rightMost += rightDelta;
    clt->insLen += leftDelta + rightDelta;

    // query_start of left-most polyT maye be lower than 0
    if (anno_arr[num_anno-1].tid == -2 && anno_arr[num_anno-1].query_start < 0) {
        anno_arr[num_anno-1].extra = -anno_arr[num_anno-1].query_start;
        anno_arr[num_anno-1].query_start = 0;
    }

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (input_bam != NULL) {sam_close(input_bam); input_bam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
}

/// @brief Find TSD and refine breakpoint
int setTsd(Cluster *clt, int localStart, int leftEnd, int rightStart)
{
    if (leftEnd < 0 && rightStart < 0)
        return 0;

    int tmpStart, tmpEnd, return_value;
    uint32_t hasTsd = 0;
    if (leftEnd < 0 || rightStart < 0) {
        tmpEnd = localStart;
        tmpEnd += (leftEnd < 0) ? rightStart : leftEnd;
        tmpStart = tmpEnd - 1;
        return_value = (leftEnd < 0) ? 2 : 1;
    }

    if (rightStart >= 0 && leftEnd >= 0) {
        tmpStart = tmpEnd = localStart;
        tmpStart += (rightStart < leftEnd) ? rightStart : leftEnd;
        tmpEnd += (rightStart < leftEnd) ? leftEnd : rightStart;
        hasTsd = ((rightStart < leftEnd) && ((leftEnd - rightStart) < 50)) ? CLT_TSD : 0;
        return_value = 3;
    }

    clt->flag |= hasTsd;
    clt->ref_start = tmpStart;
    clt->ref_end = tmpEnd;
    return return_value;
}

/// @brief Set ins-seq structure based on annotations
void set_ins_structure(Cluster *clt, Annotation *anno_arr, int num_anno, uint32_t *class_arr, int *size_arr, int *ltr_arr)
{
    if (num_anno == 0)
        return;

    qsort(anno_arr, num_anno, sizeof(Annotation), compare);
    checkGap(anno_arr, num_anno, clt, size_arr);
    checkPolyA(anno_arr, num_anno, clt);
    checkEnd(anno_arr, num_anno, clt);
    checkTEClass(anno_arr, num_anno, clt, class_arr);
    checkFlankPolyA(anno_arr, num_anno, clt);
    checkSoloLtr(anno_arr, num_anno, clt, size_arr, ltr_arr);

    if (isRightFlankMapped(clt->flag))
        adjustAnno(anno_arr, num_anno, 50);
}

/// @brief Compare function for sorting annotations
int compare(const void *a, const void *b)
{
    Annotation *pa = (Annotation *)a;
    Annotation *pb = (Annotation *)b;
    if (pa->query_start != pb->query_start)
        return pa->query_start - pb->query_start;
    else
        return pa->query_end - pb->query_end;
}

/// @brief Check whether left-/right-most annotation is close to insSeq end
void checkGap(Annotation *anno_arr, int num_anno, Cluster *clt, int *size_arr)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(anno_arr, num_anno);
    int rightIdx = getRightIdx(anno_arr, num_anno);

    int leftCutOff = (int)MAX(0.2 * size_arr[anno_arr[leftIdx].tid], 100);
    if (anno_arr[0].query_start < leftCutOff)
        clt->flag |= CLT_LEFT_NEAR_END;

    int rightCutOff = (int)MAX(0.2 * size_arr[anno_arr[rightIdx].tid], 100);
    if ((clt->insLen - anno_arr[num_anno-1].query_end) < rightCutOff)
        clt->flag |= CLT_RIGHT_NEAR_END;
}

/// @brief Get index of the left-most TE annotation
int getLeftIdx(Annotation *anno_arr, int num_anno)
{
    for (int i = 0; i < num_anno; i++) {
        if (anno_arr[i].tid != -2)
            return i;
    }
    return 0;
}

/// @brief Get index of the right-most TE annotation
int getRightIdx(Annotation *anno_arr, int num_anno)
{
    for (int i = num_anno-1; i >= 0; i--) {
        if (anno_arr[i].tid != -1)
            return i;
    }
    return num_anno - 1;
}

/// @brief Check whether the ins-seq contains valid polyA
void checkPolyA(Annotation *anno_arr, int num_anno, Cluster *clt)
{
    if (num_anno == 1)
        return;

    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(anno_arr, num_anno);
    int rightIdx = getRightIdx(anno_arr, num_anno);

    // Define valid polyT
    int hasPolyT = (anno_arr[0].tid == -2) ? 1 : 0;
    if (hasPolyT && hasFull3P(anno_arr[leftIdx].flag)) {
        int gap1 = anno_arr[leftIdx].query_start - anno_arr[leftIdx-1].query_end;
        int gap2 = anno_arr[0].query_start;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }

    // Define valid polyA
    int hasPolyA = (anno_arr[num_anno-1].tid == -1) ? 1 : 0;
    if (hasPolyA && hasFull3P(anno_arr[rightIdx].flag)) {
        int gap1 = anno_arr[rightIdx+1].query_start - anno_arr[rightIdx].query_end;
        int gap2 = clt->insLen - anno_arr[num_anno-1].query_end;
        clt->flag |= (gap1 < 20 && gap2 < 100) ? CLT_POLYA : 0;
    }
}

/// @brief Check whether the ins-seq contains complete ends
void checkEnd(Annotation *anno_arr, int num_anno, Cluster *clt)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(anno_arr, num_anno);
    int rightIdx = getRightIdx(anno_arr, num_anno);

    if (!isRevAnno(anno_arr[leftIdx]) && hasFull5P(anno_arr[leftIdx].flag))
        clt->flag |= CLT_5P_FULL;
    if (isRevAnno(anno_arr[rightIdx]) && hasFull5P(anno_arr[rightIdx].flag))
        clt->flag |= CLT_5P_FULL;

    if (isRevAnno(anno_arr[leftIdx]) && hasFull3P(anno_arr[leftIdx].flag))
        clt->flag |= CLT_3P_FULL;
    if (!isRevAnno(anno_arr[rightIdx]) && hasFull3P(anno_arr[rightIdx].flag))
        clt->flag |= CLT_3P_FULL;

    if ((clt->flag & CLT_5P_FULL) == 0) {
        clt->flag |= (!isRevAnno(anno_arr[leftIdx]) && isRightFlankMapped(clt->flag)) ? CLT_5P_UNKNOWN : 0;
        clt->flag |= (isRevAnno(anno_arr[rightIdx]) && isLeftFlankMapped(clt->flag)) ? CLT_5P_UNKNOWN : 0;
    }

    if ((clt->flag & CLT_3P_FULL) == 0) {
        clt->flag |= (isRevAnno(anno_arr[leftIdx]) && isRightFlankMapped(clt->flag)) ? CLT_3P_UNKNOWN : 0;
        clt->flag |= (!isRevAnno(anno_arr[rightIdx]) && isLeftFlankMapped(clt->flag)) ? CLT_3P_UNKNOWN : 0;
    }
}

/// @brief Check which TE class the insertion belongs to
void checkTEClass(Annotation *anno_arr, int num_anno, Cluster *clt, uint32_t *class_arr)
{
    // Find left-/right-most TE annotation
    int leftIdx = getLeftIdx(anno_arr, num_anno);
    int rightIdx = getRightIdx(anno_arr, num_anno);

    int leftLen = anno_arr[leftIdx].ref_end - anno_arr[leftIdx].ref_start;
    int rightLen = anno_arr[rightIdx].ref_end - anno_arr[rightIdx].ref_start;
    int teTid = (leftLen > rightLen) ? anno_arr[leftIdx].tid : anno_arr[rightIdx].tid;

    clt->flag |= class_arr[teTid];
}

/// @brief Check whether left-/right- assm-flank-seq contains valid polyT/A
void checkFlankPolyA(Annotation *anno_arr, int num_anno, Cluster *clt)
{
    if (hasPolyA(clt->flag))
        return;

    char assembly_fn[100] = {'\0'};
    sprintf(assembly_fn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assembly_fn);
    hts_pos_t leftLen = 0, rightLen = 0;
    char *leftSeq = NULL, *rightSeq = NULL;

    if (isLeftFlankMapped(clt->flag) || isBothFlankMapped(clt->flag))
        leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 20), (clt->leftMost - 1), &leftLen);
    
    if (isRightFlankMapped(clt->flag) || isBothFlankMapped(clt->flag))
        rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 19), &rightLen);

    int existPolyT = searchFlankPolyA(leftSeq, 0, leftLen);
    int existPolyA = searchFlankPolyA(rightSeq, 1, rightLen);
    int leftIdx = getLeftIdx(anno_arr, num_anno);
    int rightIdx = getRightIdx(anno_arr, num_anno);

    // 1) Left-most TE is reverse AND 2) Has full-3-Prime AND 3) Is retro-TE --> T-Rich is valid
    if (existPolyT && isRevAnno(anno_arr[leftIdx]) && hasFull3P(anno_arr[leftIdx].flag) && isRetroTE(anno_arr[leftIdx].flag))
        clt->flag |= CLT_AT_RICH;

    // 1) Right-most TE is forward AND 2) Has full-3-Prime AND 3) Is retro-TE --> A-Rich is valid
    if (existPolyA && !isRevAnno(anno_arr[rightIdx]) && hasFull3P(anno_arr[rightIdx].flag) && isRetroTE(anno_arr[rightIdx].flag))
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

    int thisLen = 0, max_len = 0;
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
            max_len = thisLen;
            maxSum = thisSum;
            position = i;
        }

        if (thisSum < 0 || numOther > 10)
            thisSum = thisLen = numOther = 0;
    }

    int gapToStart = isA ? position : (seqLen - position - max_len);
    double fdr = (seqLen - max_len + 1) * pow(0.25, max_len);

    if (gapToStart > 5 || max_len < 5 || fdr > 0.01)
        return 0;

    return 1;
}

/// @brief Check whether the insertion is SOLO LTR
void checkSoloLtr(Annotation *anno_arr, int num_anno, Cluster *clt, int *size_arr, int *ltr_arr)
{
    if (!isLTR(clt->flag) || !hasSingleTE(clt->flag))
        return;

    int idx = getLeftIdx(anno_arr, num_anno);
    Annotation anno = anno_arr[idx];
    int ltrLen = ltr_arr[anno.tid], teLen = size_arr[anno.tid];

    // failed to define LTR length
    if (ltrLen == 0)
        return;

    // boundary close to left-LTR
    if (anno.ref_start <= 20 && abs(anno.ref_end - ltrLen) <= 20)
        clt->flag |= CLT_SOLO_LTR;

    // boundary close to right-LTR
    if (abs(anno.ref_start - (teLen - ltrLen)) <= 20 && abs(anno.ref_end - teLen) <= 20)
        clt->flag |= CLT_SOLO_LTR;
}

/**********************
 *** Annotation I/O ***
 **********************/

/// @brief Output formated annotation records
void output_annotations(Annotation *anno_arr, int num_anno, int start_idx, const char *te_fn)
{   
    faidx_t *teFa = fai_load(te_fn);
    int numTe = faidx_nseq(teFa), strandFlag = 0, prevIdx = anno_arr[0].idx;
    int *teTable = malloc(numTe * sizeof(int));
    memset(teTable, 0, numTe * sizeof(int));

    char queryTmp[500] = {'\0'}, refTmp[500] = {'\0'};
    char queryStr[100000] = {'\0'}, refStr[100000] = {'\0'};
    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_annoFormated.txt", start_idx);
    FILE *fp = fopen(output_fn, "w");

    for (int i = 0; i < num_anno; i++)
    {
        if (anno_arr[i].idx != prevIdx) {
            writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, anno_arr[i-1]);
            prevIdx = anno_arr[i].idx;
            strandFlag = 0;
            memset(teTable, 0, numTe * sizeof(int));
            memset(queryStr, '\0', 100000);
            memset(refStr, '\0', 100000);
        }
        formatSingleAnno(anno_arr[i], queryTmp, refTmp, teFa, teTable, &strandFlag);
        strcat(queryStr, queryTmp);
        strcat(refStr, refTmp);
    }
    // Output final anno record
    writeSingleCltAnno(strandFlag, numTe, teTable, teFa, queryStr, refStr, fp, anno_arr[num_anno-1]);
    fclose(fp);

    if (teFa != NULL) {fai_destroy(teFa); teFa = NULL;}
    if (teTable != NULL) {free(teTable); teTable = NULL;}
}

/// @brief Change single annotation record into specified format
void formatSingleAnno(Annotation anno, char *queryTmp, char *refTmp, faidx_t *teFa, int *teTable, int *strandFlag)
{
    char strand = isRevAnno(anno) ? '-' : '+';
    if (anno.tid == -2 && anno.extra > 0)
        sprintf(queryTmp, "%c:(%d)%d-%d,", strand, anno.extra, anno.query_start, anno.query_end);
    else
        sprintf(queryTmp, "%c:%d-%d,", strand, anno.query_start, anno.query_end);

    if (anno.tid == -1)
        sprintf(refTmp, "PolyA,");
    else if (anno.tid == -2)
        sprintf(refTmp, "PolyT,");
    else {
        sprintf(refTmp, "%s:%d-%d,", faidx_iseq(teFa, anno.tid), anno.ref_start, anno.ref_end);
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