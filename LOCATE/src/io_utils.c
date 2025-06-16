#include "io_utils.h"

/****************************
 *** Segment Sequence IO  ***
 ****************************/

/// @brief select a segment to output as assembled contig
int get_ouput_segidx(Cluster *clt, Segment *seg_arr, Args args)
{
    int clipIdx = -1;
    int insIdx = -1;
    int maxReadLen = 0;
    int maxClipLen = 0;

    for (int i = clt->start_idx; i < clt->end_idx; i++)
    {
        Segment *segment = &seg_arr[i];
        if (overhang_too_short(segment, args.overhang))
            continue;

        if (isMidInsert(segment)) {
            if (segment->read_len <= maxReadLen)
                continue;

            maxReadLen = segment->read_len;
            insIdx = i;
            continue;
        }

        int clipLen = segment->query_end - segment->query_start;
        if (clipLen > maxClipLen) {
            maxClipLen = clipLen;
            clipIdx = i;
        }
    }

    if (insIdx > 0) // select longest spanning read
        return insIdx;

    return clipIdx; // or clip read with largest clipSize
}

/// @brief get extended region of the segment
void setTrimRegion(Segment *segment, int *start, int *end, int flank_size)
{
    *start = 0;
    *end = segment->read_len;

    if (segment->query_start - flank_size > 0)
        *start = segment->query_start - flank_size;
    if (segment->query_end + flank_size < segment->read_len)
        *end = segment->query_end + flank_size;
}


/*************************
 *** Flank Sequence IO ***
 *************************/

/// @brief Initiate FlankRegion
FlankRegion initFlankRegion()
{
    FlankRegion region;
    region.start1 = 0;
    region.start2 = 0;
    region.end1 = 0;
    region.end2 = 0;
    return region;
}

/// @brief output flank-seq and local-seq for all clusters
void extract_ref_flankseq(char *ref_fn, Cluster *clt_arr, int start_idx, int end_idx)
{
    FlankRegion region = initFlankRegion();
    faidx_t *refFa = fai_load((const char *)ref_fn);

    for (int i = start_idx; i < end_idx; i++)
    {
        Cluster *clt = &clt_arr[i];
        setFlankRegion(clt, &region);
        outputFlank(clt, refFa, region);
        outputLocal(clt, refFa, region);
    }

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
}

/// @brief define flank region on ref-genome
void setFlankRegion(Cluster *clt, FlankRegion *region)
{
    region->start1 = clt->ref_start - 499;
    region->end1 = clt->ref_start;
    region->start2 = clt->ref_end;
    region->end2 = clt->ref_end + 499;
    return;
}

/// @brief output flank-seq for single cluster
void outputFlank(Cluster *clt, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    const char *chrom = faidx_iseq(refFa, clt->tid);
    char *leftSeq = faidx_fetch_seq64(refFa, chrom, region.start1, region.end1, &seqLen);
    char *rightSeq = faidx_fetch_seq64(refFa, chrom, region.start2, region.end2, &seqLen);

    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_%d_flank.fa", clt->tid, clt->idx);
    FILE *fp = fopen(output_fn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);

    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
}

/// @brief output +-500bp local-seq around cluster position for tsd annotation
void outputLocal(Cluster *clt, faidx_t *refFa, FlankRegion region)
{
    hts_pos_t seqLen;
    int start = (region.start1 < 0) ? 0 : region.start1;
    const char *chrom = faidx_iseq(refFa, clt->tid);
    char *localSeq = faidx_fetch_seq64(refFa, chrom, start, region.end2, &seqLen);

    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_%d_local.fa", clt->tid, clt->idx);
    FILE *fp = fopen(output_fn, "w");
    fprintf(fp, ">%d\n%s\n", start, localSeq);
    fclose(fp);

    if (localSeq != NULL) {free(localSeq); localSeq=NULL;}
}


/*****************************
 *** Insertion Sequence IO ***
 *****************************/

/// @brief Initiate InsRegion
InsRegion initInsRegion()
{
    InsRegion region;
    region.rightMost = 0;
    region.leftMost = 0;
    region.tid1 = -1;
    region.tid2 = -1;
    region.cigar1 = 0;
    region.cigar2 = 0;
    region.flag = 0;
    return region;
}

/// @brief output insertion-seq from contig
void extract_ins_seq(Cluster *clt)
{
    InsRegion region = initInsRegion();
    setInsRegion(clt, &region);
    if (!isFlankMapped(clt->flag))
        return;

    char assembly_fn[100] = {'\0'};
    sprintf(assembly_fn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assembly_fn);

    if (isBothFlankMapped(clt->flag))
        outputInsSeq(assmFa, clt);
    else
        outputAssmFlank(assmFa, clt);

    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
    return;
}

/// @brief define insertion-seq region by Flank-To-Assm alignments
void setInsRegion(Cluster *clt, InsRegion *region)
{
    char input_fn[100] = {'\0'};
    sprintf(input_fn, "tmp_anno/%d_%d_FlankToAssm.bam", clt->tid, clt->idx);
    htsFile *input_bam = sam_open(input_fn, "rb");
    sam_hdr_t *header = sam_hdr_read(input_bam);
    bam1_t *bam = bam_init1();

    int minRightClip = INT_MAX;
    int minLeftClip = INT_MAX;
    int clipSize = 0;
    while (1)
    {
        int return_value = bam_read1(input_bam->fp.bgzf, bam);
        if (return_value < 0)
            break;
        if (bam_is_invalid(bam) || bam_is_rev(bam))
            continue;

        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArr = bam_get_cigar(bam);

        if (isLeftFlank(bam)) {
            clipSize = isClipInFlank(cigarArr[numCigar - 1], 0) ? bam_cigar_oplen(cigarArr[numCigar - 1]) : 0;
            if (clipSize > minRightClip)
                continue;
            minRightClip = clipSize;
            region->leftMost = bam_endpos(bam);
            region->tid1 = bam->core.tid;
            region->cigar1 = cigarArr[0];
        } else {
            clipSize = isClipInFlank(cigarArr[0], 0) ? bam_cigar_oplen(cigarArr[0]) : 0;
            if (clipSize > minLeftClip)
                continue;
            minLeftClip = clipSize;
            bam_cigar_oplen(cigarArr[0]);
            region->rightMost = bam->core.pos;
            region->tid2 = bam->core.tid;
            region->cigar2 = cigarArr[numCigar - 1];
        }
    }

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (input_bam != NULL) {sam_close(input_bam); input_bam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}

    adjustInsRegion(region);
    clt->tid1 = region->tid1;
    clt->leftMost = region->leftMost;
    clt->tid2 = region->tid2;
    clt->rightMost = region->rightMost;
    clt->flag |= region->flag;
}

/// @brief adjust region->flag
void adjustInsRegion(InsRegion *region)
{
    if (region->tid1 < 0 && region->tid2 < 0) // No flank mapped
        return;

    // Large clip in both flank
    if (isClipInFlank(region->cigar1, 400) && isClipInFlank(region->cigar2, 400))
        return;

    if (region->tid2 < 0) { // Only left lank mapped
        region->flag |= CLT_LEFT_FLANK_MAP;
        return;
    }

    if (region->tid1 < 0) { // Only right lank mapped
        region->flag |= CLT_RIGHT_FLANK_MAP;
        return;
    }

    if (isClipInFlank(region->cigar2, 400)) { // Large flank in right flank
        region->flag |= CLT_LEFT_FLANK_MAP;
        return;
    }

    if (isClipInFlank(region->cigar1, 400)) { // Large flank in left flank
        region->flag |= CLT_RIGHT_FLANK_MAP;
        return;
    }

    // Two flanks map to same contig
    if (region->tid1 == region->tid2 && region->leftMost <= (region->rightMost-100)) {
        region->flag |= CLT_SAME_FLANK_MAP;
        return;
    }

    if (region->tid1 != region->tid2) // Two flanks map to different contig
        region->flag |= CLT_DIFF_FLANK_MAP;
}

/// @brief output insertion-seq for annotation
void outputInsSeq(faidx_t *assmFa, Cluster *clt)
{
    char *insSeq = getInsSeq(assmFa, clt);
    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_%d_insertion.fa", clt->tid, clt->idx);

    FILE *fp = fopen(output_fn, "w");
    fprintf(fp, ">%d_%d_%d_%d\n%s\n", clt->tid1, clt->leftMost, clt->tid2, clt->rightMost, insSeq);
    fclose(fp);

    if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
}

/// @brief get insertion-seq
char *getInsSeq(faidx_t *assmFa, Cluster *clt)
{
    hts_pos_t seqLen;
    if ((clt->flag & CLT_LEFT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, INT_MAX, &seqLen);

    if ((clt->flag & CLT_RIGHT_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, (clt->rightMost-1), &seqLen);

    if ((clt->flag & CLT_SAME_FLANK_MAP) != 0)
        return faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, (clt->rightMost-1), &seqLen);

    if ((clt->flag & CLT_DIFF_FLANK_MAP) != 0) {
        char *seq1 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost, INT_MAX, &seqLen);
        char *seq2 = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, (clt->rightMost-1), &seqLen);

        int len1 = strlen(seq1), len2 = strlen(seq2);
        seqLen = len1 + len2 + 500;
        char *insSeq = (char *)malloc((seqLen + 1) * sizeof(char));
        char *temp = insSeq;
        memcpy(temp, seq1, len1); temp += len1;
        memset(temp, 'N', 500); temp += 500;
        memcpy(temp, seq2, len2);
        insSeq[seqLen] = '\0';

        if (seq1 != NULL) {free(seq1); seq1=NULL;}
        if (seq2 != NULL) {free(seq2); seq2=NULL;}
        return insSeq;
    }

    return NULL;
}

/// @brief output flank-seqs of insertion-seq from assembled-contig for re-defining insertion region
void outputAssmFlank(faidx_t *assmFa, Cluster *clt)
{
    int insLen = isLeftFlankMapped(clt->flag) ? (faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1)) - clt->leftMost) : clt->rightMost;
    if (insLen < 250) {
        outputInsSeq(assmFa, clt);
        return;
    }

    char *leftSeq = NULL, *rightSeq = NULL;
    hts_pos_t seqLen;

    if (isLeftFlankMapped(clt->flag)) {
        int assmLen = faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1));
        leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), clt->leftMost-200, clt->leftMost-1, &seqLen);
        rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), assmLen-200, assmLen-1, &seqLen);
    }

    if (isRightFlankMapped(clt->flag)) {
        leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), 0, 199, &seqLen);
        rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, clt->rightMost+199, &seqLen);
    }

    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_%d_assmFlank.fa", clt->tid, clt->idx);
    FILE *fp = fopen(output_fn, "w");
    fprintf(fp, ">0\n%s\n", leftSeq);
    fprintf(fp, ">1\n%s\n", rightSeq);
    fclose(fp);

    if (leftSeq != NULL) {free(leftSeq); leftSeq=NULL;}
    if (rightSeq != NULL) {free(rightSeq); rightSeq=NULL;}
}

/// @brief refine and ouput insertion-seq for single-flank-mapped cases
void re_extract_ins_seq(Cluster *clt)
{
    char assembly_fn[100] = {'\0'};
    sprintf(assembly_fn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assembly_fn);
    reSetInsRegion(clt, assmFa);
    outputInsSeq(assmFa, clt);

    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
    return;
}

/// @brief re-set insertion-seq region
void reSetInsRegion(Cluster *clt, faidx_t *assmFa)
{
    char input_fn[100] = {'\0'};
    sprintf(input_fn, "tmp_anno/%d_%d_AssmFlankToLocal.bam", clt->tid, clt->idx);
    htsFile *input_bam = sam_open(input_fn, "rb");
    sam_hdr_t *header = sam_hdr_read(input_bam);
    bam1_t *bam = bam_init1();
    int leftEnd = 0, rightStart = 0, leftLen = 0, rightLen = 0;
    int localStart = atoi(sam_hdr_tid2name(header, 0));
    uint32_t leftCigar1 = 0, leftCigar2 = 0;
    uint32_t rightCigar1 = 0, rightCigar2 = 0;

    while (1)
    {
        int return_value = bam_read1(input_bam->fp.bgzf, bam);
        if (return_value < 0)
            break;
        if (bam_is_invalid(bam) || bamIsSup(bam) || bam_is_rev(bam))
            continue;

        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArr = bam_get_cigar(bam);

        if (isLeftFlank(bam)) {
            leftEnd = bam_endpos(bam);
            leftLen = bam->core.l_qseq;
            leftCigar1 = cigarArr[0];
            leftCigar2 = cigarArr[numCigar - 1];
        } else {
            rightStart = bam->core.pos;
            rightLen = bam->core.l_qseq;
            rightCigar1 = cigarArr[0];
            rightCigar2 = cigarArr[numCigar - 1];
        }
    }
    
    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (input_bam != NULL) {sam_close(input_bam); input_bam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}

    if (leftEnd == 0 || rightStart == 0 )
        return;

    if (isClipInFlank(leftCigar1, 50) || isClipInFlank(rightCigar2, 50))
        return;

    int leftDelta = MIN(abs(leftEnd + localStart - clt->ref_start), abs(leftEnd + localStart - clt->ref_end));
    int rightDelta = MIN(abs(rightStart + localStart - clt->ref_start), abs(rightStart + localStart - clt->ref_end));
    if (leftDelta > 50 || rightDelta > 50)
        return;

    if (abs(rightStart - leftEnd) > 50)
        return;

    if (isLeftFlankMapped(clt->flag)) {
        clt->tid2 = clt->tid1;
        clt->rightMost = faidx_seq_len64(assmFa, faidx_iseq(assmFa, clt->tid1));
        clt->rightMost -= (rightLen - bam_cigar_oplen(rightCigar1));
    }

    if (isRightFlankMapped(clt->flag)) {
        clt->tid1 = clt->tid2;
        clt->leftMost = leftLen - bam_cigar_oplen(leftCigar2);
    }

    clt->flag &= 0xfffffff3;
    clt->flag |= CLT_SAME_FLANK_MAP;
}

/*******************
 *** Cluster I/O ***
 *******************/

/// @brief Output formated cluster records
void output_clusters(Cluster *clt_arr, int start_idx, int end_idx, const char *ref_fn, const char *te_fn)
{
    faidx_t *refFa = fai_load(ref_fn);
    faidx_t *teFa = fai_load(te_fn);
    char output_fn[100] = {'\0'};
    sprintf(output_fn, "tmp_anno/%d_cltFormated.txt", start_idx);

    int isAssembled;
    char *tsdSeq = NULL, *insSeq = NULL, *leftSeq = NULL, *rightSeq = NULL;
    FILE *fp = fopen(output_fn, "w");
    for (int i = start_idx; i < end_idx; i++)
    {
        Cluster *clt = &clt_arr[i];
        if (!isTEMapped(clt->flag))
            continue;

        isAssembled = ((clt->flag & CLT_ASSEMBLED) != 0);
        tsdSeq = fetchTsdSeq(refFa, clt);
        fetchSeqs(clt, &insSeq, &leftSeq, &rightSeq);

        fprintf(fp, "%d-%d\t%s\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%f\n",
                clt->tid, clt->idx, faidx_iseq(refFa, clt->tid), clt->ref_start, clt->ref_end,
                clt->probability, clt->numSegRaw, clt->numLeft, clt->numMiddle, clt->numRight,
                isAssembled, tsdSeq, insSeq, leftSeq, rightSeq, clt->flag, clt->frequency);

        free(tsdSeq); free(insSeq); free(leftSeq); free(rightSeq);
    }
    fclose(fp);

    if (refFa != NULL) {fai_destroy(refFa); refFa=NULL;}
    if (teFa != NULL) {fai_destroy(teFa); teFa=NULL;}
}

/// @brief Fetch TSD sequence from reference genome
char *fetchTsdSeq(faidx_t *refFa, Cluster *clt)
{
    hts_pos_t seqLen;
    char *tsdSeq = NULL;
    if (hasTSD(clt->flag))
        tsdSeq = faidx_fetch_seq64(refFa, faidx_iseq(refFa, clt->tid), clt->ref_start, clt->ref_end-1, &seqLen);
    else {
        tsdSeq = faidx_fetch_seq64(refFa, faidx_iseq(refFa, clt->tid), 0, 0, &seqLen);
        tsdSeq[0] = '.';
    }
    return tsdSeq;
}

/// @brief Fetch flank sequence from temporary file
void fetchSeqs(Cluster *clt, char **insSeq, char **leftSeq, char **rightSeq)
{
    hts_pos_t seqLen;
    char assembly_fn[100] = {'\0'};
    sprintf(assembly_fn, "tmp_assm/%d_%d_assembled.fa", clt->tid, clt->idx);
    faidx_t *assmFa = fai_load((const char *)assembly_fn);

    if (isLeftFlankMapped(clt->flag)) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        **rightSeq = '.';
        goto END;
    }

    if (isRightFlankMapped(clt->flag)) {
        *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, 0), 0, 0, &seqLen);
        *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
        **leftSeq = '.';
        goto END;
    }
    
    *leftSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid1), (clt->leftMost - 500), (clt->leftMost - 1), &seqLen);
    *rightSeq = faidx_fetch_seq64(assmFa, faidx_iseq(assmFa, clt->tid2), clt->rightMost, (clt->rightMost + 499), &seqLen);
    goto END;

    END:
    *insSeq = getAdjustSeq(assmFa, clt);
    if (assmFa != NULL) {fai_destroy(assmFa); assmFa=NULL;}
}

/// @brief get adjusted insertion-seq
char *getAdjustSeq(faidx_t *assmFa, Cluster *clt)
{
    char *insSeq = getInsSeq(assmFa, clt);
    if (isBothFlankMapped(clt->flag))
        return insSeq;

    int insLen = strlen(insSeq);
    int adjustLen = strlen(insSeq) + 50;
    char *adjustSeq = (char *)malloc((adjustLen + 1) * sizeof(char));
    char *temp = adjustSeq;

    if (isLeftFlankMapped(clt->flag)) {
        memcpy(temp, insSeq, insLen); temp += insLen;
        memset(temp, 'N', 50); temp += 50;
        adjustSeq[adjustLen] = '\0';

        if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
        return adjustSeq;
    }

    if (isRightFlankMapped(clt->flag)) {
        memset(temp, 'N', 50); temp += 50;
        memcpy(temp, insSeq, insLen); temp += insLen;
        adjustSeq[adjustLen] = '\0';

        if (insSeq != NULL) {free(insSeq); insSeq=NULL;}
        return adjustSeq;
    }

    return NULL;
}