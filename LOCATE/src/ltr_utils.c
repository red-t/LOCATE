#include "ltr_utils.h"


/// @brief Define and return LTR sizes for all LTR-class TE
void defineLTR(const char *teFn, const char *teClassFn)
{
    faidx_t *teFa = fai_load(teFn);
    int numTe = faidx_nseq(teFa);
    int *ltrArr = malloc(numTe * sizeof(int));
    int *classArr = getClassArr(numTe, teClassFn);

    for (int i = 0; i < numTe; i++)
    {
        if (classArr[i] != 1) {
            ltrArr[i] = 0;
            continue;
        }

        outputSeq(i, teFa);
        mapEachOther(i);
        ltrArr[i] = getLtrLen(i);
    }

    char outFn[100] = {'\0'};
    sprintf(outFn, "tmp_anno/ltrSize.txt");
    FILE *outFp = fopen(outFn, "w");
    for (int i = 0; i < numTe; i++) {
        fprintf(outFp, "%d\n", ltrArr[i]);
    }
    fclose(outFp);

    if (teFa != NULL) {fai_destroy(teFa); teFa = NULL;}
    if (ltrArr != NULL) {free(ltrArr); ltrArr = NULL;}
    if (classArr != NULL) {free(classArr); classArr = NULL;}
}

/// @brief Load TE class information
int *getClassArr(int numTe, const char *teClassFn)
{
    FILE *classFile = fopen(teClassFn, "r");
    int *classArr = malloc(numTe * sizeof(int));

    int count = 0;
    char buffer[1024], *name, *class;

    while (fgets(buffer, 1024, classFile))
    {
        name = strtok(buffer, "\t");
        class = strtok(NULL, "\t");
        if(name == NULL)
            continue;
        
        size_t length = strlen(class);
        class[length - 1] = (class[length - 1] == '\n') ? '\0' : class[length - 1];
        
        if (strcmp(class, "LTR") == 0)
            classArr[count] = 1;
        else
            classArr[count] = 0;

        count++;
    }

    fclose(classFile);
    return classArr;
}

/// @brief Output left-/right- half of LTR-class TE separately
void outputSeq(int tid, faidx_t *teFa)
{
    hts_pos_t seqLen;
    int teLen = faidx_seq_len(teFa, faidx_iseq(teFa, tid));
    char *leftHalf = faidx_fetch_seq64(teFa, faidx_iseq(teFa, tid), 0, (teLen/2), &seqLen);
    char *rightHalf = faidx_fetch_seq64(teFa, faidx_iseq(teFa, tid), (teLen/2), (teLen-1), &seqLen);

    char leftFn[100] = {'\0'};
    sprintf(leftFn, "tmp_anno/%d_left.fa", tid);
    FILE *leftFp = fopen(leftFn, "w");
    fprintf(leftFp, ">0\n%s\n", leftHalf);
    fclose(leftFp);

    char rightFn[100] = {'\0'};
    sprintf(rightFn, "tmp_anno/%d_right.fa", tid);
    FILE *rightFp = fopen(rightFn, "w");
    fprintf(rightFp, ">1\n%s\n", rightHalf);
    fclose(rightFp);
}

/// @brief Map right-half to left-half
void mapEachOther(int tid)
{
    char cmd[200] = {'0'};
    sprintf(cmd, "minimap2 -aY -x map-ont tmp_anno/%d_left.fa tmp_anno/%d_right.fa | samtools view -bhS -o tmp_anno/%d_rightToLeft.bam - >/dev/null", tid, tid, tid);
    system(cmd);
}

/// @brief Compute LTR region length
int getLtrLen(int tid)
{
    char inputFn[100] = {'\0'};
    sprintf(inputFn, "tmp_anno/%d_rightToLeft.bam", tid);
    htsFile *inputBam = sam_open(inputFn, "rb");
    sam_hdr_t *header = sam_hdr_read(inputBam);
    bam1_t *bam = bam_init1();

    int ltrLen = 0;
    while (1)
    {
        int retValue = bam_read1(inputBam->fp.bgzf, bam);
        if (retValue < 0)
            break;
        if (bamIsInvalid(bam) || bam_is_rev(bam))
            continue;
        if (bam->core.pos > 10)
            continue;

        int numCigar = bam->core.n_cigar;
        uint32_t *cigarArr = bam_get_cigar(bam);
        if (isClipInFlank(cigarArr[numCigar - 1], 10))
            continue;

        ltrLen = bam_endpos(bam) - bam->core.pos;
    }

    if (ltrLen == 0)
        printf("Warning: failed to define LTR, TE tid = %d\n", tid);

    if (bam != NULL) {bam_destroy1(bam); bam=NULL;}
    if (inputBam != NULL) {sam_close(inputBam); inputBam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    return ltrLen;
}

/// @brief Test
int main() {
    defineLTR("/home/zhongrenhu/test/data/dm3.transposon_for_simulaTE.fa", "/home/zhongrenhu/test/data/dm3.repeat.class");
}