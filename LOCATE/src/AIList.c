#include "AIList.h"
#include "htslib/sam.h"

/***************************
 *** AiList Construction ***
 ***************************/
AiList *initAiList(void)
{
    AiList *ail = malloc(1*sizeof(AiList));
    if (!ail) { fprintf(stderr, "Error: Out of memory in initAiList\n"); exit(EXIT_FAILURE); }
    ail->contigList = malloc(1*sizeof(Contig));
    if (!ail->contigList) { fprintf(stderr, "Error: Out of memory in initAiList\n"); exit(EXIT_FAILURE); }
    Contig *firstContig = &ail->contigList[0];

    firstContig->numInterval = 0;
    firstContig->maxIntervals = 64;
    firstContig->intervalList = malloc(firstContig->maxIntervals*sizeof(Interval));
    if (!firstContig->intervalList) { fprintf(stderr, "Error: Out of memory in initAiList\n"); exit(EXIT_FAILURE); }
    return ail;
}

void destroyAiList(AiList *ail)
{
    if (ail == 0) return;
    free(ail->contigList[0].intervalList);
    free(ail->contigList[0].maxEndList);
    free(ail->contigList);
    free(ail);
}

void addInterval(AiList *ail, int start, int end, int repTid)
{
    if(start > end) return;

    Contig *firstContig = &ail->contigList[0];
    if(firstContig->numInterval == firstContig->maxIntervals)
        EXPAND(firstContig->intervalList, firstContig->maxIntervals);

    Interval *newInterval = &firstContig->intervalList[firstContig->numInterval++];
    newInterval->start = start;
    newInterval->end = end;
    newInterval->repTid = repTid;
    return;
}

void readBED(AiList *ail, const char* bed_fn, const char* targetChrom)
{
    samFile *te_bam = sam_open("tmp_build/tmp.bam", "rb");
    sam_hdr_t *header = sam_hdr_read(te_bam);
    FILE *fileHandle = fopen(bed_fn, "r");
    if (!fileHandle) { sam_hdr_destroy(header); sam_close(te_bam); return; }

    char buffer[1024];
    char *chrom, *start, *end, *name;
    while (fgets(buffer, 1024, fileHandle)) {
        chrom = strtok(buffer, "\t");
        start = strtok(NULL, "\t");
        end = strtok(NULL, "\t");
        name = strtok(NULL, "\t");

        if(chrom == NULL)
            continue;
        
        int ret = strcmp(chrom, targetChrom);
        if (ret != 0)
            continue;

        size_t length = strlen(name);
        name[length-1] = (name[length-1] == '\n') ? '\0' : name[length-1];
        addInterval(ail, atol(start), atol(end), sam_hdr_name2tid(header, name));
    }

    if (te_bam != NULL) {sam_close(te_bam); te_bam=NULL;}
    if (header != NULL) {sam_hdr_destroy(header); header=NULL;}
    fclose(fileHandle);
}

void constructAiList(AiList *ail, int minCoverageLen)
{
    int minCoverageLen1 = minCoverageLen / 2;
    int minComponentLen = MAX(64, minCoverageLen);
    int numCovered = 0;
    minCoverageLen += minCoverageLen1;
    int numRemaining, bufLen, iter, compIdx, outputIdx, componentStart, idx;

    Contig *contig = &ail->contigList[0];
    Interval *outputIntervals = contig->intervalList;  // to be rebuilt in-place
    int numInterval = contig->numInterval;

    if (numInterval <= minComponentLen) {
        contig->numComp = 1;
        contig->lenComp[0] = numInterval;
        contig->idxComp[0] = 0;
    } else {
        // 1. Decomposition: split intervals into components by coverage
        Interval *inputIntervals = malloc(numInterval * sizeof(Interval));
        Interval *bufferIntervals = malloc(numInterval * sizeof(Interval));
        memcpy(inputIntervals, outputIntervals, numInterval * sizeof(Interval));

        iter = 0;
        outputIdx = 0;
        componentStart = 0;
        numRemaining = numInterval;

        while (iter < MAXC && numRemaining > minComponentLen) {
            bufLen = 0;
            for (idx = 0; idx < numRemaining - minCoverageLen; idx++) {
                int endT = inputIntervals[idx].end;
                compIdx = 1;
                numCovered = 1;
                while (compIdx < minCoverageLen && numCovered < minCoverageLen1) {
                    if (inputIntervals[compIdx + idx].end >= endT) numCovered++;
                    compIdx++;
                }
                if (numCovered < minCoverageLen1)
                    memcpy(&bufferIntervals[bufLen++], &inputIntervals[idx], sizeof(Interval));
                else
                    memcpy(&outputIntervals[outputIdx++], &inputIntervals[idx], sizeof(Interval));
            }
            memcpy(&outputIntervals[outputIdx], &inputIntervals[numRemaining - minCoverageLen],
                minCoverageLen * sizeof(Interval));
            outputIdx += minCoverageLen;
            numRemaining = bufLen;

            contig->idxComp[iter] = componentStart;
            contig->lenComp[iter] = outputIdx - componentStart;
            componentStart = outputIdx;
            iter++;

            if (numRemaining <= minComponentLen || iter == MAXC - 2) {
                if (numRemaining > 0) {
                    memcpy(&outputIntervals[outputIdx], bufferIntervals, numRemaining * sizeof(Interval));
                    contig->idxComp[iter] = outputIdx;
                    contig->lenComp[iter] = numRemaining;
                    iter++;
                }
                contig->numComp = iter;
            } else {
                memcpy(inputIntervals, bufferIntervals, numRemaining * sizeof(Interval));
            }
        }
        free(bufferIntervals);
        free(inputIntervals);
    }

    // 2. Augmentation: build max-end index for each component
    contig->maxEndList = malloc(numInterval * sizeof(int));
    for (compIdx = 0; compIdx < contig->numComp; compIdx++) {
        componentStart = contig->idxComp[compIdx];
        int compEnd = componentStart + contig->lenComp[compIdx];
        int maxEnd = outputIntervals[componentStart].end;
        contig->maxEndList[componentStart] = maxEnd;
        for (idx = componentStart + 1; idx < compEnd; idx++) {
            if (outputIntervals[idx].end > maxEnd) maxEnd = outputIntervals[idx].end;
            contig->maxEndList[idx] = maxEnd;
        }
    }
}


/********************
 *** AiList Query ***
 ********************/
#define isOverlap1(query_start, query_end, interval) ((interval)->start < (query_end) && (interval)->end > (query_start))
#define isOverlap2(query_start, interval) ((interval)->end > (query_start))

int binarySearch(Interval *intervalList, int startIndex, int endIndex, int query_end)
{
    //find targetEnd: index of the rightmost interval satisfying (start < query_end)
    int left = startIndex, right = endIndex-1, middle, targetEnd = -1;

    if(intervalList[right].start < query_end) return right;
    else if(intervalList[left].start >= query_end) return -1;

    while(left < right-1) {
        middle = (left + right) / 2;
        if(intervalList[middle].start >= query_end)
            right = middle - 1;
        else
            left = middle;
    }
    
    if(intervalList[right].start < query_end)
        targetEnd = right;
    else if(intervalList[left].start < query_end)
        targetEnd = left;
        
    return targetEnd; 
}

static inline int getMinDistance(int queryPoint, Interval *interval)
{ return MIN(abs(queryPoint - interval->start), abs(queryPoint - interval->end)); }

static inline void updateMinDistPoint(Interval *interval, int queryPoint, int *numOverlap, int *minDistance)
{
    int prevMinDist = getMinDistance(queryPoint, interval);
    *minDistance = MIN(*minDistance, prevMinDist);
    (*numOverlap)++;
}

void ailistQueryPoint(AiList *ailist, int queryPoint, int flankSize, int *numOverlap, int *minDistance)
{
    int query_start = (queryPoint < flankSize) ? 0 : queryPoint - flankSize;
    int query_end = queryPoint + flankSize;
    Contig *contig = &ailist->contigList[0];
    
    // when there're no intervals, skip
    if (contig->lenComp[0] <= 0)
        return;

    // search intervals in each component one by one
    for(int i = 0; i < contig->numComp; i++){
        int compStart = contig->idxComp[i];
        int compEnd = compStart + contig->lenComp[i];

        if(contig->lenComp[i] <= 15) {
            for(int j = compStart; j < compEnd; j++) {
                Interval *interval = &contig->intervalList[j];
                if(isOverlap1(query_start, query_end, interval))
                    updateMinDistPoint(interval, queryPoint, numOverlap, minDistance);
            }
            continue;
        }

        // j-th interval is the right-most interval with (start < query_end)
        int j = binarySearch(contig->intervalList, compStart, compEnd, query_end);
        while(j >= compStart && contig->maxEndList[j] > query_start) {
            Interval *interval = &contig->intervalList[j];
            if(isOverlap2(query_start, interval))
                updateMinDistPoint(interval, queryPoint, numOverlap, minDistance);
            j--;
        }
    }
}

static inline void updateMinDistInterval(Interval *interval, int start, int end, int *numOverlap, int *minDistance)
{
    int prevMinDistL = getMinDistance(start, interval);
    int prevMinDistR = getMinDistance(end, interval);
    *minDistance = MIN(*minDistance, MIN(prevMinDistL, prevMinDistR));
    (*numOverlap)++;
}

int ailistQueryInterval(AiList *ailist, int start, int end, int flankSize, int *numOverlap, int *minDistance)
{
    int query_start = (start < flankSize) ? 0 : start - flankSize;
    int query_end = end + flankSize;
    int repTid = -1;
    Contig *contig = &ailist->contigList[0];

    // when there're no intervals, skip
    if (contig->lenComp[0] <= 0)
        return repTid;

    // search intervals in each component one by one
    for(int i = 0; i < contig->numComp; i++){
        int compStart = contig->idxComp[i];
        int compEnd = compStart + contig->lenComp[i];

        if(contig->lenComp[i] <= 15) {
            for(int j = compStart; j < compEnd; j++) {
                Interval *interval = &contig->intervalList[j];
                if(isOverlap1(query_start, query_end, interval)) {
                    updateMinDistInterval(interval, start, end, numOverlap, minDistance);
                    repTid = interval->repTid;
                }
            }
            continue;
        }

        int j = binarySearch(contig->intervalList, compStart, compEnd, query_end);
        while(j >= compStart && contig->maxEndList[j] > query_start) {
            Interval *interval = &contig->intervalList[j];
            if(isOverlap2(query_start, interval)) {
                updateMinDistInterval(interval, start, end, numOverlap, minDistance);
                repTid = interval->repTid;
            }
            j--;
        }
    }

    return repTid;
}