#include "AIList.h"
#include "htslib/sam.h"

/***************************
 *** AiList Construction ***
 ***************************/
AiList *initAiList(void)
{
	AiList *ail = malloc(1*sizeof(AiList));
	ail->contigList = malloc(1*sizeof(Contig));
    Contig *firstContig = &ail->contigList[0];

    firstContig->numInterval = 0;
    firstContig->maxIntervals = 64;
	firstContig->intervalList = malloc(firstContig->maxIntervals*sizeof(Interval));
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
{   //New continueous memory?
    int minCoverageLen1=minCoverageLen/2, j1, minL = MAX(64, minCoverageLen);
    minCoverageLen += minCoverageLen1;
    int lenT, len, iter, j, k, k0, t;
	//1. Decomposition
	Contig   *p  = &ail->contigList[0];
	Interval *L1 = p->intervalList;                         //L1: to be rebuilt
	int  numInterval = p->numInterval;
    if(numInterval<=minL){
        p->numComp = 1, p->lenComp[0] = numInterval, p->idxComp[0] = 0;
    }
    else{
        Interval *L0 = malloc(numInterval*sizeof(Interval)); 	//L0: serve as input list
        Interval *L2 = malloc(numInterval*sizeof(Interval));   //L2: extracted list
        memcpy(L0, L1, numInterval*sizeof(Interval));
        iter = 0;	k = 0;	k0 = 0;
        lenT = numInterval;
        while(iter<MAXC && lenT>minL){
            len = 0;
            for(t=0; t<lenT-minCoverageLen; t++){
                int tt = L0[t].end;
                j=1;    j1=1;
                while(j<minCoverageLen && j1<minCoverageLen1){
                    if(L0[j+t].end>=tt) j1++;
                    j++;
                }
                if(j1<minCoverageLen1) memcpy(&L2[len++], &L0[t], sizeof(Interval));
                else memcpy(&L1[k++], &L0[t], sizeof(Interval));
            }
            memcpy(&L1[k], &L0[lenT-minCoverageLen], minCoverageLen*sizeof(Interval));
            k += minCoverageLen, lenT = len;
            p->idxComp[iter] = k0;
            p->lenComp[iter] = k-k0;
            k0 = k, iter++;
            if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                if(lenT>0){
                    memcpy(&L1[k], L2, lenT*sizeof(Interval));
                    p->idxComp[iter] = k;
                    p->lenComp[iter] = lenT;
                    iter++;
                }
                p->numComp = iter;
            }
            else memcpy(L0, L2, lenT*sizeof(Interval));
        }
        free(L2),free(L0);
    }
    //2. Augmentation
    p->maxEndList = malloc(numInterval*sizeof(int));
    for(j=0; j<p->numComp; j++){
        k0 = p->idxComp[j];
        k  = k0 + p->lenComp[j];
        int tt = L1[k0].end;
        p->maxEndList[k0] = tt;
        for(t=k0+1; t<k; t++){
            if(L1[t].end > tt) tt = L1[t].end;
            p->maxEndList[t] = tt;
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