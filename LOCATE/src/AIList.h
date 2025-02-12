//=====================================================================================
//Common structs, parameters, functions
//by Jianglin Feng  09/5/2018
//modified by Zhongren Hu	14/3/2023
//-------------------------------------------------------------------------------------
#ifndef __AILIST_H__
#define __AILIST_H__
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "khash.h"


/*********************
 * Convenient macros *
 *********************/
#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	} while (0)

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAXC 10	//max number of components


/******************
 *** Structures ***
 ******************/
typedef struct {
    int start;
    int end;
	int repTid;
} Interval;

typedef struct{
	int64_t numInterval, maxIntervals;
	Interval *intervalList;
	int numComp, lenComp[MAXC], idxComp[MAXC]; 
	int *maxEndList;
} Contig;

/// Augmented Interval List.
typedef struct {
	Contig *contigList; // default only 1 contig
} AiList;


/***************************
 *** AiList Construction ***
 ***************************/
AiList *initAiList(void);
void destroyAiList(AiList *ail);
void readBED(AiList *ail, const char *bedFn, const char *targetChrom);
void constructAiList(AiList *ail, int minCoverageLen);


/********************
 *** AiList Query ***
 ********************/
/// Compute the minimum distance between query position to the nearest interval
void ailistQueryPoint(AiList *ailist, int queryPoint, int flankSize, int *numOverlap, int *minDistance);

/// Compute the minimum distance between query interval to the nearest interval
int ailistQueryInterval(AiList *ailist, int start, int end, int flankSize, int *numOverlap, int *minDistance);

#endif