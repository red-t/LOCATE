#ifndef LOCATE_CONSTANTS_H
#define LOCATE_CONSTANTS_H

/****************************************
 *** Named Constants for LOCATE Pipeline ***
 ****************************************/

/* Flank sizes */
#define FLANK_SIZE_50            50
#define FLANK_SIZE_100           100
#define FLANK_SIZE_200           200
#define FLANK_SIZE_500           500
#define FLANK_SIZE_3000          3000

/* PolyA / PolyT detection */
#define TRUNC_PROPORTION         0.08
#define TRUNC_MIN_THRESHOLD      100
#define MIN_POLYA_LEN            5
#define POLYA_FDR_THRESHOLD      0.01
#define POLYA_BASE_PROB          0.25
#define POLYA_NUM_OTHER_MAX      10
#define POLYA_GAP_THRESHOLD      20

/* Alignment location types */
#define LOCATION_NORMAL          1
#define LOCATION_BOUNDARY        2
#define LOCATION_INSIDE          4
#define LOCATION_ONE_BOUNDARY    8
#define LOCATION_ONE_INSIDE      16

/* Cluster types */
#define CLT_TYPE_GERMLINE        0
#define CLT_TYPE_SOMATIC_1       1
#define CLT_TYPE_SOMATIC_2       2

/* Assembly */
#define SAMTOOLS_FILTER_FLAG     3332
#define ASSM_PRIMARY_RETRIES     5
#define ASSM_DEFAULT_NODE_LEN    256
#define ASSM_MAX_POLYMER_LEN     20

/* Clip / flank adjustment thresholds */
#define CLIP_LARGE_THRESHOLD     400
#define FLANK_BREAKPOINT_DELTA   50
#define TSD_BREAKPOINT_DELTA     50
#define FLANK_SIZE_20            20

/* Dynamic array resizing */
#define RESIZE_FACTOR            1.5
#define THRESHOLD_FRACTION       0.9

/* Default capacities for dynamic arrays */
#define CAPACITY_INIT_SEG        10000
#define CAPACITY_INIT_CLT        10000
#define CAPACITY_INIT_TE         10000
#define CAPACITY_HIGH_QUAL       2000
#define CAPACITY_ANNO            3000
#define CAPACITY_CLASS           200
#define CAPACITY_SIZE            200

/* AiList construction */
#define AI_MIN_COMPONENT_LEN     64
#define AI_MIN_COVERAGE_LEN      20

/* TE alignment fractions */
#define TE_ALIGNED_FRAC_CUTOFF   0.8
#define GAP_CUTOFF_FRACTION      0.2

/* Iteration ranges */
#define ITER_RANGE_FREQ          50000
#define SEG_FLANK_SIZE           3000

/* Minimum segment / overhang */
#define MIN_OVERHANG             100
#define MIN_SEG_LEN              100

/* Insertion gap filler */
#define DIFF_FLANK_GAP_FILLER    500
#define SINGLE_FLANK_ADJUST      50
#define MIN_INS_FLANK_LEN        250

/* File path buffer sizes */
#define FN_BUF_SIZE              100
#define STR_BUF_SIZE_Q           500
#define STR_BUF_SIZE_L           100000

#endif /* LOCATE_CONSTANTS_H */
