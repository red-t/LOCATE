#ifndef LTR_UTILS_H
#define LTR_UTILS_H

#include <stdlib.h>
#include "io_utils.h"

/// @brief Define and return LTR sizes for all LTR-class TE
void define_ltr(const char *te_fn, const char *teClassFn);

/// @brief Load TE class information
int *get_class_arr(int numTe, const char *teClassFn);

/// @brief Output left-/right- half of LTR-class TE separately
void outputSeq(int tid, faidx_t *teFa);

/// @brief Map right-half to left-half
void mapEachOther(int tid);

/// @brief Compute LTR region length
int getLtrLen(int tid);

#endif // LTR_UTILS_H