from .FileIO cimport *
from .Cluster cimport bamIsInvalid, getMapLenAndDiv, getHighQualClts

cdef extern from "src/ltr_utils.h":
    void defineLTR(const char *teFn, const char *teClassFn)

cpdef object runInParallel(object cmdArgs)