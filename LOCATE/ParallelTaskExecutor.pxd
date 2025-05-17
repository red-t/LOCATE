from .FileIO cimport *
from .Cluster cimport bam_is_invalid, get_mapping_length_and_divergence, get_high_quality_clusters

cdef extern from "src/ltr_utils.h":
    void define_ltr(const char *te_fn, const char *class_fn)

cpdef object run_in_parallel(object cmd_args)