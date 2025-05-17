from .FileIO cimport *
from cpython cimport PyBytes_FromStringAndSize

cpdef assemble_cluster(Cluster[::1] clt_view, dict cluster_data_by_tid, tuple block, object cmd_args, int extra_thread)