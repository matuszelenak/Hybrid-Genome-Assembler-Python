from typing import NamedTuple


class ClusterMergeCase(NamedTuple):
    cluster_x_size: int
    cluster_y_size: int
    cluster_x_kmers: int
    cluster_y_kmers: int
    shared_kmers: int
