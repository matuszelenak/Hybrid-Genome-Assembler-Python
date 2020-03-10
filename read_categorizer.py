import argparse
import itertools
import math
import random
from typing import List, Tuple, Dict

from file_utils import GenomeReadDataReader
from structures import GenomeReadCluster, GenomeReadData
from utils import iter_with_progress, split_into_chunks


def read_cluster_shared_kmers_count(cluster_a: GenomeReadCluster, cluster_b: GenomeReadCluster) -> int:
    return (cluster_a.characteristic_kmers & cluster_b.characteristic_kmers).estimated_item_count


class Unionizer:
    parents: Dict[GenomeReadCluster, GenomeReadCluster]

    def _get_parent(self, cluster: GenomeReadCluster) -> GenomeReadCluster:
        if self.parents[cluster] == cluster:
            return cluster

        self.parents[cluster] = self._get_parent(self.parents[cluster])
        return self.parents[cluster]

    def _unite(self, cluster_x: GenomeReadCluster, cluster_y: GenomeReadCluster) -> bool:
        parent_x, parent_y = self._get_parent(cluster_x), self._get_parent(cluster_y)
        if parent_x == parent_y:
            return False

        if parent_x.size > parent_y.size:
            bigger, smaller = parent_x, parent_y
        else:
            bigger, smaller = parent_y, parent_x

        bigger |= smaller
        self.parents[smaller] = bigger
        smaller.clear()

        return True

    @staticmethod
    def _calculate_cross_cluster_scores(read_clusters: List[GenomeReadCluster]) -> List[Tuple[int, GenomeReadCluster, GenomeReadCluster]]:
        for cluster_x, cluster_y in itertools.combinations(read_clusters, 2):# iter_with_progress(itertools.combinations(read_clusters, 2), total_length=len(read_clusters) * len(read_clusters), start_message='Calculating cross-cluster scores'):
            yield read_cluster_shared_kmers_count(cluster_x, cluster_y), cluster_x, cluster_y

    def run(self, read_clusters: List[GenomeReadCluster], final_cluster_count: int, minimum_shared_kmers: int = 4) -> List[GenomeReadCluster]:
        group_scores = list(self._calculate_cross_cluster_scores(read_clusters))
        group_scores.sort(key=lambda s: s[0], reverse=True)
        group_scores = filter(lambda x: x[0] >= minimum_shared_kmers, group_scores)

        self.parents = {cluster: cluster for cluster in read_clusters}

        num_of_components = len(read_clusters)
        for score, cluster_x, cluster_y in iter_with_progress(group_scores):
            x_labels, y_labels = cluster_x.labels, cluster_y.labels
            if self._unite(cluster_x, cluster_y):
                # print(f'Joined clusters {x_labels} and {y_labels} on {score} shared')
                num_of_components -= 1

            if num_of_components <= final_cluster_count:
                break

        return [cluster for cluster in read_clusters if cluster.size > 0]


class ReadClusterer:
    def __init__(self, read_data: List[GenomeReadData]):
        self.read_data = [data for data in read_data if data.characteristic_kmers_count > 0]

    def run_clustering(self):
        # Sort reads by the number of characteristic kmers they contain
        read_data_sorted_by_occurrences = sorted(
            self.read_data,
            key=lambda read: read.characteristic_kmers_count,
            reverse=True
        )
        # Since many reads will contain 0 or only a few, take only the top 30% of the reads
        read_data_sorted_by_occurrences = read_data_sorted_by_occurrences[:int(len(read_data_sorted_by_occurrences) * 0.3)]

        reference_id_generator = itertools.count(0)
        read_clusters = [
            GenomeReadCluster(reference_id=next(reference_id_generator), reads=[read], characteristic_kmers=read.characteristic_kmers)
            for read in read_data_sorted_by_occurrences
        ]

        # Split the reads into groups of size sqrt(Number of reads) (to utilize divide and conquer method)
        first_level_group_size = math.ceil(math.sqrt(len(read_clusters)))
        # first_level_group_size = 3000
        cluster_groups = split_into_chunks(read_clusters, first_level_group_size)

        for group in list(cluster_groups):
            u = Unionizer()
            clusters = u.run(group, 2)
            for cluster in sorted(clusters, key=lambda c: c.size, reverse=True)[:10]:
                print(cluster.stats)


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='File path of the generated read file')
args = parser.parse_args()

with GenomeReadDataReader(args.filename) as r:
    read_data = list(r.get_read_data())
    cl = ReadClusterer(read_data)
    cl.run_clustering()
