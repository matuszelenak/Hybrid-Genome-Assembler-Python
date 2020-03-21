import argparse
import csv
import itertools
import math
from collections import defaultdict
from typing import Iterator, Dict, Set, List, Tuple

import numpy as np
# from keras.models import load_model

from file_utils import GenomeReadDataReader
from structures import GenomeReadData, GenomeReadCluster, ClusterMergeCase
from utils import iter_with_progress, Cache

KmerIndex = Dict[str, Set[int]]
ClusterConnection = Tuple[int, int, int]


# model = load_model('ecoli_model.hdf5')


cache = Cache()


def get_clusters_from_reads(reads: Iterator[GenomeReadData]) -> List[GenomeReadCluster]:
    cluster_id_counter = itertools.count()
    clusters = []
    for read in reads:
        if not read.characteristic_kmers:
            continue

        clusters.append(
            GenomeReadCluster(
                reference_id=next(cluster_id_counter),
                characteristic_kmers=read.characteristic_kmers.copy(),
                reads=[read]
            )
        )
    print(f'Made {len(clusters)} clusters')
    return clusters


@cache.checkpoint()
def get_kmer_index(clusters: List[GenomeReadCluster]) -> KmerIndex:
    clusters_containing_kmer = defaultdict(set)
    for cluster in iter_with_progress(clusters, total_length=len(clusters), start_message='Building k-mer index for clusters'):
        for kmer in cluster.characteristic_kmers:
            clusters_containing_kmer[kmer].add(cluster.reference_id)

    return clusters_containing_kmer


def get_cluster_connections(clusters: List[GenomeReadCluster], kmer_index: KmerIndex) -> List[ClusterConnection]:
    connections: List[ClusterConnection] = []
    for pivot_cluster in iter_with_progress(clusters, start_message='Calculating connections...'):
        shared_kmer_counts = defaultdict(int)
        for kmer in pivot_cluster.characteristic_kmers:
            for merge_candidate_id in kmer_index[kmer]:
                shared_kmer_counts[merge_candidate_id] += 1

        shared_kmer_counts.pop(pivot_cluster.reference_id)

        connections.extend([(shared_kmers, cluster_id, pivot_cluster.reference_id) for cluster_id, shared_kmers in shared_kmer_counts.items()])

    return sorted(filter(lambda x: x[0] > 1, connections), reverse=True, key=lambda s: s[0])


def clustering_round(clusters: List[GenomeReadCluster], clusters_containing_kmer: KmerIndex) -> Tuple[List[GenomeReadCluster], Set[ClusterMergeCase]]:
    merge_test_cases = []
    cluster_id_to_cluster = {cluster.reference_id: cluster for cluster in clusters}

    cluster_connections = get_cluster_connections(clusters, clusters_containing_kmer)
    # TODO create metric that filters out good connections

    # model_inputs = np.array([
    #     (
    #         cluster_id_to_cluster[cluster_x].size,
    #         cluster_id_to_cluster[cluster_y].size,
    #         len(cluster_id_to_cluster[cluster_x].characteristic_kmers),
    #         len(cluster_id_to_cluster[cluster_x].characteristic_kmers),
    #         score
    #     )
    #     for score, cluster_x, cluster_y in cluster_connections
    # ])
    # should_merge = model.predict(model_inputs).flatten().tolist()

    for score, cluster_x_id, cluster_y_id in iter_with_progress(cluster_connections, start_message='Gathering merge cases...'): #[:math.ceil(len(cluster_connections) * 0.1)]:
        cluster_x: GenomeReadCluster = cluster_id_to_cluster.get(cluster_x_id)
        cluster_y: GenomeReadCluster = cluster_id_to_cluster.get(cluster_y_id)

        merge_case_kwargs = dict(
            cluster_x_size=cluster_x.size,
            cluster_y_size=cluster_y.size,
            cluster_x_kmers=len(cluster_x.characteristic_kmers),
            cluster_y_kmers=len(cluster_y.characteristic_kmers),
            shared_kmers=score,
        )
        if cluster_x.categories != cluster_y.categories:
            merge_test_cases.append(ClusterMergeCase(merge=0, **merge_case_kwargs))
        else:
            merge_test_cases.append(ClusterMergeCase(merge=1, **merge_case_kwargs))

    for score, cluster_x_id, cluster_y_id in iter_with_progress(cluster_connections[:math.ceil(len(cluster_connections) * 0.1)], start_message='Merging clusters...'): #[:math.ceil(len(cluster_connections) * 0.1)]:
        cluster_x: GenomeReadCluster = cluster_id_to_cluster.get(cluster_x_id)
        cluster_y: GenomeReadCluster = cluster_id_to_cluster.get(cluster_y_id)

        if not (cluster_x and cluster_y):
            continue
        #
        # if cluster_x.categories != cluster_y.categories:
        #     continue

        # if merge < 0.8:
        #     continue

        if cluster_x.size > cluster_y.size:
            bigger, smaller = cluster_x, cluster_y
        else:
            bigger, smaller = cluster_y, cluster_x

        for kmer in smaller.characteristic_kmers:
            clusters_containing_kmer[kmer].remove(smaller.reference_id)
            clusters_containing_kmer[kmer].add(bigger.reference_id)

        bigger |= smaller
        cluster_id_to_cluster.pop(smaller.reference_id)

    return list(cluster_id_to_cluster.values()), set(merge_test_cases)


def run_clustering(clusters: List[GenomeReadCluster], kmer_index: KmerIndex) -> List[ClusterMergeCase]:
    num_of_components = len(clusters)
    all_merge_cases = []
    while True:
        clusters, merge_cases = clustering_round(clusters, kmer_index)
        print([cluster.consistency for cluster in sorted(clusters, key=lambda cl: cl.size, reverse=True)[:100]])
        all_merge_cases.extend(merge_cases)

        if len(clusters) == num_of_components:
            break

        num_of_components = len(clusters)

    print(f'Ended up with {len(clusters)} clusters')

    return all_merge_cases


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='File path of the generated read file')
args = parser.parse_args()

cache.global_kwargs['read_data'] = args.filename.replace('/', '')

with GenomeReadDataReader(args.filename) as r:
    c = get_clusters_from_reads(r.get_read_data())
    index = get_kmer_index(c)
    cases = run_clustering(c, index)

with open(f'{args.filename}__merge_cases.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(map(tuple, cases))
