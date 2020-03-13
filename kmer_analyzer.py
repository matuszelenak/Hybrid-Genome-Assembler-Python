import argparse
import json
import math
from collections import defaultdict, OrderedDict
from functools import lru_cache
from typing import Iterable, List, Tuple, Dict, Generator

import numpy as np
from matplotlib import pyplot as plt

from file_utils import GenomeReader, GenomeReadDataWriter
from structures import GenomeRead, BASES, GenomeReadData
from utils import iter_with_progress, split_into_chunks


class GenomeReadSet:
    reader: GenomeReader = None

    def __init__(self, reader: GenomeReader):
        self.reader = reader
        self.meta = self.reader.meta

    @property
    def reads(self) -> Iterable[GenomeRead]:
        yield from self.reader.get_reads()

    @lru_cache(maxsize=10)
    def kmer_occurrences(self, k: int) -> Dict[str, Dict[str, int]]:
        kmer_occurrences = defaultdict(lambda: defaultdict(int))
        for genome_read in iter_with_progress(self.reads, total_length=self.meta.num_of_reads, start_message=f'Computing {k}-mer occurrences'):
            for kmer in genome_read.iterate_kmers(k):
                kmer_occurrences[kmer][genome_read.category] += 1

        return kmer_occurrences

    def get_annotated_read_data(self, k: int, lower: int, upper: int) -> Generator[GenomeReadData, None, None]:
        characteristic_kmers = set(sequence for sequence, counts in self.kmer_occurrences(k).items() if lower <= sum(counts.values()) <= upper)

        for genome_read in iter_with_progress(self.reads, total_length=self.meta.num_of_reads, start_message=f'Indexing characteristic {k}-mers in reads'):
            read_characteristic_kmers = set()
            for kmer in genome_read.iterate_kmers(k):
                if kmer in characteristic_kmers:
                    read_characteristic_kmers.add(kmer)

            data = GenomeReadData(label=genome_read.label, category=genome_read.category, characteristic_kmers=read_characteristic_kmers)
            yield data


def get_kmer_specificity_levels(kmer_occurrences: Dict[str, Dict[str, int]], max_coverage: int):
    breakpoints = [50, 70, 80, 90, 95, 99, 100, 101]
    colors = ['k', 'tab:purple', 'r', 'xkcd:orange', 'y', 'b', 'g']

    def get_bracket_key(value):
        for low_point, high_point in zip(breakpoints, breakpoints[1:]):
            if low_point <= value < high_point:
                return f'{low_point} <= x < {high_point}'

    coverage_to_specificity_bins: Dict[int, Dict[str, List[str]]] = defaultdict(lambda: OrderedDict({
        f'{low} <= x < {up}': []
        for low, up in zip(breakpoints, breakpoints[1:])
    }))
    for kmer, counts in kmer_occurrences.items():
        total_count = sum(counts.values())
        prevalent_count = max(counts.values())
        specificity = round(prevalent_count / total_count * 100)
        coverage_to_specificity_bins[total_count][get_bracket_key(specificity)].append(kmer)

    indices = np.arange(max_coverage)
    bottoms = np.zeros((max_coverage,))
    bracket_set = set(get_bracket_key(x) for x in breakpoints[:-1])
    plots = []
    for bracket_key, color in zip(bracket_set, colors):
        specificity_level_data = np.array(
            [
                len(coverage_to_specificity_bins[coverage].get(bracket_key, []))
                for coverage in range(max_coverage)
            ]
        )

        plots.append(plt.bar(indices, specificity_level_data, width=0.3, color=color))
        bottoms += specificity_level_data

    plt.legend(plots, list(bracket_set))
    plt.show()


def plot_occurrence_histograms(datasets: List[Tuple[int, List[int]]]):
    rounded_up_num_of_graphs = math.ceil(len(datasets) / 3) * 3

    fig, axs = plt.subplots(rounded_up_num_of_graphs, sharex='col', sharey='row')
    axs = axs.reshape((rounded_up_num_of_graphs // 3, 3))

    for row_num, dataset_row in enumerate(split_into_chunks(datasets, 3)):
        for col_num, (k, histogram_data) in enumerate(dataset_row):
            axs[row_num, col_num].set_title(f'Occurrences of characteristic {k}-mers in reads')
            axs[row_num, col_num].hist(histogram_data, bins=len(set(histogram_data)) or 1)

    plt.show()


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='File path of the generated read file')
args = parser.parse_args()

with GenomeReader(args.filename) as r:
    initial_k = math.ceil(math.log(r.meta.genome_size, len(BASES))) + 1
    read_set = GenomeReadSet(r)

    for k_option in range(initial_k, initial_k + 2):
        occ = read_set.kmer_occurrences(k_option)
        get_kmer_specificity_levels(occ, max_coverage=r.meta.coverage * 2)

    k_length, cov_low, cov_high = [int(x) for x in input().split()]
    read_data = list(read_set.get_annotated_read_data(k_length, cov_low, cov_high))

    with GenomeReadDataWriter(f'{args.filename}__read_data.json') as w:
        w.write(read_data)

    with open(f'{args.filename}__kmers.json', 'w') as f:
        json.dump({key: counts for key, counts in read_set.kmer_occurrences(k_length).items() if cov_low <= sum(counts.values()) <= cov_high}, f, indent=4)
