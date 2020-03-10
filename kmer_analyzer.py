import argparse
import math
from collections import defaultdict
from functools import lru_cache
from typing import Iterable, List, Tuple, Dict, Generator

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
    def kmer_occurrences(self, k: int) -> Dict[str, int]:
        kmer_occurrences = defaultdict(int)
        for genome_read in iter_with_progress(self.reads, total_length=self.meta.num_of_reads, start_message=f'Computing {k}-mer occurrences'):
            for kmer in genome_read.iterate_kmers(k):
                kmer_occurrences[kmer.sequence] += 1

        return kmer_occurrences

    def get_annotated_read_data(self, k: int, lower: int, upper: int) -> Generator[GenomeReadData, None, None]:
        characteristic_kmers = set(sequence for sequence, count in self.kmer_occurrences(k).items() if lower <= count <= upper)

        for genome_read in iter_with_progress(self.reads, total_length=self.meta.num_of_reads, start_message=f'Indexing characteristic {k}-mers in reads'):
            data = GenomeReadData(
                label=genome_read.label
            )
            data.index_characteristic_kmers(genome_read, characteristic_kmers, len(characteristic_kmers) // 10)
            yield data


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
    initial_k = math.ceil(math.log(r.meta.genome_size, len(BASES)))
    read_set = GenomeReadSet(r)

    plot_occurrence_histograms([
        (k_option, list(filter(lambda value: value <= r.meta.coverage * 3, read_set.kmer_occurrences(k_option).values())))
        for k_option in range(initial_k, initial_k + 2)
    ])

    k_length, cov_low, cow_high = [int(x) for x in input().split()]

    read_data = list(read_set.get_annotated_read_data(k_length, cov_low, cow_high))

    with GenomeReadDataWriter(f'{args.filename}__read_data.json') as w:
        w.write(read_data)
