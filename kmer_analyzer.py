import argparse
import math
from collections import defaultdict
from functools import lru_cache, partial
from typing import Iterable, Tuple, Dict, Generator

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt

from file_utils import GenomeReadDataWriter
from structures import GenomeReadData, GenomeReadsMetaData
from utils import iter_with_progress, iterate_kmers_in_record, minimum_and_average, Cache, dd

KmerOccurrences = Dict[str, Dict[str, int]]
KmerSpecificity = Dict[int, Dict[Tuple[int, int], int]]

MIN_QUALITY = 10
AVG_MIN_QUALITY = 25
min_and_avg = partial(partial(minimum_and_average, MIN_QUALITY), AVG_MIN_QUALITY)
cache = Cache()


class KmerHistogramPlotter:
    brackets = (
        (50, 70, 'k'),
        (70, 85, 'r'),
        (85, 90, 'xkcd:orange'),
        (90, 95, 'y'),
        (95, 100, 'b'),
        (100, 101, 'g')
    )

    @classmethod
    @lru_cache(maxsize=150)
    def get_value_bracket(cls, value) -> Tuple[int, int]:
        for low, up, _ in cls.brackets:
            if low <= value < up:
                return low, up

    @classmethod
    def get_kmer_specificity(cls, occurrences: KmerOccurrences) -> KmerSpecificity:
        coverage_to_specificity_bins: KmerSpecificity = defaultdict(lambda: defaultdict(int))
        for kmer, counts in occurrences.items():
            total_count = sum(counts.values())
            prevalent_count = max(counts.values())
            specificity = prevalent_count / total_count * 100
            coverage_to_specificity_bins[total_count][cls.get_value_bracket(specificity)] += 1

        return coverage_to_specificity_bins

    @classmethod
    def render_subplot(cls, subplot, coverage_specificity: KmerSpecificity, max_coverage: int):
        bars = []
        indices = np.arange(max_coverage)
        bottoms = np.zeros((max_coverage,))
        for low, high, color in cls.brackets:
            specificity_level_data = np.array(
                [
                    coverage_specificity[coverage].get((low, high), 1) if coverage > 2 else 1
                    for coverage in range(max_coverage)
                ]
            )
            bars.append(subplot.bar(indices, specificity_level_data, bottom=bottoms, width=0.3, color=color))
            bottoms += specificity_level_data

        # subplot.set_yscale('log')
        subplot.legend(bars, [f'{low} ≤ x < {high}' for low, high, _ in cls.brackets])

    @classmethod
    def plot_histograms(cls, k_value_occurrences: Dict[int, KmerOccurrences], max_coverage: int):
        fig, axs = plt.subplots(len(list(k_value_occurrences.keys())), sharex='col', sharey='row')
        k_values = sorted(list(k_value_occurrences.keys()))
        for row, k_value in enumerate(k_values):
            coverage_specificity = cls.get_kmer_specificity(k_value_occurrences[k_value])
            axs[row].set_title(f'Occurrences of characteristic {k_value}-mers in reads ({len(k_value_occurrences[k_value])} total)')
            cls.render_subplot(axs[row], coverage_specificity, max_coverage)

        plt.show()


@cache.checkpoint(min_quality=MIN_QUALITY, avg_quality=AVG_MIN_QUALITY)
def kmer_occurrences(reads: Iterable[SeqRecord], k: int, read_count: int = None) -> KmerOccurrences:
    occurrences = defaultdict(dd)
    cnt = 0
    total = 0
    for read_record in iter_with_progress(reads, total_length=read_count, start_message=f'Computing occurrences of {k}-mers'):
        for kmer in iterate_kmers_in_record(read_record, k, min_and_avg):
            occurrences[kmer][read_record.description.split(' ')[-1]] += 1
            cnt += 1

        total += (len(read_record) - (k - 1))

    print(f'Gathered {cnt} out of {total} kmers')

    return occurrences


@cache.checkpoint(min_quality=MIN_QUALITY, avg_quality=AVG_MIN_QUALITY)
def get_annotated_read_data(reads: Iterable[SeqRecord], occurrences: KmerOccurrences, lower: int, upper: int, read_count: int) -> Generator[GenomeReadData, None, None]:
    characteristic_kmers = set(sequence for sequence, counts in occurrences.items() if lower <= sum(counts.values()) <= upper)
    k = len(next(iter(characteristic_kmers)))

    for read_record in iter_with_progress(reads, total_length=read_count, start_message='Annotating read data...'):
        read_characteristic_kmers = set()
        for kmer in iterate_kmers_in_record(read_record, k, min_and_avg):
            if kmer in characteristic_kmers:
                read_characteristic_kmers.add(kmer)

        yield GenomeReadData(label=str(read_record.id), category=read_record.description.split(' ')[-1], characteristic_kmers=read_characteristic_kmers)


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='File path of the generated read file')
args = parser.parse_args()

cache.global_kwargs['filename'] = args.filename.replace('/', '')

meta = GenomeReadsMetaData.from_file(f'{args.filename}_meta.json')

k_guess = math.ceil(math.log(meta.genome_size, len(meta.alphabet))) + 3
k_value_occurrences = {
    k_option: kmer_occurrences(SeqIO.parse(args.filename, 'fastq'), k_option, meta.num_of_reads)
    for k_option in range(k_guess, k_guess + 2)
}

KmerHistogramPlotter.plot_histograms(k_value_occurrences, meta.coverage * 4)

print('Enter the selected k, minimal and maximal coverage')
response = input()
if response == 'q':
    exit()

k_length, cov_low, cov_high = [int(x) for x in response.split()]
read_data = list(get_annotated_read_data(SeqIO.parse(args.filename, 'fastq'), k_value_occurrences[k_length], cov_low, cov_high, meta.num_of_reads))

with GenomeReadDataWriter(f'{args.filename}__read_data.json') as w:
    w.write(read_data)
