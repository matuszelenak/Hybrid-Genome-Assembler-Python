import argparse
import random
from typing import List, Tuple, Generator

import numpy as np

from file_utils import GenomeWriter
from structures import GenomeMetaData, BASES, GenomeRead


class GenomeReadGenerator:
    def __init__(self, genome_size: int, coverage: int, read_length: int, difference_rate: float, circular: bool = True):
        self.meta = GenomeMetaData(
            genome_size=genome_size,
            coverage=coverage,
            read_length=read_length,
            num_of_reads=self.get_num_of_reads(genome_size, coverage, read_length),
            difference=difference_rate,
            circular=circular
        )

    @staticmethod
    def generate_sequences(length: int = 1000, difference_rate: float = 0) -> Tuple[List[str], List[str]]:
        sequence_a, sequence_b = [], []
        for base, mutation_roll in zip(random.choices(BASES, k=length), np.random.sample(length)):
            if mutation_roll < difference_rate:
                mutated_base = BASES[(BASES.index(base) + random.randint(1, len(BASES) - 1)) % len(BASES)]
            else:
                mutated_base = base

            sequence_a.append(base)
            sequence_b.append(mutated_base)

        return sequence_a, sequence_b

    @staticmethod
    def get_num_of_reads(seq_length: int, coverage: int, read_length: int):
        return int(coverage * seq_length / read_length)

    def reads_from_sequence(self, sequence: List[str], coverage: int, read_length: int, label_prefix: str, circular: bool = True) -> Generator[GenomeRead, None, None]:
        num_of_reads = self.get_num_of_reads(len(sequence), coverage, read_length)
        if circular:
            read_beginnings = np.random.randint(0, len(sequence), size=num_of_reads)
        else:
            read_beginnings = np.random.randint(0, len(sequence) - read_length, size=num_of_reads)

        for read_id, beginning in enumerate(read_beginnings):
            read_seq = sequence[beginning:beginning + read_length]
            if beginning + read_length > len(sequence):
                read_seq += sequence[:(beginning + read_length) % len(sequence)]
            yield GenomeRead(label=f'{label_prefix}-{read_id}', sequence=''.join(read_seq))

    def get_reads(self):
        sequence_a, sequence_b = self.generate_sequences(length=self.meta.genome_size, difference_rate=self.meta.difference)
        reads_a = self.reads_from_sequence(sequence_a, self.meta.coverage // 2, self.meta.read_length, 'A', self.meta.circular)
        reads_b = self.reads_from_sequence(sequence_b, self.meta.coverage // 2, self.meta.read_length, 'B', self.meta.circular)
        a_exhausted, b_exhausted = False, False
        while True:
            if random.randint(0, 1) == 0:
                try:
                    yield next(reads_a)
                except StopIteration:
                    a_exhausted = True
            else:
                try:
                    yield next(reads_b)
                except StopIteration:
                    b_exhausted = True

            if a_exhausted and b_exhausted:
                return


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str, help='File path of the generated read file')
parser.add_argument('-n', dest='genome_size', nargs='?', default=10000, type=int, help='Size of the two reference genomes')
parser.add_argument('-c', dest='coverage', nargs='?', default=30, type=int, help='Coverage of the reads')
parser.add_argument('-r', dest='read_length', nargs='?', default=150, type=int, help='Read length')
parser.add_argument('-d', dest='difference', nargs='?', default=0.05, type=float, help='Difference rate between the two sequences')
parser.add_argument('-C', action='store_true', help='Circular genome')
args = parser.parse_args()

generator = GenomeReadGenerator(args.genome_size, args.coverage, args.read_length, args.difference)
with GenomeWriter(args.filename, generator.meta) as writer:
    writer.write(generator.get_reads())
