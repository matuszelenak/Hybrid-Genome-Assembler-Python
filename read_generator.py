import argparse
import os
import random
import itertools
import subprocess
from typing import List, Tuple, Generator, Optional, Iterable

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna as DNA

from file_utils import get_file_type
from structures import GenomeReadsMetaData


class ReadGeneratorBackend:
    @staticmethod
    def get_num_of_reads(seq_length: int, coverage: int, read_length: int) -> int:
        return int(coverage * seq_length / read_length)

    @staticmethod
    def reads_from_sequence(record: SeqRecord, coverage: int, read_length: int, label_prefix: Optional[str] = None) -> Generator[SeqRecord, None, None]:
        raise NotImplementedError


class ArtIlluminaBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, coverage: int, read_length: int, label_prefix: Optional[str] = None) -> Generator[SeqRecord, None, None]:
        temporary_file = f'{record.id}_tmp.fasta'
        output_file = f'{record.id}_tmp'
        SeqIO.write([record], temporary_file, 'fasta')

        command = ['art_illumina', '-ss', 'HS25', '-l', str(read_length), '-f', str(coverage), '-na', '-i', temporary_file, '-o', output_file]
        subprocess.Popen(command).wait()

        for read_record in SeqIO.parse(f'{output_file}.fq', 'fastq'):
            read_record.description = record.id
            yield read_record

        os.remove(temporary_file)
        os.remove(f'{output_file}.fq')


class CustomReadGenerator(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, coverage: int, read_length: int, label_prefix: Optional[str] = None) -> Generator[SeqRecord, None, None]:
        num_of_reads = ReadGeneratorBackend.get_num_of_reads(len(record), coverage, read_length)
        read_beginnings = np.random.randint(0, len(record) - read_length, size=num_of_reads)

        label_prefix = label_prefix or record.id
        for read_id, beginning in enumerate(read_beginnings):
            yield SeqRecord(
                record.seq[beginning: beginning + read_length],
                id=f'{label_prefix}-{read_id}',
                description=record.id
            )


get_backend = {'custom': CustomReadGenerator, 'art': ArtIlluminaBackend}.get


def generate_mutated_sequence_pair(genome_size: int, difference_rate: float) -> Tuple[SeqRecord, SeqRecord]:
    def indices_to_strand(indices: List[int]) -> str:
        return ''.join(DNA.letters[index] for index in indices)

    original_sequence_indices = np.random.choice(len(DNA.letters), genome_size)
    # Choose which positions to mutate and shift those positions by 1 to 3 slots
    mutation_index_shift = np.random.binomial(1, difference_rate, genome_size) * np.random.randint(1, 4, size=genome_size)

    mutated_sequence_indices = (original_sequence_indices + mutation_index_shift) % len(DNA.letters)

    original_sequence = indices_to_strand(original_sequence_indices)
    mutated_sequence = indices_to_strand(mutated_sequence_indices)
    return SeqRecord(Seq(original_sequence, DNA), id='Original'), SeqRecord(Seq(mutated_sequence, DNA), id='Mutated')


def mix_reads(*read_lists: List[SeqRecord], shuffle=False) -> Iterable[SeqRecord]:
    all_reads = itertools.chain(*read_lists)
    if shuffle:
        all_reads = list(all_reads)
        random.shuffle(all_reads)
    return all_reads


parser = argparse.ArgumentParser()
# Reference genomes parameters
parser.add_argument('-i', dest='input_files', type=str, nargs=2, help='Paths to files with reference genomes', required=False)
parser.add_argument('-n', dest='genome_size', nargs='?', type=int, help='Size of the two reference genomes')
parser.add_argument('-d', dest='difference_rate', nargs='?', type=float, help='Difference rate between the two sequences')

# Read parameters
parser.add_argument('-b', dest='backend', default='Custom', choices=['custom', 'art'], help='Read generation backend: Either Art Illumina or Custom')
parser.add_argument('-c', dest='coverage', nargs='?', default=30, type=int, help='Coverage of the reads')
parser.add_argument('-r', dest='read_length', nargs='?', default=150, type=int, help='Read length')

parser.add_argument('-o', dest='output_filename', type=str, help='Name of the output file with reads', default='out.fasta')
args = parser.parse_args()

if not args.input_files and not all([args.genome_size, args.difference_rate]):
    raise ValueError('Either specify 2 input files or parameters of the genome for mutation')


meta_kwargs = {}
if args.input_files:
    sequences = [SeqIO.read(filename, format=get_file_type(filename)) for filename in args.input_files]
    meta_kwargs.update({
        'genome_size': max(len(record) for record in sequences)
    })
else:
    original, mutated = generate_mutated_sequence_pair(args.genome_size, args.difference_rate)
    sequences = [original, mutated]
    meta_kwargs.update({
        'genome_size': args.genome_size,
        'difference': args.difference_rate,
    })

backend: ReadGeneratorBackend = get_backend(args.backend)
reads = itertools.chain.from_iterable(backend.reads_from_sequence(seq, args.coverage, args.read_length) for seq in sequences)
read_count = SeqIO.write(reads, args.output_filename, 'fasta')

meta_kwargs.update({
    'read_length': args.read_length,
    'coverage': args.coverage,
    'num_of_reads': read_count,
    'categories': [record.id for record in sequences],
    'alphabet': DNA.letters
})

meta = GenomeReadsMetaData(**meta_kwargs)
with open(f'{args.output_filename}_meta.json', 'w') as f:
    f.write(str(meta))
