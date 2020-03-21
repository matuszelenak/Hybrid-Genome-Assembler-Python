import argparse
import os
import subprocess
from typing import List, Tuple

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna as DNA

from file_utils import get_file_type
from structures import GenomeReadsMetaData
from utils import sequence_complement


SIMLORD_BIN = '/home/whiskas/miniconda3/envs/simlord/bin/simlord'


class ReadGeneratorBackend:
    @staticmethod
    def get_num_of_reads(seq_length: int, coverage: int, read_length: int) -> int:
        return int(coverage * seq_length / read_length)

    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True) -> int:
        raise NotImplementedError


class ArtIlluminaBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True) -> int:
        temporary_file = f'{record.id}_tmp.fasta'
        SeqIO.write([record], temporary_file, 'fasta')

        command = ['art_illumina', '-ss', 'HS25', '-l', str(read_length), '-f', str(coverage), '-na', '-i', temporary_file, '-o', output_prefix]

        print(command)
        subprocess.Popen(command).wait()

        os.remove(temporary_file)

        return sum(1 for _ in SeqIO.parse(f'{output_prefix}.fq', 'fastq'))


class CustomReadGenerator(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True) -> int:
        num_of_reads = ReadGeneratorBackend.get_num_of_reads(len(record), coverage, read_length)
        read_beginnings = np.random.randint(0, len(record) - read_length, size=num_of_reads)

        strands = [str(record.seq), sequence_complement(str(record.seq))]
        strand_choices = np.random.binomial(1, 0.5, len(read_beginnings))

        SeqIO.write(
            (
                SeqRecord(
                    Seq(strands[strand_choice][beginning: beginning + read_length]),
                    id=f'{record.id}-{read_id}',
                    description='',
                    letter_annotations={'phred_quality': [41 for _ in range(read_length)]}
                )
                for read_id, (beginning, strand_choice) in enumerate(zip(read_beginnings, strand_choices))
            ),
            f'{output_prefix}.fq',
            'fastq'
        )

        return num_of_reads


class PacBioBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True) -> int:
        temporary_file = f'{record.id}_tmp.fasta'
        SeqIO.write([record], temporary_file, 'fasta')

        cmd = [SIMLORD_BIN, '-rr', temporary_file, '--no-sam', '-fl', str(read_length), '--without-ns', '-c', str(coverage), '-pi', '0.001', '-pd', '0.001', '-ps', '0.01', output_prefix]
        subprocess.Popen(cmd).wait()

        os.rename(f'{output_prefix}.fastq', f'{output_prefix}.fq')

        os.remove(temporary_file)

        return sum(1 for _ in SeqIO.parse(f'{output_prefix}.fq', 'fastq'))


backends = (
    ('custom', CustomReadGenerator),
    ('art', ArtIlluminaBackend),
    ('pacbio', PacBioBackend)
)
get_backend = dict(backends).get


def generate_mutated_sequence_pair(genome_size: int, difference_rate: float) -> Tuple[SeqRecord, SeqRecord]:
    def indices_to_strand(indices: List[int]) -> str:
        return ''.join(DNA.letters[index] for index in indices)

    original_sequence_indices = np.random.choice(len(DNA.letters), genome_size)
    # Choose which positions to mutate and shift those positions by 1 to 3 slots
    mutation_index_shift = np.random.binomial(1, difference_rate, genome_size) * np.random.randint(1, 4, size=genome_size)

    mutated_sequence_indices = (original_sequence_indices + mutation_index_shift) % len(DNA.letters)

    original_sequence = indices_to_strand(original_sequence_indices)
    original_record = SeqRecord(
        Seq(original_sequence, DNA),
        id=f'Original:size={genome_size}',
    )
    mutated_sequence = indices_to_strand(mutated_sequence_indices)
    mutated_record = SeqRecord(
        Seq(mutated_sequence, DNA),
        id=f'Mutated:size={genome_size}'
    )

    SeqIO.write(original_record, f'sequence_data/original.fasta', 'fasta')
    SeqIO.write(mutated_record, f'sequence_data/mutated.fasta', 'fasta')

    return original_record, mutated_record


parser = argparse.ArgumentParser()
# Reference genomes parameters
parser.add_argument('-i', dest='input_files', type=str, nargs=2, help='Paths to files with reference genomes', required=False)
parser.add_argument('-n', dest='genome_size', nargs='?', type=int, help='Size of the two reference genomes')
parser.add_argument('-d', dest='difference_rate', nargs='?', type=float, help='Difference rate between the two sequences')

# Read parameters
parser.add_argument('-b', dest='backend', default='Custom', choices=[b[0] for b in backends], help='Read generation backend: Either Art Illumina or Custom')
parser.add_argument('-c', dest='coverage', nargs='?', default=30, type=int, help='Coverage in reads for one sequence')
parser.add_argument('-r', dest='read_length', nargs='?', default=150, type=int, help='Read length')

args = parser.parse_args()

if not args.input_files and not all([args.genome_size, args.difference_rate]):
    raise ValueError('Either specify 2 input files or parameters of the genome for mutation')

if args.input_files:
    sequences = [SeqIO.read(filename, format=get_file_type(filename)) for filename in args.input_files]
else:
    original, mutated = generate_mutated_sequence_pair(args.genome_size, args.difference_rate)
    sequences = [original, mutated]

backend: ReadGeneratorBackend = get_backend(args.backend)
for seq in sequences:
    prefix = f'read_data/{seq.id}_{args.backend}_{args.read_length}b_{args.coverage}x'
    read_count = backend.reads_from_sequence(
        seq,
        prefix,
        args.coverage,
        args.read_length
    )
    meta = GenomeReadsMetaData(read_length=args.read_length, coverage=args.coverage, num_of_reads=read_count, genome_size=len(seq))
    with open(f'{prefix}_meta.json', 'w') as f:
        f.write(str(meta))
