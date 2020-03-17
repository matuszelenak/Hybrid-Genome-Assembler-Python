import argparse
import math
from collections import Counter, defaultdict
from typing import List, Generator, Dict, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna as DNA

from utils import iterate_kmers_in_record


def get_expected_characteristic_kmers(sequence_a: SeqRecord, sequence_b: SeqRecord, k: int) -> Tuple[List[str], List[str]]:
    total, diff = len(sequence_a), 0
    for base_a, base_b in zip(sequence_a.seq, sequence_b.seq):
        if base_a != base_b:
            diff += 1

    print(diff / total * 100)

    char_a, char_b = [], []
    for kmer_a, kmer_b in zip(iterate_kmers_in_record(sequence_a.seq, k), iterate_kmers_in_record(sequence_b.seq, k)):
        if kmer_a != kmer_b:
            char_a.append(kmer_a)
            char_b.append(kmer_b)

    return char_a, char_b


parser = argparse.ArgumentParser()
# Reference genomes parameters
parser.add_argument('input_files', type=str, nargs=2, help='Paths to files with reference genomes')
args = parser.parse_args()

seq_a = SeqIO.read(args.input_files[0], 'fasta')
seq_b = SeqIO.read(args.input_files[1], 'fasta')

k_guess = math.ceil(math.log(len(seq_a), 4)) + 3
print(k_guess)

expected_a, expected_b = get_expected_characteristic_kmers(seq_a, seq_b, k_guess)
c = Counter(expected_a).most_common(100)
print(len(expected_a), len(set(expected_a)))
print(c)

print('\n\n\n')
c = Counter(expected_b).most_common(100)
print(len(expected_b), len(set(expected_b)))
print(c)

print(len(set(expected_a).intersection(set(expected_b))))
