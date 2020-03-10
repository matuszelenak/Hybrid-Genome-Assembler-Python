import json
from dataclasses import dataclass
from typing import List, Set, Optional, Generator

from utils import BloomFilter

BASES = ['A', 'C', 'G', 'T']


@dataclass
class Kmer:
    sequence: str

    def __str__(self):
        return f'Kmer {self.sequence}'

    def __repr__(self):
        return self.sequence

    def __hash__(self):
        return self.sequence.__hash__()


@dataclass
class GenomeRead:
    label: str
    sequence: str

    def iterate_kmers(self, k: int) -> Generator[Kmer, None, None]:
        for start in range(len(self.sequence) - (k - 1)):
            yield Kmer(self.sequence[start:start + k])

    def __str__(self):
        return f'>{self.label}\n{"".join(self.sequence)}'


@dataclass
class GenomeMetaData:
    genome_size: int
    num_of_reads: int
    coverage: int = 1
    difference: float = 0
    read_length: int = 0
    circular: bool = True

    def __str__(self):
        return json.dumps(dict(
            genome_size=self.genome_size,
            coverage=self.coverage,
            num_of_reads=self.num_of_reads,
            read_length=self.read_length,
            difference=self.difference
        ))

    @classmethod
    def from_string(cls, string):
        try:
            data = json.loads(string)
            return cls(**data)
        except json.JSONDecodeError as e:
            print(e)
            return
        except TypeError:
            return


class GenomeReadData:
    label: str
    characteristic_kmers: BloomFilter
    characteristic_kmers_count: int

    def __init__(self, label):
        self.label = label

    def index_characteristic_kmers(self, genome_read: GenomeRead, characteristic_kmer_set: Set[str], estimated_size: int) -> int:
        k: int = len(next(iter(characteristic_kmer_set)))
        self.characteristic_kmers = BloomFilter(estimated_size, 0.01)
        characteristic, total = 0, 0
        for kmer in genome_read.iterate_kmers(k):
            if kmer.sequence in characteristic_kmer_set:
                self.characteristic_kmers.add(kmer.sequence)
                characteristic += 1

            total += 1

        self.characteristic_kmers_count = characteristic
        return characteristic

    def __hash__(self):
        return self.label

    def __repr__(self):
        return f'{self.label}: {self.characteristic_kmers_count} characteristic kmers, {self.characteristic_kmers}'

    @classmethod
    def from_json(cls, data) -> 'GenomeReadData':
        obj = cls(data['label'])
        obj.characteristic_kmers_count = data['characteristic_kmer_count']
        return obj

    def to_json(self) -> dict:
        return dict(
            label=self.label,
            characteristic_kmer_count=self.characteristic_kmers_count,
            characteristic_kmers=self.characteristic_kmers.bit_array.to01() if self.characteristic_kmers else ''
        )


@dataclass
class GenomeReadCluster:
    reference_id: int
    characteristic_kmers: Optional[BloomFilter]
    reads: List[GenomeReadData]

    @property
    def size(self):
        return len(self.reads)

    def clear(self):
        self.characteristic_kmers = None
        self.reads = []

    @property
    def labels(self):
        return ','.join([read.label for read in self.reads])[:200]

    @property
    def stats(self):
        a_reads = len([x for x in self.reads if x.label[0] == 'A'])
        b_reads = len([x for x in self.reads if x.label[0] != 'A'])
        return f'#{self.reference_id} size {self.size} A {a_reads}/{b_reads} B'

    def __ior__(self, other: 'GenomeReadCluster'):
        self.reads.extend(other.reads)
        self.characteristic_kmers |= other.characteristic_kmers
        return self

    def __eq__(self, other: 'GenomeReadCluster'):
        return self.reference_id == other.reference_id

    def __hash__(self) -> int:
        return self.reference_id

    def __repr__(self):
        return f'Cluster #{self.reference_id} of size {self.size} kmers {self.characteristic_kmers}'
