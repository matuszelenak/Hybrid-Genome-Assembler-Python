import json
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Set, Generator, Tuple

BASES = ['A', 'C', 'G', 'T']


@dataclass
class GenomeRead:
    label: str
    category: str
    sequence: str

    def iterate_kmers(self, k: int) -> Generator[str, None, None]:
        for start in range(len(self.sequence) - (k - 1)):
            yield self.sequence[start:start + k]

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
    categories: Tuple[str, str] = ('A', 'B')

    def __str__(self):
        return json.dumps(dict(
            genome_size=self.genome_size,
            coverage=self.coverage,
            num_of_reads=self.num_of_reads,
            read_length=self.read_length,
            difference=self.difference,
            categories=self.categories
        ), indent=4)

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


@dataclass
class GenomeReadData:
    label: str
    category: str
    characteristic_kmers: Set[str]

    @property
    def characteristic_kmers_count(self):
        return len(self.characteristic_kmers)

    def __hash__(self):
        return self.label

    def __repr__(self):
        return f'{self.label}: {self.characteristic_kmers_count} characteristic kmers'

    def __copy__(self):
        return GenomeReadData(label=self.label, category=self.category, characteristic_kmers=self.characteristic_kmers.copy())

    @classmethod
    def from_json(cls, data) -> 'GenomeReadData':
        return cls(
            label=data['label'],
            category=data['category'],
            characteristic_kmers=set(data['characteristic_kmers'])
        )

    def to_json(self) -> dict:
        return dict(
            label=self.label,
            category=self.category,
            characteristic_kmers=list(self.characteristic_kmers)
        )


@dataclass
class GenomeReadCluster:
    reference_id: int
    reads: List[GenomeReadData]
    characteristic_kmers: Set[str]

    @property
    def size(self):
        return len(self.reads)

    @property
    def consistency(self):
        category_counts = defaultdict(int)
        for read in self.reads:
            category_counts[read.category] += 1
        return '/'.join(map(str, category_counts.values()))

    @property
    def categories(self) -> Set[str]:
        return set([read.category for read in self.reads])

    def copy(self):
        return GenomeReadCluster(
            reference_id=self.reference_id,
            characteristic_kmers=self.characteristic_kmers.copy(),
            reads=self.reads[::]
        )

    def __ior__(self, other: 'GenomeReadCluster'):
        self.reads.extend(other.reads)
        self.characteristic_kmers |= other.characteristic_kmers
        return self

    def __or__(self, other: 'GenomeReadCluster'):
        copy = self.copy()
        copy |= other
        return copy

    def __eq__(self, other: 'GenomeReadCluster'):
        return self.reference_id == other.reference_id

    def __hash__(self) -> int:
        return self.reference_id

    def __repr__(self):
        return f'#{self.reference_id}({len(self.characteristic_kmers)})[{self.consistency}]'
