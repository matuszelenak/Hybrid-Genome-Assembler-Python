import json
from typing import Iterable, Generator

from structures import GenomeMetaData, GenomeRead, GenomeReadData
from utils import BloomFilter


class GenomeReader:
    meta: GenomeMetaData

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.file = open(self.filename)
        self.meta = GenomeMetaData.from_string(self.file.readline().strip())
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def get_reads(self) -> Generator[GenomeRead, None, None]:
        self.file.seek(0)
        self.file.readline()

        while True:
            label = self.file.readline().strip()[1:]
            sequence = self.file.readline().strip()
            if not (label and sequence):
                return

            yield GenomeRead(label=label, sequence=sequence)


class GenomeWriter:
    def __init__(self, filename: str, meta: GenomeMetaData):
        self.filename = filename
        self.meta = meta

    def __enter__(self):
        self.file = open(self.filename, 'w')
        self.file.write(str(self.meta))
        self.file.write('\n')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def write(self, reads: Iterable[GenomeRead]):
        for read in reads:
            self.file.write(str(read))
            self.file.write('\n')


class GenomeReadDataReader:
    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.file = open(self.filename)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def get_read_data(self) -> Generator[GenomeReadData, None, None]:
        j = json.load(self.file)
        bloom_fp_prob = j['bloom_config']['fp_prob']
        bloom_items_count = j['bloom_config']['items_count']
        for read_json in j['reads']:
            obj = GenomeReadData.from_json(read_json)
            obj.characteristic_kmers = BloomFilter(bloom_items_count, bloom_fp_prob)
            obj.characteristic_kmers.set_content(read_json['characteristic_kmers'])

            yield obj


class GenomeReadDataWriter:
    def __init__(self, filename: str):
        self.filename = filename

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def write(self, read_data: Iterable[GenomeReadData]):
        bloom_filter = next(iter(read_data)).characteristic_kmers
        j = {
            'bloom_config': {
                'fp_prob': bloom_filter.fp_prob,
                'items_count': bloom_filter.items_count
            },
            'reads': [
                read.to_json() for read in read_data
            ]
        }
        json.dump(j, self.file, indent=4)
