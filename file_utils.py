import json
import os
from typing import Iterable, Generator

from structures import GenomeMetaData, GenomeRead, GenomeReadData


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
            category = label.split('-')[0]
            if not (label and sequence):
                return

            yield GenomeRead(label=label, sequence=sequence, category=category)


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
        for read_json in j['reads']:
            yield GenomeReadData.from_json(read_json)


class GenomeReadDataWriter:
    def __init__(self, filename: str):
        self.filename = filename

    def __enter__(self):
        self.file = open(self.filename, 'w')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def write(self, read_data: Iterable[GenomeReadData]):
        json.dump({
            'reads': [
                read.to_json() for read in read_data
            ]
        }, self.file)


def get_file_type(filename):
    _, ext = os.path.splitext(filename)
    if ext.lower() in ('.fasta', '.fa'):
        return 'fasta'
    if ext.lower() in ('.fq', '.fastq'):
        return 'fastq'
    return None
