import json
import os
from typing import Iterable, Generator

from structures import GenomeReadData


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
        }, self.file, indent=2)


def get_file_type(filename):
    _, ext = os.path.splitext(filename)
    if ext.lower() in ('.fasta', '.fa'):
        return 'fasta'
    if ext.lower() in ('.fq', '.fastq'):
        return 'fastq'
    return None
