import os
import pickle
from functools import lru_cache
from typing import Iterable, List, Generator

from Bio.SeqRecord import SeqRecord


def minimum_base_quality(minimum: int, qualities: List[int]):
    if min(qualities) < minimum:
        return False

    return True


def minimum_and_average(minimum, average, qualities: List[int]):
    if min(qualities) < minimum or sum(qualities) / len(qualities) < average:
        return False
    return True


@lru_cache(maxsize=100000)
def sequence_complement(sequence: str):
    c = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get
    return ''.join(map(c, reversed(sequence)))


def iterate_kmers_in_record(record: SeqRecord, k: int, quality_filter=None):
    sequence = str(record.seq)
    if callable(quality_filter):
        qualities: List[int] = record.letter_annotations['phred_quality']
        for start in range(len(sequence) - (k - 1)):
            kmer_quality = qualities[start:start + k]
            if quality_filter(kmer_quality):
                yield sequence[start:start + k]
    else:
        for start in range(len(sequence) - (k - 1)):
            yield sequence[start:start + k]


def iterate_kmer_signatures(record: SeqRecord, k: int, quality_filter=None):
    for kmer in iterate_kmers_in_record(record, k, quality_filter=quality_filter):
        yield min(kmer, sequence_complement(kmer))


def update_progress(progress, total):
    print("\rProgress: [{0:50s}] {1:.1f}% ({2}/{3})".format('#' * int(progress / total * 50), progress / total * 100, progress, total), end="", flush=True)
    if progress == total:
        print()


def iter_with_progress(iterable: Iterable, total_length: int = None, start_message: str = None):
    if not total_length:
        iterable = list(iterable)

    if start_message:
        print(start_message)

    for processed, item in enumerate(iterable):
        yield item
        update_progress(processed + 1, total_length or len(iterable))


def split_into_chunks(iterable: list, chunk_length: int):
    for beginning in range(0, len(iterable), chunk_length):
        yield iterable[beginning:beginning + chunk_length]


def with_category(*reads_by_category: Iterable[SeqRecord]) -> Generator[SeqRecord, None, None]:
    for category_id, category_reads in enumerate(reads_by_category):
        for record in category_reads:
            record.category_id = category_id
            yield record


class Cache:
    global_kwargs = {}

    def get_key_from_arguments(self, *args, function_specific_kwargs=None, **kwargs):
        function_specific_kwargs = function_specific_kwargs or {}
        primitives = (int, str, bool, float)

        valid_args = []
        valid_kwargs = {}
        for arg in args:
            if isinstance(arg, primitives):
                valid_args.append(arg)

        for kw, value in kwargs.items():
            if isinstance(value, primitives):
                valid_kwargs[kw] = value

        valid_kwargs.update({
            kw: str(val) for kw, val in function_specific_kwargs.items()
        })

        valid_kwargs.update({
            kw: str(val) for kw, val in self.global_kwargs.items()
        })

        return ';'.join(map(str, valid_args)) + '|' + ';'.join(f'{kw}={str(valid_kwargs[kw])}' for kw in sorted(valid_kwargs.keys()))

    def checkpoint(self, **cache_kwargs):
        def decorator(function):
            def wrapper(*args, **kwargs):
                key = self.get_key_from_arguments(*args, function_specific_kwargs=cache_kwargs, **kwargs)
                path = os.path.join('cache', function.__name__ + ':' + key)
                if os.path.exists(path):
                    with open(path, 'rb') as f:
                        return pickle.load(f)

                result = function(*args, **kwargs)

                if not os.path.exists('cache'):
                    os.mkdir('cache')

                with open(path, 'wb') as f:
                    pickle.dump(result, f)

                return result

            return wrapper

        return decorator
