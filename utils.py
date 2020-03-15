from typing import Iterable

from Bio import Seq


def iterate_kmers_in_sequence(sequence: Seq, k: int):
    for start in range(len(sequence) - (k - 1)):
        yield sequence[start:start + k]


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
