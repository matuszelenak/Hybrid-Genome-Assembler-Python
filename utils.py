import math
from typing import Iterable

import mmh3
from bitarray import bitarray


class BloomFilter:
    def __init__(self, items_count, fp_prob):
        self.items_count = items_count
        self.fp_prob = fp_prob

        self.fp_prob = fp_prob
        self.size = self.get_size(items_count, fp_prob)
        self.hash_count = self.get_hash_count(self.size, items_count)
        self.bit_array = bitarray(self.size)
        self.bit_array.setall(0)

    def add(self, item):
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            self.bit_array[digest] = True

    def check(self, item):
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if not self.bit_array[digest]:
                return False
        return True

    @classmethod
    def get_size(cls, n, p):
        m = -(n * math.log(p)) / (math.log(2) ** 2)
        return int(m)

    @classmethod
    def get_hash_count(cls, m, n):
        k = (m / n) * math.log(2)
        return int(k)

    def __len__(self):
        return round((self.size / self.hash_count) * math.log(self.size / self.bit_array.count(False)))

    def __iand__(self, other: 'BloomFilter'):
        assert (self.bit_array.length(), self.hash_count) == (other.bit_array.length(), other.hash_count)
        self.bit_array &= other.bit_array
        return self

    def __and__(self, other: 'BloomFilter'):
        assert (self.bit_array.length(), self.hash_count) == (other.bit_array.length(), other.hash_count)
        result = BloomFilter(self.items_count, self.fp_prob)
        result.bit_array = self.bit_array & other.bit_array
        return result

    def __ior__(self, other: 'BloomFilter'):
        assert (self.bit_array.length(), self.hash_count) == (other.bit_array.length(), other.hash_count)
        self.bit_array |= other.bit_array
        return self

    def __or__(self, other: 'BloomFilter'):
        assert (self.bit_array.length(), self.hash_count) == (other.bit_array.length(), other.hash_count)
        result = BloomFilter(self.items_count, self.fp_prob)
        result.bit_array = self.bit_array | other.bit_array
        return result

    def __repr__(self):
        return f'{self.bit_array.to01()}'

    def copy(self):
        f = BloomFilter(self.items_count, self.fp_prob)
        f.bit_array = self.bit_array.copy()
        return f

    def set_content(self, bitstring: str):
        self.bit_array = bitarray(bitstring)

    def bulk_add(self, items):
        for item in items:
            self.add(item)


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
