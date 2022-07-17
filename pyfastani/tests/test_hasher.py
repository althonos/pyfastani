import os
import pickle
import sys
import unittest
import warnings

import pyfastani

class _TestHasher(object):

    _Hasher = None

    def test_hash_block_dimensions(self):
        data = b"ATGC" * 4

        k = 16
        hasher = self._Hasher(k)
        hashes = hasher.hash_block(data)
        self.assertEqual(len(hashes), 1)

        data += b"ATGC"
        hashes = hasher.hash_block(data)
        self.assertEqual(len(hashes), 5)

    def test_hash_block_consistency_8(self):
        data = b"ATGC" * 8

        k = 8
        hasher = self._Hasher(k)
        hashes = hasher.hash_block(data)

        for i, h in enumerate(hashes):
            self.assertEqual(hex(h), hex(hasher.hash(data[i:i+k])))

    def test_hash_block_consistency_16(self):
        data = b"ATGC" * 8

        k = 16
        hasher = self._Hasher(k)
        hashes = hasher.hash_block(data)

        for i, h in enumerate(hashes):
            self.assertEqual(hex(h), hex(hasher.hash(data[i:i+k])))


class TestMurmur3Hasher(_TestHasher, unittest.TestCase):

    _Hasher = pyfastani._fastani.Murmur3Hasher


class TestCuteHasher(_TestHasher, unittest.TestCase):

    _Hasher = pyfastani._fastani.CuteHasher

    # def test_reverse_hash(self):
    #     fwd = b"ATGC" * 4
    #     bwd =  b"GCAT" * 4
    #
    #     k = 16
    #     hasher = self._Hasher(k)
    #     h_fwd = hasher.hash(fwd)
    #     h_bwd = hasher.hash(bwd)
    #
    #     h_fwd_invert  = pyfastani._fastani._invert_hash(h_fwd)
    #     h_fwd_invert ^= 0xAAAAAAAA
    #
    #     self.assertEqual(h_bwd, h_fwd_invert, f"{h_bwd:032b} != {h_fwd_invert:032b}")
