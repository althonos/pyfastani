import os
import pickle
import sys
import unittest
import warnings

import pyfastani


class TestSketch(unittest.TestCase):

    def test_init_errors(self):
        """Check that constructor parameters are properly validated.
        """
        self.assertRaises(TypeError, pyfastani.Sketch, k="1")
        self.assertRaises(TypeError, pyfastani.Sketch, fragment_length="1")
        self.assertRaises(TypeError, pyfastani.Sketch, minimum_fraction="0.5")

        self.assertRaises(OverflowError, pyfastani.Sketch, k=2**32)
        self.assertRaises(ValueError, pyfastani.Sketch, k=0)
        self.assertRaises(ValueError, pyfastani.Sketch, p_value=-1.0)
        self.assertRaises(ValueError, pyfastani.Sketch, percentage_identity=-1.0)
        self.assertRaises(ValueError, pyfastani.Sketch, percentage_identity=200.0)

    def test_reinit(self):
        """Check that calling `__init__` more than once does not crash.
        """
        sketch = pyfastani.Sketch(fragment_length=100)
        sketch.add_genome("test", "ATGC"*100)
        self.assertEqual(sketch.names, ["test"])
        self.assertEqual(sketch.fragment_length, 100)

        sketch.__init__(fragment_length=200)
        self.assertEqual(sketch.names, [])
        self.assertEqual(sketch.fragment_length, 200)

    def test_add_draft_warnings(self):
        """Check that `Sketch.add_draft` raises warnings as expected.
        """
        sketch = pyfastani.Sketch()
        with warnings.catch_warnings(record=True) as catch:
            sketch.add_draft("short_seq", ["ATGC"*1000, "ATGC"])
            self.assertEqual(len(catch), 1)  # second sequence is too short

    def test_add_sequence_short(self):
        """Check that a sequence too short to be hashed is still recorded.
        """
        sketch = pyfastani.Sketch()
        with warnings.catch_warnings(record=True) as catch:
            warnings.simplefilter("ignore")
            sketch.add_draft("short_seq", ["ATGC"*1000, "ATGC"])
            self.assertEqual(sketch.names, ["short_seq"])

    def test_pickle(self):
        """Check that `Sketch` objects can be pickled and unpickled.
        """
        sketch = pyfastani.Sketch()
        sketch.add_genome("short_seq", "ATGC"*1000)

        sketch2 = pickle.loads(pickle.dumps(sketch))
        self.assertEqual(sketch2.names, ["short_seq"])
