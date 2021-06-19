import os
import sys
import unittest

from .. import Sketch
from .._fasta import Parser

PROJECT_PATH = os.path.realpath(os.path.join(__file__, "..", "..", ".."))
FASTANI_PATH = os.path.join(PROJECT_PATH, "vendor", "FastANI")

ECOLI = os.path.join(FASTANI_PATH, "data", "Escherichia_coli_str_K12_MG1655.fna")
SFLEXNERI = os.path.join(FASTANI_PATH, "data", "Shigella_flexneri_2a_01.fna")

class _TestANI(object):

    def _load_fasta():
        pass

    @unittest.skipUnless(os.path.exists(ECOLI), "missing FastANI data files")
    @unittest.skipUnless(os.path.exists(SFLEXNERI), "missing FastANI data files")
    def test_fastani_example(self):
        """Check that we get the same results as FastANI on their example data.
        """
        # The example in the FastANI README is the following
        # $ ./fastANI -q data/Shigella_flexneri_2a_01.fna -r data/Escherichia_coli_str_K12_MG1655.fna
        # data/Shigella_flexneri_2a_01.fna data/Escherichia_coli_str_K12_MG1655.fna 97.7507 1303 1608

        sketch = Sketch()

        ref = self._load_fasta(ECOLI)
        sketch.add_genome("Escherichia_coli_str_K12_MG1655", self._get_sequence(ref[0]))

        mapper = sketch.index()

        contigs = self._load_fasta(SFLEXNERI)
        hits = mapper.query_draft(map(self._get_sequence, contigs))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Escherichia_coli_str_K12_MG1655")
        self.assertEqual(hits[0].matches, 1303)
        self.assertEqual(hits[0].fragments, 1608)
        self.assertAlmostEqual(hits[0].identity, 97.7507, places=4)


class TestANIString(_TestANI, unittest.TestCase):

    def _load_fasta(self, path):
        return list(Parser(path))

    def _get_sequence(self, record):
        return record.seq.decode("ascii")


class TestANIBytes(_TestANI, unittest.TestCase):

    def _load_fasta(self, path):
        return list(Parser(path))

    def _get_sequence(self, record):
        return record.seq


try:
    from skbio import io as skbio_io
except ImportError:
    skbio_io = None

@unittest.skipUnless(skbio_io, "Scikit-bio is required for this test suite")
class TestANISkbio(_TestANI, unittest.TestCase):

    def _load_fasta(self, path):
        return list(skbio_io.read(path, "fasta"))

    def _get_sequence(self, sequence):
        return sequence.values.view('B')

try:
    import Bio.SeqIO
except ImportError:
    Bio = None

@unittest.skipUnless(Bio, "Biopython is required for this test suite")
class TestANIBiopython(_TestANI, unittest.TestCase):

    def _load_fasta(self, path):
        return list(Bio.SeqIO.parse(path, "fasta"))

    def _get_sequence(self, record):
        version = tuple(map(int, Bio.__version__.split(".")))
        return record.seq.encode() if version < (1, 79) else bytes(record.seq)
