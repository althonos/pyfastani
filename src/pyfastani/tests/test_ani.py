import os
import pickle
import sys
import unittest

from .. import Sketch
from .._fasta import Parser

# test data folder
DATA_FOLDER = os.path.realpath(os.path.join(__file__, "..", "data"))

# vendored test files from the FastANI sources
ECOLI = os.path.join(DATA_FOLDER, "Escherichia_coli_str_K12_MG1655.fna")
SFLEXNERI = os.path.join(DATA_FOLDER, "Shigella_flexneri_2a_01.fna")

# local test files from MIBiG
BGC0001425 = os.path.join(DATA_FOLDER, "BGC0001425.faa")
BGC0001427 = os.path.join(DATA_FOLDER, "BGC0001427.faa")
BGC0001428 = os.path.join(DATA_FOLDER, "BGC0001428.faa")


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
        sketch.add_draft(
            "Escherichia_coli_str_K12_MG1655",
            [self._get_sequence(r) for r in self._load_fasta(ECOLI)]
        )

        mapper = sketch.index()

        contigs = self._load_fasta(SFLEXNERI)
        hits = mapper.query_draft(map(self._get_sequence, contigs))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Escherichia_coli_str_K12_MG1655")
        self.assertEqual(hits[0].matches, 1303)
        self.assertEqual(hits[0].fragments, 1608)
        self.assertAlmostEqual(hits[0].identity, 97.7507, places=4)

    @unittest.skipUnless(os.path.exists(ECOLI), "missing FastANI data files")
    def test_escherichia_minimizers(self):
        """Check that we extract as many minimizers as FastANI on their data.
        """
        contigs = [self._get_sequence(r) for r in self._load_fasta(ECOLI)]

        sketch = Sketch()
        self.assertEqual(sketch.window_size, 24)
        sketch.add_draft("Escherichia_coli_str_K12_MG1655", contigs)
        self.assertEqual(len(sketch.minimizers), 371301)
        mapper = sketch.index()
        self.assertEqual(len(mapper.lookup_index), 361568)

        hits = mapper.query_draft(contigs)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Escherichia_coli_str_K12_MG1655")
        self.assertEqual(hits[0].matches, 1547)
        self.assertEqual(hits[0].fragments, 1547)
        self.assertAlmostEqual(hits[0].identity, 100.0)

    @unittest.skipUnless(os.path.exists(SFLEXNERI), "missing FastANI data files")
    def test_shigella_minimizers(self):
        """Check that we extract as many minimizers as FastANI on their data.
        """
        contigs = [self._get_sequence(r) for r in self._load_fasta(SFLEXNERI)]

        sketch = Sketch()
        self.assertEqual(sketch.window_size, 24)
        sketch.add_draft("Shigella_flexneri_2a_01", contigs)
        self.assertEqual(len(sketch.minimizers), 386387)
        mapper = sketch.index()
        self.assertEqual(len(mapper.lookup_index), 347908)

        hits = mapper.query_draft(contigs)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Shigella_flexneri_2a_01")
        self.assertEqual(hits[0].matches, 1600)
        self.assertEqual(hits[0].fragments, 1608)
        self.assertAlmostEqual(hits[0].identity, 100.0)

    @unittest.skipUnless(os.path.exists(BGC0001425), "missing test data files")
    @unittest.skipUnless(os.path.exists(BGC0001427), "missing test data files")
    @unittest.skipUnless(os.path.exists(BGC0001428), "missing test data files")
    def test_myxochromide_bgcs(self):
        """Check that we get expected hits between homologous BGCs.
        """
        sketch = Sketch(protein=True, fragment_length=100)
        bgc1 = self._load_fasta(BGC0001425)
        sketch.add_draft("BGC0001425", map(self._get_sequence, bgc1))
        bgc2 = self._load_fasta(BGC0001427)
        sketch.add_draft("BGC0001427", map(self._get_sequence, bgc1))

        mapper = sketch.index()
        bgc3 = self._load_fasta(BGC0001428)
        hits = mapper.query_draft(map(self._get_sequence, bgc3))

        self.assertEqual(len(hits), 2)
        self.assertEqual(hits[0].name, "BGC0001425")
        self.assertEqual(hits[0].matches, 130)
        self.assertEqual(hits[0].fragments, 176)
        self.assertEqual(hits[1].name, "BGC0001427")
        self.assertEqual(hits[1].matches, 130)
        self.assertEqual(hits[1].fragments, 176)


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

    def test_sketch_pickling(self):
        """Check that pickling before indexing produces consistent results.
        """
        sketch = Sketch()

        ref = self._load_fasta(ECOLI)
        sketch.add_genome("Escherichia_coli_str_K12_MG1655", self._get_sequence(ref[0]))

        mapper = pickle.loads(pickle.dumps(sketch)).index()

        contigs = self._load_fasta(SFLEXNERI)
        hits = mapper.query_draft(map(self._get_sequence, contigs))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Escherichia_coli_str_K12_MG1655")
        self.assertEqual(hits[0].matches, 1303)
        self.assertEqual(hits[0].fragments, 1608)
        self.assertAlmostEqual(hits[0].identity, 97.7507, places=4)

    def test_mapper_pickling(self):
        """Check that pickling after indexing produces consistent results.
        """
        sketch = Sketch()

        ref = self._load_fasta(ECOLI)
        sketch.add_genome("Escherichia_coli_str_K12_MG1655", self._get_sequence(ref[0]))

        mapper = pickle.loads(pickle.dumps(sketch.index()))

        contigs = self._load_fasta(SFLEXNERI)
        hits = mapper.query_draft(map(self._get_sequence, contigs))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Escherichia_coli_str_K12_MG1655")
        self.assertEqual(hits[0].matches, 1303)
        self.assertEqual(hits[0].fragments, 1608)
        self.assertAlmostEqual(hits[0].identity, 97.7507, places=4)


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
