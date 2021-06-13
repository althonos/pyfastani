import unittest

import pyfastani


class _TestMapper(object):

    def _load_fasta():
        pass

    def test_fastani_example(self):
        # The example in the FastANI README is the following
        # $ ./fastANI -q data/Shigella_flexneri_2a_01.fna -r data/Escherichia_coli_str_K12_MG1655.fna
        # data/Shigella_flexneri_2a_01.fna data/Escherichia_coli_str_K12_MG1655.fna 97.7507 1303 1608

        pass

    def test_fastani_example_reversed(self):
        # Same as the FastANI README example, but swapping query and reference
        # $ ./fastANI -r data/Shigella_flexneri_2a_01.fna -q data/Escherichia_coli_str_K12_MG1655.fna -o fastani.txt
        # $ cat fastani.txt
        # data/Escherichia_coli_str_K12_MG1655.fna	data/Shigella_flexneri_2a_01.fna	97.664	1322	1547

        mapper = pyfastani.Mapper()

        contigs = self._load_fasta("vendor/FastANI/data/Shigella_flexneri_2a_01.fna")
        mapper.add_draft("Shigella_flexneri_2a_01", map(self._get_sequence, contigs))
        mapper.index()

        query = self._load_fasta("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna")
        hits = mapper.query_genome(self._get_sequence(query[0]))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, "Shigella_flexneri_2a_01")
        self.assertEqual(hits[0].matches, 1322)
        self.assertEqual(hits[0].fragments, 1547)
        self.assertAlmostEqual(hits[0].identity, 97.664, places=4)


try:
    from skbio import io as skbio_io
except ImportError:
    skbio_io = None

@unittest.skipUnless(skbio_io, "Scikit-bio is required for this test suite")
class TestMapperSkbio(_TestMapper, unittest.TestCase):

    def _load_fasta(self, path):
        return list(skbio_io.read(path, "fasta"))

    def _get_sequence(self, sequence):
        return sequence.values.view('B')

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

@unittest.skipUnless(SeqIO, "Biopython is required for this test suite")
class TestMapperBiopython(_TestMapper, unittest.TestCase):

    def _load_fasta(self, path):
        return list(SeqIO.parse(path, "fasta"))

    def _get_sequence(self, sequence):
        return sequence.seq.encode()
