import unittest

import Bio.SeqIO
import pyfastani


class TestMapper(unittest.TestCase):

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

        reference = list(Bio.SeqIO.parse("vendor/FastANI/data/Shigella_flexneri_2a_01.fna", "fasta"))
        mapper.add_draft(reference[0].id, (str(r.seq) for r in reference))
        mapper.index()

        query = Bio.SeqIO.read("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna", "fasta")
        hits = mapper.query_genome(str(query.seq))

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].name, reference[0].id)
        self.assertEqual(hits[0].matches, 1322)
        self.assertEqual(hits[0].fragments, 1547)
        self.assertAlmostEqual(hits[0].identity, 97.664, places=4)
