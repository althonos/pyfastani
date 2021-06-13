# 🐍⏩🧬 pyFastANI [![Stars](https://img.shields.io/github/stars/althonos/pyfastani.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyfastani/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [FastANI](https://github.com/ParBLiSS/FastANI/), a method for fast whole-genome similarity estimation.*

[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyfastani/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyfastani.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyfastani/issues)

## 🗺️ Overview

FastANI is a method published in 2018 by Jain *et al.* for high-throughput
computATION of whole-genome [Average Nucleotide Identity (ANI)](https://img.jgi.doe.gov/docs/ANI.pdf).
It uses [MashMap](https://github.com/marbl/MashMap) to compute orthologous mappings
without the need for expensive alignments.


`pyfastani` is a Python module, implemented using the [Cython](https://cython.org/)
language, that provides bindings to FastANI. It directly interacts with the
FastANI internals, which has the following advantages over CLI wrappers:

- **single dependency**: If your software or your analysis pipeline is
  distributed as a Python package, you can add `pyfastani` as a dependency to
  your project, and stop worrying about the FastANI binaries being properly
  setup on the end-user machine.
- **no intermediate files**: Everything happens in memory, in Python objects
  you have control on, making it easier to pass your sequences to FastANI
  without needing to write them to a temporary file.
- **no input formatting**: The memory model of FastANI can be interacted with
  to add new sequences directly without needing to write them to dedicated
  files. This is useful if your sequences are already loaded in memory,
  for instance because you obtained them from another Python library (such as
  or [Biopython](https://biopython.org/)).

*This library is still a work-in-progress, and in an experimental stage,
but it should already pack enough features to run one-to-one computations.*


## 💡 Example

The following snippets show how to compute the ANI between two genomes,
with the reference being a draft genome. For one-to-many or many-to-many
searches, simply add additional references with `m.add_draft` before indexing.
*Note that any name can be given to the reference sequences, this will just
affect the `name` attribute of the hits returned for a query.*

### 🔬 [Biopython](https://github.com/biopython/biopython)

Biopython does not let us access to the sequence directly, so we need to
convert it to bytes first with the `bytes` builtin function. For older
versions of Biopython (earlier than 1.79), use `record.seq.encode()`
instead of `bytes(record.seq).`

```python
import pyfastani
import Bio.SeqIO

m = pyfastani.Mapper()

# add a single draft genome to the mapper, and index it
ref = list(Bio.SeqIO.parse("vendor/FastANI/data/Shigella_flexneri_2a_01.fna", "fasta"))
m.add_draft("Shigella_flexneri_2a_01", (bytes(record.seq) for record in ref))
m.index()

# read the query and query the mapper
query = Bio.SeqIO.read("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna", "fasta")
hits = m.query_sequence(bytes(query.seq))

for hit in hits:
    print("Escherichia_coli_str_K12_MG1655", hit.name, hit.identity, hit.matches, hit.fragments)
```

### 🧪 [Scikit-bio](https://github.com/biocore/scikit-bio)

Scikit-bio lets us access to the sequence directly as a `numpy` array, but
shows the values as byte strings by default. To make them readable as
`char` (for compatibility with the C code), they must be cast with
`seq.values.view('B')`.

```python
import pyfastani
import skbio.io

m = pyfastani.Mapper()

ref = list(skbio.io.read("vendor/FastANI/data/Shigella_flexneri_2a_01.fna", "fasta"))
m.add_draft("Shigella_flexneri_2a_01", (seq.values.view('B') for seq in ref))
m.index()

# read the query and query the mapper
query = next(skbio.io.read("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna", "fasta"))
hits = m.query_genome(query.values.view('B'))

for hit in hits:
    print("Escherichia_coli_str_K12_MG1655", hit.name, hit.identity, hit.matches, hit.fragments)
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pyFastANI/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyFastANI/blob/master/CONTRIBUTING.md)
for more details.


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The FastANI code was written by [Chirag Jain](https://github.com/cjain7)
and is distributed under the terms of the
[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/) license,
unless otherwise specified from vendored sources.
See `vendor/FastANI/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original FastANI authors](https://github.com/cjain7). It was developed by
[Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*