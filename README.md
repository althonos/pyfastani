# 🐍⏩🧬 PyFastANI [![Stars](https://img.shields.io/github/stars/althonos/pyfastani.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyfastani/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [FastANI](https://github.com/ParBLiSS/FastANI/), a method for fast whole-genome similarity estimation.*

[![Actions](https://img.shields.io/github/workflow/status/althonos/pyfastani/Test/main.svg?logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyfastani/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyfastani/branch/main.svg?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyfastani/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/pyfastani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyfastani)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyfastani?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyfastani)
[![AUR](https://img.shields.io/aur/version/python-pyfastani?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyfastani)
[![Wheel](https://img.shields.io/pypi/wheel/pyfastani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyfastani/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyfastani.svg?style=flat-square&maxAge=600)](https://pypi.org/project/pyfastani/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyfastani.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyfastani/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyfastani/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pyfastani/)
[![Issues](https://img.shields.io/github/issues/althonos/pyfastani.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyfastani/issues)
[![Docs](https://img.shields.io/readthedocs/pyfastani/latest?style=flat-square&maxAge=600)](https://pyfastani.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyfastani/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyfastani)](https://pepy.tech/project/pyfastani)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.4940237-purple?style=flat-square&maxAge=86400)](https://doi.org/10.5281/zenodo.4940237)


## 🗺️ Overview

FastANI is a method published in 2018 by Jain *et al.* for high-throughput
computation of whole-genome [Average Nucleotide Identity (ANI)](https://img.jgi.doe.gov/docs/ANI.pdf).
It uses [MashMap](https://github.com/marbl/MashMap) to compute orthologous mappings
without the need for expensive alignments.


`pyfastani` is a Python module, implemented using the [Cython](https://cython.org/)
language, that provides bindings to FastANI. It directly interacts with the
FastANI internals, which has the following advantages over CLI wrappers:

- **simpler compilation**: FastANI requires several additional libraries,
  which make compilation of the original binary non-trivial. In PyFastANI,
  libraries that were needed for threading or I/O are provided as stubs,
  and `Boost::math` headers are vendored so you can build the package without
  hassle. Or even better, just install from one of the provided wheels!
- **single dependency**: If your software or your analysis pipeline is
  distributed as a Python package, you can add `pyfastani` as a dependency to
  your project, and stop worrying about the FastANI binary being present on
  the end-user machine.
- **sans I/O**: Everything happens in memory, in Python objects you control,
  making it easier to pass your sequences to FastANI
  without needing to write them to a temporary file.

*This library is still a work-in-progress, and in an experimental stage,
but it should already pack enough features to be used in a standard pipeline.*


## 🔧 Installing

PyFastANI can be installed directly from [PyPI](https://pypi.org/project/pyfastani/),
which hosts some pre-built CPython wheels for x86-64 Unix platforms, as well
as the code required to compile from source with Cython:
```console
$ pip install pyfastani
```

In the event you have to compile the package from source, all the required
libraries are vendored in the source distribution, so you'll only need a
C/C++ compiler.

Otherwise, PyFastANI is also available as a [Bioconda](https://pyfastani.github.io/)
package:
```console
$ conda install -c bioconda pyfastani
```

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
instead of `bytes(record.seq)`.

```python
import pyfastani
import Bio.SeqIO

sketch = pyfastani.Sketch()

# add a single draft genome to the mapper, and index it
ref = list(Bio.SeqIO.parse("vendor/FastANI/data/Shigella_flexneri_2a_01.fna", "fasta"))
sketch.add_draft("S. flexneri", (bytes(record.seq) for record in ref))

# index the sketch and get a mapper
mapper = sketch.index()

# read the query and query the mapper
query = Bio.SeqIO.read("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna", "fasta")
hits = mapper.query_sequence(bytes(query.seq))

for hit in hits:
    print("E. coli K12 MG1655", hit.name, hit.identity, hit.matches, hit.fragments)
```

### 🧪 [Scikit-bio](https://github.com/biocore/scikit-bio)

Scikit-bio lets us access to the sequence directly as a `numpy` array, but
shows the values as byte strings by default. To make them readable as
`char` (for compatibility with the C code), they must be cast with
`seq.values.view('B')`.

```python
import pyfastani
import skbio.io

sketch = pyfastani.Sketch()

ref = list(skbio.io.read("vendor/FastANI/data/Shigella_flexneri_2a_01.fna", "fasta"))
sketch.add_draft("Shigella_flexneri_2a_01", (seq.values.view('B') for seq in ref))

mapper = sketch.index()

# read the query and query the mapper
query = next(skbio.io.read("vendor/FastANI/data/Escherichia_coli_str_K12_MG1655.fna", "fasta"))
hits = mapper.query_genome(query.values.view('B'))

for hit in hits:
    print("E. coli K12 MG1655", hit.name, hit.identity, hit.matches, hit.fragments)
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pyfastani/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyfastani/blob/master/CONTRIBUTING.md)
for more details.


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).

The FastANI code was written by [Chirag Jain](https://github.com/cjain7)
and is distributed under the terms of the
[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/),
unless otherwise specified in vendored sources. See `vendor/FastANI/LICENSE`
for more information.
The `cpu_features` code was written by [Guillaume Chatelet](https://github.com/gchatelet)
and is distributed under the terms of the [Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/).
See `vendor/cpu_features/LICENSE` for more information.
The `Boost::math` headers were written by [Boost Libraries](https://www.boost.org/) contributors
and is distributed under the terms of the [Boost Software License](https://choosealicense.com/licenses/bsl-1.0/).
See `vendor/boost-math/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original FastANI authors](https://github.com/cjain7). It was developed by
[Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
