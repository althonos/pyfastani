PyFastANI |Stars|
=================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyfastani.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyfastani/stargazers
   :class: dark-light

`Cython <https://cython.org/>`_ *bindings and Python interface to* `FastANI <ttps://github.com/ParBLiSS/FastANI/>`_,
*a method for fast whole-genome similarity estimation.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyfastani/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyfastani/actions
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyfastani/branch/main.svg?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyfastani/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/pyfastani.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyfastani
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyfastani?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyfastani
   :class: dark-light

.. |AUR| image:: https://img.shields.io/aur/version/python-pyfastani?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyfastani
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyfastani?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfastani/#files
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyfastani.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfastani/#files
   :class: dark-light

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyfastani.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyfastani/#files
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfastani/
   :class: dark-light

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyfastani/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyfastani.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyfastani/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/pyfastani?style=flat-square&maxAge=3600
   :target: http://pyfastani.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfastani/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyfastani?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyfastani
   :class: dark-light


Overview
--------

FastANI is a method published in 2018 by Jain *et al.* for high-throughput
computation of whole-genome `Average Nucleotide Identity (ANI) <https://img.jgi.doe.gov/docs/ANI.pdf>`_.
It uses `MashMap <https://github.com/marbl/MashMap>`_ to compute orthologous mappings
without the need for expensive alignments.

``pyfastani`` is a Python module, implemented using the `Cython <https://cython.org/>`_
language, that provides bindings to FastANI. It directly interacts with the
FastANI internals, which has the following advantages over CLI wrappers:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyfastani`` as a ``pip`` or ``conda`` dependency, no need
      for the ``fastani`` binary or any external dependency.

   .. grid-item-card:: :fas:`hammer` Easy compilation

      Required libraries that were needed for threading or I/O are provided 
      as stubs, and `Boost::math` headers are vendored to build the package 
      without any system dependency.

   .. grid-item-card:: :fas:`file-circle-xmark` Sans I/O

      Everything happens in memory, making it easier to pass your sequences 
      to FastANI without needing to write them to a temporary file.

   .. grid-item-card:: :fas:`server` Multi-threaded

      Genome query resolves the fragment mapping step in parallel, leading to 
      shorter querying times even with a single genome.

   .. grid-item-card:: :fas:`microchip` Portable

      Get SIMD-acceleration on any supported platform without having
      to build the package from scratch.

   .. grid-item-card:: :fas:`magnifying-glass` Introspectable

      The genome sketches can be accessed from the Python API, allowing
      to view the minimizers for a genome database.


Setup
-----

PyFastANI is available for all modern Python versions (3.7+).

Run ``pip install pyfastani`` in a shell to download the latest release 
from PyPi, or have a look at the :doc:`Installation page <guide/install>` to find 
other ways to install ``pyfastani``.

Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. grid:: 1 3 5 5
   :gutter: 1

   .. grid-item-card:: :fas:`diamond` PyHMMER
      :link: https://pyhmmer.readthedocs.io

      Profile Hidden Markov Models (with HMMER).

   .. grid-item-card:: :fas:`fire` Pyrodigal
      :link: https://pyrodigal.readthedocs.io

      Prokaryotic Gene Finding (with Prodigal).

   .. grid-item-card:: :fas:`virus-covid` Pyrodigal-gv
      :link: https://github.com/althonos/pyrodigal-gv

      Pyrodigal for Giant Viruses.

   .. grid-item-card:: :fas:`align-center` PyFAMSA
      :link: https://pyfamsa.readthedocs.io

      Multiple Sequence Alignment (with FAMSA).

   .. grid-item-card:: :fas:`scissors` PytrimAl
      :link: https://pytrimal.readthedocs.io

      Alignment Trimming (with trimAl).

   .. grid-item-card:: :fas:`music` LightMotif
      :link: https://lightmotif.readthedocs.io

      Platform-accelerated motif scoring.

   .. grid-item-card:: :fas:`knife;fa-custom` Diced
      :link: https://diced.readthedocs.io

      CRISPR Detection (with MinCED).

   .. grid-item-card:: :fas:`table-cells` Scoring Matrices
      :link: https://scoring-matrices.readthedocs.io

      Scoring matrices for Cython.

   .. grid-item-card:: :fas:`chain` Pyskani
      :link: https://pyskani.readthedocs.io

      Average Nucleotide Identity (with skani).

   .. grid-item-card:: :fas:`forward-fast` PyFastANI
      :link: https://pyfastani.readthedocs.io

      Average Nucleotide Identity (with FastANI).

   .. grid-item-card:: :fas:`magnifying-glass` PyJess
      :link: https://pyjess.readthedocs.io

      Geometric Template Matching (with Jess).

   .. grid-item-card:: :fas:`repeat` PyTantan
      :link: https://pytantan.readthedocs.io

      Tandem Repeat Masking (with Tantan).

   .. grid-item-card:: :fas:`gem` PyOpal
      :link: https://pyopal.readthedocs.io

      Query/Database Aligner (with Opal).

   .. grid-item-card:: :fas:`sword;fa-custom` PySWRD
      :link: https://pyswrd.readthedocs.io

      Database Heuristic Filtering (with SWORD).

   .. grid-item-card:: :fas:`rocket` Mini3di
      :link: https://github.com/althonos/mini3di

      Protein structure to 3di in pure Python.

   .. grid-item-card:: :fas:`calculator` ``peptides.py``
      :link: https://peptides.readthedocs.io

      Peptide descriptors for Python.

   .. grid-item-card:: :fas:`diagram-project` Pronto
      :link: https://pronto.readthedocs.io

      Open Biomedical Ontologies for Python.

   .. grid-item-card:: :fas:`box` NAFcodec
      :link: https://nafcodec.readthedocs.io

      Nucleotide Archival Format for Python.

   .. grid-item-card:: :fas:`bank` ``gb-io.py``
      :link: https://gb-io.readthedocs.io

      Fast GenBank parser for Python (with ``gb-io``).



License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.

The fastANI code was written by `Chirag Jain <https://github.com/cjain7>`_
and is distributed under the terms of the
`Apache License 2.0 <https://choosealicense.com/licenses/apache-2.0/>`_,
unless otherwise specified in vendored sources. The ``cpu_features`` code was written by `Guillaume Chatelet <https://github.com/gchatelet>`_
and is distributed under the terms of the `Apache License 2.0 <https://choosealicense.com/licenses/apache-2.0/>`_.
The ``Boost::math`` headers were written by `Boost Libraries <https://www.boost.org/>`_ contributors
and is distributed under the terms of the `Boost Software License <https://choosealicense.com/licenses/bsl-1.0/>`_.
See the :doc:`Copyright page <guide/copyright>` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* ``fastANI`` *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
