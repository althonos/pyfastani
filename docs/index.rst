PyFastANI |Stars|
=================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyfastani.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyfastani/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `FastANI <ttps://github.com/ParBLiSS/FastANI/>`_,
*a method for fast whole-genome similarity estimation.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/workflow/status/althonos/pyfastani/Test/main?logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyfastani/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyfastani/branch/main.svg?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyfastani/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyfastani.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyfastani

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyfastani?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyfastani

.. |AUR| image:: https://img.shields.io/aur/version/python-pyfastani?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyfastani

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyfastani?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfastani/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyfastani.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfastani/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyfastani.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyfastani/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfastani/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyfastani/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyfastani.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyfastani/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyfastani?style=flat-square&maxAge=3600
   :target: http://pyfastani.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfastani/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyfastani
   :target: https://pepy.tech/project/pyfastani


Overview
--------

FastANI is a method published in 2018 by Jain *et al.* for high-throughput
computation of whole-genome `Average Nucleotide Identity (ANI) <https://img.jgi.doe.gov/docs/ANI.pdf>`_.
It uses `MashMap <https://github.com/marbl/MashMap>`_ to compute orthologous mappings
without the need for expensive alignments.

``pyfastani`` is a Python module, implemented using the `Cython <https://cython.org/>`_
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


Setup
-----

Run ``pip install pyfastani`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <install>` to find other ways to install ``pyfastani``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


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

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* ``fastANI`` *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
