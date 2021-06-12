#!/usr/bin/env python
# coding: utf-8

import collections
import configparser
import glob
import itertools
import os
import platform
import re
import sys
import subprocess
from unittest import mock

import setuptools
from distutils import log
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants --------------------------------------------------------------

PLATFORM_MACHINE   = platform.machine()
SYS_IMPLEMENTATION = sys.implementation.name


# --- Utils ------------------------------------------------------------------

def _split_multiline(value):
    value = value.strip()
    sep = max('\n,;', key=value.count)
    return list(filter(None, map(lambda x: x.strip(), value.split(sep))))


# --- `setup.py` commands ----------------------------------------------------

class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.force = self.force
        self._clib_cmd.debug = self.debug

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {"include_path": ["include", "pyfastani"], "compiler_directives": {}}
        cython_args["compile_time_env"] = {"SYS_IMPLEMENTATION": SYS_IMPLEMENTATION}
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
            cython_args["compiler_directives"]["profile"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False
            cython_args["compiler_directives"]["cdivision"] = True

        # cythonize and patch the extensions
        self.extensions = cythonize(self.extensions, **cython_args)
        for ext in self.extensions:
            ext._needs_stub = False

        # # update the compiler include and link dirs to use the
        # # temporary build folder so that the platform-specific headers
        # # and static libs can be found

        # check the libraries have been built already
        if not self.distribution.have_run["build_clib"]:
            self._clib_cmd.run()

        # build the extensions as normal
        _build_ext.run(self)

    def build_extension(self, ext):
        # update compile flags if compiling in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-Og")
                ext.extra_compile_args.append("--coverage")
                ext.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class clean(_clean):

    def run(self):

        source_dir_abs = os.path.join(os.path.dirname(__file__), "pyfastani")
        source_dir = os.path.relpath(source_dir_abs)

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                log.info("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- Cython extensions ------------------------------------------------------

extensions = [
    Extension(
        "pyfastani._fastani",
        [os.path.join("pyfastani", "_utils.cpp"), os.path.join("pyfastani", "_fastani.pyx")],
        language="c++",
        libraries=["gsl", "gslcblas", "stdc++", "z", "m"],
        include_dirs=["include", "pyfastani", os.path.join("vendor", "FastANI", "src")],
        extra_compile_args=["-std=c++11"],
        extra_link_args=["-std=c++11"],
    )
]

# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=extensions,
    # libraries=libraries,
    cmdclass=dict(
        build_ext=build_ext,
        clean=clean,
        sdist=sdist,
    ),
)
