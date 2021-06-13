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


class publicize_headers(_build_ext):
    """A custom command to make everything in C++ headers public.
    """

    user_options = _build_ext.user_options + [
        ("sources=", "s", "the folder with the source headers to patch")
    ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.sources = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
        if self.sources is not None:
            self.sources = _split_multiline(self.sources)

    def run(self):
        # patch sources
        for source_dir in self.sources:
            for dirpath, dirnames, filenames in os.walk(source_dir):
                base_dirpath = os.path.relpath(dirpath, start=source_dir)
                new_dirpath = os.path.join( self.build_temp, base_dirpath )
                for filename in filenames:
                    if filename.endswith((".h", ".hpp")):
                        filepath = os.path.join(dirpath, filename)
                        new_filepath = os.path.join( new_dirpath, filename )
                        self.make_file([filepath], new_filepath, self.patch_header, (filepath, new_filepath))

        # update the include dirs
        for ext in self.extensions:
            ext.include_dirs.append(self.build_temp)

    def patch_header(self, old_path, new_path):
        self.mkpath(os.path.dirname(new_path))
        with open(old_path, "r") as src:
            source = src.read()

        # make everything public!
        source = re.sub("private:", "public:", source)
        # patch some classes with private members at the
        # beginning of their declaration (skch::Sketch)
        source = re.sub(r"class (.+\s*)\{(\s*)", r"class \1 {\npublic:\n\2", source)
        # patch the contructors of skch::Map and skch::Seq so that they don't
        # do anything
        source = re.sub(r"(Sketch|Map)\(([^)]*)\)([^}]+)\{[^}]*}", r"\1(\2)\3 {}", source)

        with open(new_path, "w") as dst:
            dst.write(source)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    # --- Autotools-like helpers ---

    def _silent_spawn(self, cmd):
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            raise CompileError(err.stderr)

    def _needs_clang_flags(self):
        self.mkpath(self.build_temp)
        testfile = os.path.join(self.build_temp, "testc++11.cpp")

        with open(testfile, "w") as f:
            f.write('#include <chrono>\n')
        try:
            with mock.patch.object(self.compiler, "spawn", new=self._silent_spawn):
                objects = self.compiler.compile([testfile], debug=self.debug)
        except CompileError as err:
            log.warn('failed to include <chrono>, assuming we need clang flags')
            return True
        else:
            log.info('successfully built a C++11 program with default flags')
            return False

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.force = self.force
        self._clib_cmd.debug = self.debug

    def run(self):
        # make sure sources have been patched to expose private fields
        if not self.distribution.have_run.get("publicize_headers", False):
            pub_cmd = self.get_finalized_command("publicize_headers")
            pub_cmd.force = self.force
            pub_cmd.run()

        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {"include_path": ["include", "pyfastani"], "compiler_directives": {}}
        cython_args["compile_time_env"] = {"FASTANI_PRIVATE_ACCESS": 1}
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

        # make sure to build with C++11
        if self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c++11")
            ext.extra_link_args.append("/std:c++11")
        else:
            ext.extra_compile_args.append("-std=c++11")
            ext.extra_link_args.append("-std=c++11")

        # in case we are compiling with clang, make sure to use libstdc++
        if self.compiler.compiler_type == "unix" and sys.platform == "darwin":
            ext.extra_compile_args.append("-stdlib=libc++")
            ext.extra_link_args.append("-stdlib=libc++")

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
        [os.path.join("pyfastani", x) for x in ("_utils.cpp", "_fastani.pyx")],
        language="c++",
        include_dirs=["include", "pyfastani"],
        define_macros=[("USE_BOOST", 1)],
    )
]

# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=extensions,
    cmdclass=dict(
        build_ext=build_ext,
        publicize_headers=publicize_headers,
        clean=clean,
        sdist=sdist,
    ),
)
