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
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library as _Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants --------------------------------------------------------------

PLATFORM_MACHINE = platform.machine()
PLATFORM_SYSTEM  = platform.system()
UNIX = not sys.platform.startswith("win")

PROCESSOR_IS_MIPS = PLATFORM_MACHINE.startswith("mips")
PROCESSOR_IS_ARM = PLATFORM_MACHINE.startswith("arm")
PROCESSOR_IS_AARCH64 = PLATFORM_MACHINE.startswith("aarch64")
PROCESSOR_IS_X86 = PLATFORM_MACHINE.startswith(("x86_64", "AMD64", "amd64", "i386", "i486"))
PROCESSOR_IS_POWER = PLATFORM_MACHINE.startswith(("powerpc", "ppc"))


# --- Utils ------------------------------------------------------------------

def _split_multiline(value):
    value = value.strip()
    sep = max('\n,;', key=value.count)
    return list(filter(None, map(lambda x: x.strip(), value.split(sep))))

def _patch_osx_compiler(compiler):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next((i for i in range(1, len(flags)) if flags[i-1] == "-arch" and flags[i] != platform.machine()), None)
        if i is not None:
            flags.pop(i)
            flags.pop(i-1)

# --- `setuptools` classes ---------------------------------------------------

class PlatformCode:

    def __init__(self, platform, sources, extra_compile_args=None):
        self.platform = platform
        self.sources = sources
        self.extra_compile_args = extra_compile_args or []


class Library(_Library):

    def __init__(self, name, sources, *args, **kwargs):
        self.platform_code = kwargs.pop("platform_code", None) or []
        super().__init__(name, sources, *args, **kwargs)
        for platform_code in self.platform_code:
            self.depends.extend(platform_code.sources)


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


class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
    """

    # --- Silent invocation of the compiler ---

    def _silent_spawn(self, cmd):
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            raise CompileError(err.stderr)

    def _silent_compile(self, *args, **kwargs):
        with mock.patch.object(self.compiler, "spawn", new=self._silent_spawn):
            return self.compiler.compile(*args, **kwargs)

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Build code ---

    def _check_function(self, funcname, header, args="()"):
        print('checking whether function', repr(funcname), 'is available', end="... ", file=sys.stderr)
        self.mkpath(self.build_temp)

        base = "have_{}".format(funcname)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}{};
                    return 0;
                }}
            """.format(header, funcname, args))
        try:
            objects = self.compiler.compile([testfile], debug=self.debug)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
        except CompileError:
            print("no", file=sys.stderr)
            return False
        else:
            print("yes", file=sys.stderr)
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _has_dlopen(self):
        self.mkpath(self.build_temp)
        filename = os.path.join(self.build_temp, "test_dlopen.c")
        with open(filename, "w") as f:
            f.write("""
                #include <dlfcn.h>
                int main() { dlopen("", 0); }
            """)
        try:
            self.compiler.compile([filename], debug=True)
        except CompileError:
            log.warn("could not find `dlopen` function from <dlfcn.h>")
            return False
        else:
            log.info("found `dlopen` function from <dlfcn.h>")
            return True

    def _has_getauxval(self):
        self.mkpath(self.build_temp)
        filename = os.path.join(self.build_temp, "test_dlopen.c")
        with open(filename, "w") as f:
            f.write("""
                #include <sys/getauxval.h>
                int main() { getauxval(0); }
            """)
        try:
            self._silent_compile([filename], debug=True)
        except CompileError:
            log.warn("could not find `getauxval` function from <sys/getauxval.h>")
            return False
        else:
            log.info("found `getauxval` function from <sys/getauxval.h>")
            return True

    def _add_platform_defines(self, library):
        # add some compatibility defines from `cpu_features`
        if PROCESSOR_IS_X86 and sys.platform == "darwin":
            library.define_macros.append(("HAVE_SYSCTLBYNAME", 1))
        if UNIX:
            hardware_detect = False
            if self._has_dlopen():
                library.define_macros.append(("HAVE_DLFCN_H", 1))
                hardware_detect = True
            if self._has_getauxval():
                library.define_macros.append(("HAVE_STRONG_GETAUXVAL", 1))
                hardware_detect = True
            if hardware_detect:
                library.sources.append(os.path.join("vendor", "cpu_features", "src", "hwcaps.c"))

    def build_libraries(self, libraries):
        # check for functions required for libcpu_features on OSX
        if PLATFORM_SYSTEM == "Darwin":
            _patch_osx_compiler(self.compiler)
            if self._check_function("sysctlbyname", "sys/sysctl.h", args="(NULL, NULL, 0, NULL, 0)"):
                self.compiler.define_macro("HAVE_SYSCTLBYNAME", 1)

        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            sources = library.sources.copy()
            for platform_code in library.platform_code:
                sources.extend(platform_code.sources)
            self.make_file(
                sources,
                self.compiler.library_filename(library.name, output_dir=self.build_clib),
                self.build_library,
                (library,)
            )

    def build_library(self, library):
        # add specific code for `cpu_features`
        if library.name == "cpu_features":
            self._add_platform_defines(library)

        # update compile flags if compiling in debug or release mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-Og")
                library.extra_compile_args.append("--coverage")
                library.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Od")

        # attempt to build the platform specific code
        extra_objects = []
        for platform_code in library.platform_code:
            extra_preargs = library.extra_compile_args + platform_code.extra_compile_args
            try:
                extra_objects.extend(self._silent_compile(
                    platform_code.sources,
                    output_dir=self.build_temp,
                    include_dirs=library.include_dirs + [self.build_clib],
                    macros=library.define_macros,
                    debug=self.debug,
                    depends=library.depends,
                    extra_preargs=extra_preargs,
                ))
            except CompileError:
                log.warn(f"failed to compile platform-specific {platform_code.platform} code")
            else:
                log.info(f"successfully built platform-specific {platform_code.platform} code")
                self.compiler.define_macro(f"{platform_code.platform}_BUILD_SUPPORTED")

        # build objects and create a static library
        objects = self.compiler.compile(
            library.sources,
            output_dir=self.build_temp,
            include_dirs=library.include_dirs + [self.build_clib],
            macros=library.define_macros,
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )
        self.compiler.create_static_lib(
            objects + extra_objects,
            library.name,
            output_dir=self.build_clib,
            debug=self.debug,
        )


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    # --- Autotools-like helpers ---

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._pub_cmd = self.get_finalized_command("publicize_headers")
        self._clib_cmd.force = self._pub_cmd.force = self.force
        self._clib_cmd.debug = self.debug
        self.build_clib = self._clib_cmd.build_clib

    def run(self):
        # make sure sources have been patched to expose private fields
        if not self.distribution.have_run.get("publicize_headers", False):
            self._pub_cmd.run()
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include", "pyfastani"],
            "compiler_directives": {},
            "compile_time_env": {
                "FASTANI_PRIVATE_ACCESS": 1,
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
            }
        }
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

    def build_extensions(self):
        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)
        # build the extensions as normal
        _build_ext.build_extensions(self)

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

        # C++ OS-specific options
        if ext.language == "c++":
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

        # add extra static objects
        for lib in ext.libraries:
            ext.extra_objects.append(self.compiler.library_filename(
                lib,
                output_dir=self.build_clib
            ))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class clean(_clean):

    def remove_file(self, filename):
        if os.path.exists(filename):
            log.info("removing {!r}".format(filename))
            os.remove(filename)
        else:
            log.info("{!r} does not exist -- can't clean it".format(filename))

    def run(self):
        _clean.run(self)

        _build_cmd = self.get_finalized_command("build_ext")
        _build_cmd.inplace = True

        for ext in self.distribution.ext_modules:
            filename = _build_cmd.get_ext_filename(ext.name)
            if self.all:
                self.remove_file(filename)
            basename = _build_cmd.get_ext_fullname(ext.name).replace(".", os.path.sep)
            for ext in ["c", "cpp", "html"]:
                filename = os.path.extsep.join([basename, ext])
                self.remove_file(filename)


# --- Cython extensions ------------------------------------------------------



libraries = [
    Library(
        "cpu_features",
        include_dirs=[os.path.join("vendor", "cpu_features", "include")],
        define_macros=[("STACK_LINE_READER_BUFFER_SIZE", 1024)],
        sources=[
            os.path.join("vendor", "cpu_features", "src", x)
            for x in ("filesystem.c", "stack_line_reader.c", "string_view.c")
        ],
        platform_code=[
            PlatformCode(
                platform="AARCH64",
                sources=[os.path.join("vendor", "cpu_features", "src", "cpuinfo_aarch64.c")]
            ),
            PlatformCode(
                platform="ARM",
                sources=[os.path.join("vendor", "cpu_features", "src", "cpuinfo_arm.c")]
            ),
            PlatformCode(
                platform="X86",
                sources=[os.path.join("vendor", "cpu_features", "src", "cpuinfo_x86.c")]
            )
        ]
    ),
    Library(
        "sequtils",
        include_dirs=[
            os.path.join("pyfastani", "_sequtils"),
            os.path.join("vendor", "cpu_features", "include")
        ],
        sources=[os.path.join("pyfastani", "_sequtils", "sequtils.cpp")],
        libraries=["cpu_features"],
        language="c++",
        platform_code=[
            PlatformCode(
                platform="NEON",
                sources=[os.path.join("pyfastani", "_sequtils", "neon.c")],
                extra_compile_args=[] if PROCESSOR_IS_AARCH64 else ["-mfpu=neon"],
            ),
            PlatformCode(
                platform="SSE2",
                sources=[os.path.join("pyfastani", "_sequtils", "sse2.c")],
                extra_compile_args=["-msse2"],
            ),
            PlatformCode(
                platform="SSSE3",
                sources=[os.path.join("pyfastani", "_sequtils", "ssse3.c")],
                extra_compile_args=["-mssse3"],
            )
        ]
    ),
]

extensions = [
    Extension(
        "pyfastani._fastani",
        [os.path.join("pyfastani", x) for x in ("_utils.cpp", "omp.cpp", "_fastani.pyx")],
        language="c++",
        include_dirs=[
            "include",
            "pyfastani",
            os.path.join("vendor", "cpu_features", "include"),
            os.path.join("vendor", "boost-math", "include"),
            os.path.join("pyfastani", "_sequtils"),
        ],
        libraries=["sequtils"],
        define_macros=[
            ("USE_BOOST", 1), # used to compile FastANI without GSL
            ("BOOST_NO_EXCEPTIONS", 1), # don't raise overflow errors
            ("BOOST_IF_CONSTEXPR", "if"),
        ],
    ),
    Extension(
        "pyfastani._fasta",
        [os.path.join("pyfastani", "_fasta.pyx")],
        include_dirs=[
            "include",
            "pyfastani",
            os.path.join("vendor", "cpu_features", "include"),
            os.path.join("pyfastani", "simd")
        ],
        language="c",
        libraries=["sequtils"],
    )
]

# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=extensions,
    libraries=libraries,
    cmdclass=dict(
        build_clib=build_clib,
        build_ext=build_ext,
        publicize_headers=publicize_headers,
        clean=clean,
        sdist=sdist,
    ),
)
