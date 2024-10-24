# Contributing to PyFastANI

For bug fixes or new features, please file an issue before submitting a
pull request. If the change isn't trivial, it may be best to wait for
feedback.

## Setting up a local repository

Make sure you clone the repository in recursive mode, so you also get the
wrapped code of FastANI which is exposed as a ``git`` submodule:

```console
$ git clone --recursive https://github.com/althonos/pyfastani
```

## Compiling the extension

Compiling requires the `boost::math` module from [Boost](https://www.boost.org/).
Depending on your system, you may have to install them yourself.

To compile the extension, use the following command:

```console
$ python -m pip install -v -e . --no-build-isolation
```

## Running tests

Tests are written as usual Python unit tests with the `unittest` module of
the standard library. Running them requires the extension to be built
locally:

```console
$ python -m pip install -v -e . --no-build-isolation
$ python -m unittest pyfastani.tests -vv
```

## Running benchmarks

### Query fragments mapping

The `benches` folder contains benchmarks for evaluating the performance of
the node connection scoring step, essentially to make sure that the
multi-threading makes it faster.

Start by building `pyfastani` locally, and making sure it is compiled
in release mode:
```console
$ python -m pip install -v . --no-build-isolation
```

Then make sure you have the required packages and data:
```console
$ pip install --user -r benches/mapping/requirements.txt
$ python benches/data/download.py
```

Finally, run the benchmarks and plot the results:
```console
$ python benches/mapping/bench.py -d benches/data/ -o times.json
$ python benches/mapping/plot.py -i times.json --show
```

## Coding guidelines

This project targets Python 3.7 or later.

Python objects should be typed; since it is not supported by Cython,
you must manually declare types in type stubs (`.pyi` files). In Python
files, you can add type annotations to function signatures (supported in
Python 3.5) or in variable assignments (supported from Python 3.6
onward).

### Interfacing with C

When interfacing with C, and in particular with pointers, use assertions
everywhere you assume the pointer to be non-NULL.

### Interfacing with C++

When wrapping objects, use stack allocation where possible.
