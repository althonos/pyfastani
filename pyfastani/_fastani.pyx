# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp, binding=False
"""Bindings to FastANI, a method for fast whole-genome similarity estimation.

References:
    - Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S.
      *High throughput ANI analysis of 90K prokaryotic genomes reveals clear
      species boundaries*. Nat Commun. 2018 Nov 30;9(1):5114.
      :doi:`10.1038/s41467-018-07641-9`. :PMID:`30504855`. :PMC:`PMC6269478`.

"""

# --- C imports --------------------------------------------------------------

cimport cython
cimport cython.parallel
cimport libcpp11.chrono
from cpython.buffer cimport Py_buffer, PyBUF_READ
from cython.operator cimport dereference, preincrement, postincrement
from cpython.ref cimport Py_INCREF
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.memoryview cimport PyMemoryView_FromMemory, PyMemoryView_Check, PyMemoryView_GET_BUFFER
from libc.stdio cimport printf
from libc.string cimport memcpy
from libc.limits cimport INT_MAX
from libc.stdint cimport int64_t, uint64_t
from libc.stdlib cimport malloc, realloc, free
from libcpp cimport bool, nullptr
from libcpp.algorithm cimport sort, unique
from libcpp.deque cimport deque
from libcpp.utility cimport move, pair
from libcpp.functional cimport function
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp11.fstream cimport ofstream
from openmp cimport omp_lock_t, omp_init_lock, omp_set_lock, omp_unset_lock

cimport fastani.map.map_stats
from kseq cimport kseq_t, kstring_t
from fastani.cgi.compute_core_identity cimport computeCGI
from fastani.cgi.cgid_types cimport CGI_Results
from fastani.map cimport base_types
from fastani.map.compute_map cimport (
    Map as Map_t,
    L1_candidateLocus_t,
    L2_mapLocus_t
)
from fastani.map.common_func cimport addMinimizers, getHash
from fastani.map.map_parameters cimport Parameters as Parameters_t
from fastani.map.map_stats cimport estimateMinimumHitsRelaxed, recommendedWindowSize, j2md, md_lower_bound
from fastani.map.win_sketch cimport Sketch as Sketch_t
from fastani.map.base_types cimport (
    hash_t,
    seqno_t,
    offset_t,
    ContigInfo as ContigInfo_t,
    MappingResult as MappingResult_t,
    MappingResultsVector_t,
    MinimizerInfo as MinimizerInfo_t,
    MinimizerMetaData as MinimizerMetaData_t,
    MinimizerMapKeyType as MinimizerMapKeyType_t,
    MinimizerMapValueType as MinimizerMapValueType_t,
    QueryMetaData as QueryMetaData_t,
)

# HACK: we need kseq_t* as a template argument, which is not supported by
#       Cython at the moment, so we just `typedef kseq_t* kseq_ptr_t` in
#       an external C++ header to make Cython happy
from _utils cimport kseq_ptr_t, toupper, distance
from _unicode cimport *
from _sequtils cimport copy_upper, reverse_complement
from _atomic_vector cimport atomic_vector

# --- Python imports ---------------------------------------------------------

import array
import os
import multiprocessing.pool
import threading
import warnings

# --- Constants --------------------------------------------------------------

DEF _MAX_KMER_SIZE = 2048
MAX_KMER_SIZE = _MAX_KMER_SIZE


# --- Cython helpers ---------------------------------------------------------

cdef ssize_t _read_nucl(
    const int kind,
    const void* data,
    const ssize_t slen,
    const ssize_t i,
    char* fwd,
    char* bwd,
) nogil except -1:
    """Read characters of a nucleotide sequence into `fwd` and `bwd` buffers.

    Reads at most _MAX_KMER_SIZE characters, and returns the number of
    characters read. ``i`` is the position of the first letter to read
    in the string. Returns the number of characters read.

    """
    cdef ssize_t j
    cdef ssize_t length
    cdef char    nuc

    if _MAX_KMER_SIZE <= slen - i:
        length = _MAX_KMER_SIZE
    elif i < slen:
        length = slen - i
    else:
        length = 0

    # if UCS-1, bytes are next to each other, so we can use the SIMD
    # implementations to copy into uppercase
    if kind == PyUnicode_1BYTE_KIND:
        copy_upper(&fwd[_MAX_KMER_SIZE], &(<char*> data)[i], length)
    else:
        for j in range(length):
            fwd[_MAX_KMER_SIZE + j] = toupper(<int> PyUnicode_READ(kind, data, i + j))

    # reverse complement in backward buffer
    reverse_complement(&bwd[_MAX_KMER_SIZE - length], &fwd[_MAX_KMER_SIZE], length)

    return length


cdef int _add_minimizers_nucl(
    vector[MinimizerInfo_t] &minimizer_index,
    const int kind,
    const void* data,
    const ssize_t slen,
    const int kmer_size,
    const int window_size,
    const seqno_t seq_counter,
) nogil except 1:
    """Add the minimizers for a single contig to the sketcher.

    Adapted from the ``skch::commonFunc::addMinimizers`` method in
    ``commonFunc.hpp`` to work without the ``kseq`` I/O.

    """
    cdef deque[pair[MinimizerInfo_t, int64_t]] q
    cdef void*                                 last
    cdef int64_t                               i
    cdef int64_t                               current_window_id
    cdef hash_t                                hash_fwd
    cdef hash_t                                hash_bwd
    cdef hash_t                                current_kmer
    cdef MinimizerInfo_t                       info
    cdef char                                  fwd[_MAX_KMER_SIZE*2]
    cdef char                                  bwd[_MAX_KMER_SIZE*2]

    cdef uint64_t n = 0

    # initial fill of the buffer for the sequence sliding window
    # supporting any unicode sequence in canonical form (including
    # byte buffers containing ASCII characters)
    _read_nucl(kind, data, slen, 0, fwd, bwd)

    # process all windows of width `kmer_size` in the input sequence
    for i in range(slen - kmer_size + 1):
        # if reaching the end of the sliding window, slide buffer
        # left, and read new block in right side
        if i % _MAX_KMER_SIZE == 0:
            memcpy(&fwd[0], &fwd[_MAX_KMER_SIZE], _MAX_KMER_SIZE)
            memcpy(&bwd[_MAX_KMER_SIZE], &bwd[0], _MAX_KMER_SIZE)
            _read_nucl(kind, data, slen, i + _MAX_KMER_SIZE, fwd, bwd)
        # compute forward hash
        hash_fwd = getHash(<char*> &fwd[i % _MAX_KMER_SIZE], kmer_size)
        hash_bwd = getHash(<char*> &bwd[2*_MAX_KMER_SIZE - i % _MAX_KMER_SIZE - kmer_size], kmer_size)
        n += 1
        # only record asymmetric k-mers
        if hash_bwd != hash_fwd:
            # record window size for the minimizer
            current_window_id = i - window_size + 1
            # Take minimum value of kmer and its reverse complement
            current_kmer = min(hash_fwd, hash_bwd)
            # If front minimum is not in the current window, remove it
            while not q.empty() and q.front().second <= i - window_size:
                q.pop_front()
            # hashes less than equal to current_k;er can be discarded
            while not q.empty() and q.back().first.hash >= current_kmer:
                q.pop_back()
            # push current_kmer and position to back of the queue
            info.hash = current_kmer
            info.seqId = seq_counter
            info.wpos = 0
            q.push_back(pair[MinimizerInfo_t, int64_t](info, i))
            # select the minimizer from Q and put into index
            if current_window_id >= 0:
                if minimizer_index.empty() or minimizer_index.back() != q.front().first:
                    q.front().first.wpos = current_window_id
                    minimizer_index.push_back(q.front().first)


cdef ssize_t _read_prot(
    const int kind,
    const void* data,
    const ssize_t slen,
    const ssize_t i,
    char* fwd,
) nogil except -1:
    cdef ssize_t j
    cdef ssize_t length
    cdef char    nuc

    if _MAX_KMER_SIZE <= slen - i:
        length = _MAX_KMER_SIZE
    elif i < slen:
        length = slen - i
    else:
        length = 0

    # if UCS-1, bytes are next to each other, so we can use the SIMD
    # implementations to copy into uppercase
    if kind == PyUnicode_1BYTE_KIND:
        copy_upper(&fwd[_MAX_KMER_SIZE], &(<char*> data)[i], length)
    else:
        for j in range(length):
            fwd[_MAX_KMER_SIZE + j] = toupper(<int> PyUnicode_READ(kind, data, i + j))


cdef int _add_minimizers_prot(
    vector[MinimizerInfo_t] &minimizer_index,
    const int kind,
    const void* data,
    const ssize_t slen,
    const int kmer_size,
    const int window_size,
    const seqno_t seq_counter,
) nogil except 1:
    """Add the minimizers for a single protein to the sketcher.

    Adapted from the ``skch::commonFunc::addMinimizers`` method in
    ``commonFunc.hpp`` to work without the ``kseq`` I/O.

    """
    cdef deque[pair[MinimizerInfo_t, int64_t]] q
    cdef void*                                 last
    cdef int64_t                               i
    cdef int64_t                               current_window_id
    cdef hash_t                                current_kmer
    cdef MinimizerInfo_t                       info
    cdef char                                  fwd[_MAX_KMER_SIZE*2]

    cdef uint64_t n = 0

    # initial fill of the buffer for the sequence sliding window
    # supporting any unicode sequence in canonical form (including
    # byte buffers containing ASCII characters)
    _read_prot(kind, data, slen, 0, fwd)

    # process all windows of width `kmer_size` in the input sequence
    for i in range(slen - kmer_size + 1):
        # if reaching the end of the sliding window, slide buffer
        # left, and read new block in right side
        if i % _MAX_KMER_SIZE == 0:
            memcpy(&fwd[0], &fwd[_MAX_KMER_SIZE], _MAX_KMER_SIZE)
            _read_prot(kind, data, slen, i + _MAX_KMER_SIZE, fwd)
        # compute forward hash
        current_kmer = getHash(<char*> &fwd[i % _MAX_KMER_SIZE], kmer_size)
        n += 1
        # record window size for the minimizer
        current_window_id = i - window_size + 1
        # If front minimum is not in the current window, remove it
        while not q.empty() and q.front().second <= i - window_size:
            q.pop_front()
        # hashes less than equal to current_k;er can be discarded
        while not q.empty() and q.back().first.hash >= current_kmer:
            q.pop_back()
        # push current_kmer and position to back of the queue
        info.hash = current_kmer
        info.seqId = seq_counter
        info.wpos = 0
        q.push_back(pair[MinimizerInfo_t, int64_t](info, i))
        # select the minimizer from Q and put into index
        if current_window_id >= 0:
            if minimizer_index.empty() or minimizer_index.back() != q.front().first:
                q.front().first.wpos = current_window_id
                minimizer_index.push_back(q.front().first)


# --- Cython classes ---------------------------------------------------------


class _DummyPool:
    """A dummy `~ThreadPool` that runs everything in the main thread.
    """

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def map(self, func, iterable):
        for x in iterable:
            func(x)


cdef class _Map:
    """A private class wrapping a heap-allocated compute map.
    """

    # --- Attributes ---------------------------------------------------------

    cdef Map_t* _map

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._map = NULL

    def __init__(self):
        raise TypeError(f"Cannot instantiate objects of type {type(self).__name__!r}")

    def __dealloc__(self):
        del self._map


cdef class _FinalMappings:
    """A private class wrapping a vector of L2 mapping results.
    """

    # --- Attributes ---------------------------------------------------------

    cdef atomic_vector[MappingResult_t] _vec

    # --- Magic methods ------------------------------------------------------

    def __init__(self):
        raise TypeError(f"Cannot instantiate objects of type {type(self).__name__!r}")


cdef class _Parameterized:
    """A base class for types wrapping a `skch::Parameters` C++ object.
    """

    # --- Attributes ---------------------------------------------------------

    cdef Parameters_t _param

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._param.threads = 1
        self._param.reportAll = True
        self._param.visualize = False
        self._param.matrixOutput = False

    def __getstate__(self):
        return {
            "kmerSize": self._param.kmerSize,
            "windowSize": self._param.windowSize,
            "minReadLength": self._param.minReadLength,
            "minFraction": self._param.minFraction,
            "threads": self._param.threads,
            "alphabetSize": self._param.alphabetSize,
            "referenceSize": self._param.referenceSize,
            "percentageIdentity": self._param.percentageIdentity,
            "p_value": self._param.p_value,
        }

    def __setstate__(self, state):
        self._param.kmerSize = state["kmerSize"]
        self._param.windowSize = state["windowSize"]
        self._param.minReadLength = state["minReadLength"]
        self._param.minFraction = state["minFraction"]
        self._param.threads = state["threads"]
        self._param.alphabetSize = state["alphabetSize"]
        self._param.referenceSize = state["referenceSize"]
        self._param.percentageIdentity = state["percentageIdentity"]
        self._param.p_value = state["p_value"]

    # --- Properties ---------------------------------------------------------

    @property
    def k(self):
        """`int`: The k-mer size used for sketching.
        """
        return self._param.kmerSize

    @property
    def window_size(self):
        """`int`: The window size used for sketching.
        """
        return self._param.windowSize

    @property
    def fragment_length(self):
        """`int`: The minimum read length to use for mapping.
        """
        return self._param.minReadLength

    @property
    def minimum_fraction(self):
        """`float`: The minimum genome fraction required to trust ANI values.
        """
        return self._param.minFraction

    @property
    def percentage_identity(self):
        """`float`: The identity threshold for similarity when estimating hits.
        """
        return self._param.percentageIdentity

    @property
    def p_value(self):
        """`float`: The p-value threshold for similarity when estimating hits.
        """
        return self._param.p_value

    @property
    def protein(self):
        """`bool`: Whether or not the object expects peptides or nucleotides.
        """
        return self._param.alphabetSize == 20


@cython.final
cdef class Sketch(_Parameterized):
    """An index computing minimizers over the reference genomes.

    Use this class to add reference genomes with the `add_genome` or
    `add_draft` methods, then call the `index` method to obtain a `Mapper`
    that can be used to map query genomes.

    Attributes:
        minimizers (`~pyfastani.Minimizers`): A view over the minimizers
            currently recorded in the sketch.

    """

    # --- Attributes ---------------------------------------------------------

    cdef          Sketch_t*        _sk          # the internal Sketch_t
    cdef          size_t           _counter     # the number of contigs (not genomes)
    cdef          vector[uint64_t] _lengths     # array mapping each genome to its length
    cdef          list             _names       # list mapping each genome to its name
    cdef readonly Minimizers        minimizers  # a view over the minimizers

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        # create a new Sketch with the parameters
        self._sk = new Sketch_t(self._param)
        # create the minimizers view
        self.minimizers = Minimizers()
        self.minimizers._owner = self
        self.minimizers._minimizers = &self._sk.minimizerIndex
        # create a new list of names
        self._names = []

    def __init__(
        self,
        *,
        unsigned int k=16,
        unsigned int fragment_length=3000,
        float minimum_fraction=0.2,
        double p_value=1e-03,
        float percentage_identity=80.0,
        uint64_t reference_size=5_000_000,
        bint protein=False,
    ):
        """__init__(self, *, k=16, fragment_length=3000, minimum_fraction=0.2, p_value=1e-03, percentage_identity=80, reference_size=5e9, protein=False)\n--

        Create a new FastANI sequence sketch.

        Keyword Arguments:
            k (`int`): The size of the k-mers. FastANI authors recommend
                a size of at most 16, but any positive number below up to
                `pyfastani.MAX_KMER_SIZE` will work.
            fragment_length (`int`): The lengths the blocks should have
                when splitting the query. Queries smaller than this number
                won't be processed.
            minimum_fraction (`float`): The minimum fraction of genome that
                must be shared for a hit to be reported. If reference and
                query genome size differ, smaller one among the two is
                considered.
            p_value (`float`): The p-value cutoff. *Used to determine the
                recommended window size.*
            percentage_identity (`int`): An identity percentage above which
                ANI values between two sequences can be trusted. *Used to
                to determine the recommended window size.*
            reference_size (`int`): An estimate of the reference length.
                *Used to determine the recommended window size.*
            protein (`bool`): Whether or not protein sequences are expected.
                If `True`, the alphabet size is changed from 4 to 20,
                minimizers are not computed on the "reverse" strand, and the
                window size is set to 1.

        """
        if minimum_fraction > 1 or minimum_fraction < 0:
            raise ValueError(f"minimum_fraction must be between 0 and 1, got {minimum_fraction!r}")
        if fragment_length <= 0:
            raise ValueError(f"fragment_length must be strictly positive, got {fragment_length!r}")
        if p_value <= 0:
            raise ValueError(f"p_value must be positive, got {p_value!r}")
        if percentage_identity > 100 or percentage_identity < 0:
            raise ValueError(f"percentage_identity must be between 0 and 100, got {percentage_identity!r}")
        if k <= 0:
            raise ValueError(f"k must be strictly positive, got {k!r}")
        elif k > _MAX_KMER_SIZE:
            raise BufferError(f"k must be smaller than {_MAX_KMER_SIZE}, got {k}")
        elif k > 16:
            warnings.warn(
                f"Using k-mer size greater than 16 ({k!r}), accuracy will be degraded.",
                UserWarning,
            )

        # store parameters
        self._param.kmerSize = k
        self._param.minReadLength = fragment_length
        self._param.minFraction = minimum_fraction
        self._param.p_value = p_value
        self._param.percentageIdentity = percentage_identity
        self._param.referenceSize = reference_size

        if protein:
            self._param.alphabetSize = 20
            self._param.windowSize = 1
        else:
            self._param.alphabetSize = 4
            self._param.windowSize = fastani.map.map_stats.recommendedWindowSize(
                self._param.p_value,
                self._param.kmerSize,
                self._param.alphabetSize,
                self._param.percentageIdentity,
                self._param.minReadLength,
                self._param.referenceSize
            )

        # initialize bookkeeping values and make sure self._sk is cleared
        # (in case __init__ is called more than once)
        self.clear()

    def __dealloc__(self):
        del self._sk

    def __getstate__(self):
        return {
            "parameters": _Parameterized.__getstate__(self),
            "counter": self._counter,
            "lengths": list(self._lengths),
            "names": list(self._names),
            "sketch": {
                "sequencesByFileInfo": list(self._sk.sequencesByFileInfo),
                "minimizers": self.minimizers.__getstate__(),
            }
        }

    def __setstate__(self, state):
        _Parameterized.__setstate__(self, state["parameters"])
        self._counter = state["counter"]
        self._lengths = state["lengths"]
        self._names = state["names"]

        self._sk.sequencesByFileInfo = state["sketch"]["sequencesByFileInfo"]
        self.minimizers.__setstate__(state["sketch"]["minimizers"])

    # --- Properties ---------------------------------------------------------

    @property
    def occurences_threshold(self):
        """`int`: The occurence threshold above which minimizers are ignored.
        """
        assert self._sk != nullptr
        return self._sk.getFreqThreshold()

    @property
    def names(self):
        """`list` of `str`: The names of the sequences currently sketched.
        """
        return self._names[:]

    # --- Methods ------------------------------------------------------------

    cdef int _add_draft(self, object name, object contigs) except 1:
        """Add a draft genome to the sketcher.

        Adapted from the ``skch::Sketch::build`` method in ``winSketch.hpp``
        to work without the ``kseq`` I/O.

        """
        assert self._sk != nullptr

        cdef size_t                   total  = 0
        cdef Parameters_t*            param  = &self._param
        cdef object                   contig

        # variables to index the contig as text or bytes
        cdef const unsigned char[::1] view
        cdef int                      kind
        cdef void*                    data
        cdef ssize_t                  slen

        for contig in contigs:

            # get a way to read each letter of the contig,
            # independently of it being `str`, `bytes`, `bytearray`, etc.
            if isinstance(contig, str):
                # make sure the unicode string is in canonical form,
                # --> won't be needed anymore in Python 3.12
                IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
                    PyUnicode_READY(contig)
                # get kind and data for efficient indexing
                kind = PyUnicode_KIND(contig)
                data = PyUnicode_DATA(contig)
                slen = PyUnicode_GET_LENGTH(contig)
            else:
                # attempt to view the contig as a buffer of contiguous bytes
                view = contig
                # pretend the the bytes are an ASCII (UCS-1) encoded string
                kind = PyUnicode_1BYTE_KIND
                slen = view.shape[0]
                if slen != 0:
                    data = <void*> &view[0]

            # check the sequence is large enough to compute minimizers
            if slen >= param.windowSize and slen >= param.kmerSize:
                with nogil:
                    if param.alphabetSize == 4:
                        _add_minimizers_nucl(
                            self._sk.minimizerIndex,
                            kind,
                            data,
                            slen,
                            param.kmerSize,
                            param.windowSize,
                            self._counter,
                        )
                    else:
                        _add_minimizers_prot(
                            self._sk.minimizerIndex,
                            kind,
                            data,
                            slen,
                            param.kmerSize,
                            param.windowSize,
                            self._counter,
                        )
            else:
                warnings.warn(
                  (
                    "Sketch received a short contig relative to parameters, "
                    "minimizers will not be added."
                  ),
                  UserWarning
                )

            # compute the genome length with respect to the fragment length
            total += (slen // param.minReadLength) * param.minReadLength

            # record the number of contigs currently stored
            self._counter += 1

        # record the reference genome name and total length
        self._names.append(name)
        self._lengths.push_back(total)

        # record the number of contigs in this genome
        self._sk.sequencesByFileInfo.push_back(self._counter)

    cpdef Sketch add_draft(self, object name, object contigs):
        """add_draft(self, name, contigs)\n--

        Add a reference draft genome to the sketcher.

        Using this method is fine even when the genome has a single contig,
        although `Sketch.add_genome` is easier to use in that case.

        Arguments:
            name (`object`): The name of the genome to add. When a reference
                matches this query genome, ``name`` will be exposed as the
                `Hit.name` attribute of the corresponding hit.
            contigs (iterable of `str` or `bytes`): The contigs of the genome.

        Returns:
            `Sketch`: the object itself, for method chaining.

        Hint:
            Contigs smaller than the window size and the k-mer size will
            be skipped.

        """
        # delegate to the C code
        self._add_draft(name, contigs)
        return self

    cpdef Sketch add_genome(self, object name, object sequence):
        """add_genome(self, name, sequence)\n--

        Add a reference genome to the sketcher.

        This method is a shortcut for `Sketch.add_draft` when a genome is
        complete (i.e. only contains a single contig).

        Arguments:
            name (`object`): The name of the genome to add. When a reference
                matches this query genome, ``name`` will be exposed as the
                `Hit.name` attribute of the corresponding hit.
            sequence (`str` or `bytes`): The sequence of the genome.

        Returns:
            `Sketch`: the object itself, for method chaining.

        Hint:
            Sequence must be larger than the window size and the k-mer size
            to be sketched, otherwise no minifiers will be computed.

        """
        # delegate to the C code
        self._add_draft(name, (sequence,))
        return self

    cpdef Sketch clear(self):
        """clear(self)\n--

        Reset the `Sketch`, removing any reference genome it may contain.

        Returns:
            `Sketch`: the object itself, for method chaining.

        """
        # reset self
        self._names.clear()
        self._lengths.clear()
        self._counter = 0
        # reset self._sk
        self._sk.freqThreshold = INT_MAX
        self._sk.metadata.clear()
        self._sk.sequencesByFileInfo.clear()
        self._sk.minimizerPosLookupIndex.clear()
        self._sk.minimizerIndex.clear()
        self._sk.minimizerFreqHistogram.clear()
        # return self for chaining
        return self

    cpdef Mapper index(self):
        """index(self)\n--

        Index the reference genomes for fast lookups using the minimizers.

        Once all the reference sequences have been added to the `Sketch`,
        use this method to create an efficient mapper, dropping the most
        common minifiers among the reference sequences.

        Returns:
            `~pyfastani.Mapper`: An indexed mapper that can be used
            for fast querying.

        Note:
            Calling this method will effectively transfer ownership of
            the data to the `Mapper`, and reset the internals of this
            `Sketch`. It will be essentially cleared, but should remain
            usable.

        """
        # compute the index and the frequency histogram
        self._sk.index()
        self._sk.computeFreqHist()
        # create the Mapper for this Sketch
        cdef Mapper mapper = Mapper.__new__(Mapper)
        mapper._param = self._param # copy params
        mapper._sk = self._sk
        mapper._names = self._names.copy()
        mapper._lengths.swap(self._lengths)
        # set the minimizers
        mapper.minimizers = Minimizers.__new__(Minimizers)
        mapper.minimizers._owner = mapper
        mapper.minimizers._minimizers = &mapper._sk.minimizerIndex
        # reset the current sketch
        self._sk = new Sketch_t(self._param)
        self.clear()
        # return the new mapper
        return mapper


@cython.final
cdef class Mapper(_Parameterized):
    """A genome mapper using Murmur3 hashes and k-mers to compute ANI.

    Attributes:
        minimizers (`~pyfastani.Minimizers`): A view over the minimizers
            recorded in the mapper.

    """

    # --- Attributes ---------------------------------------------------------

    cdef          Sketch_t*        _sk
    cdef          vector[uint64_t] _lengths
    cdef          list             _names
    cdef readonly Minimizers       minimizers

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        # create a new Sketch with the parameters
        self._sk = new Sketch_t(self._param)
        # create the minimizers view
        self.minimizers = Minimizers()
        self.minimizers._owner = self
        self.minimizers._minimizers = &self._sk.minimizerIndex

    def __init__(self, *args, **kwargs):
        raise TypeError("Mapper cannot be instantiated, use `Sketch.index` instead.")

    def __dealloc__(self):
        del self._sk

    def __getstate__(self):
        return {
            "parameters": _Parameterized.__getstate__(self),
            "lengths": list(self._lengths),
            "names": list(self._names),
            "sketch": {
                "sequencesByFileInfo": list(self._sk.sequencesByFileInfo),
                "minimizers": self.minimizers.__getstate__(),
                "minimizerFreqHistogram": dict(self._sk.minimizerFreqHistogram),
                "minimizerPosLookupIndex": self.lookup_index,
            }
        }

    def __setstate__(self, state):
        _Parameterized.__setstate__(self, state["parameters"])
        self._lengths = state["lengths"]
        self._names = state["names"]

        self._sk.minimizerFreqHistogram = state["sketch"]["minimizerFreqHistogram"]
        self._sk.sequencesByFileInfo = state["sketch"]["sequencesByFileInfo"]
        self.minimizers.__setstate__(state["sketch"]["minimizers"])

        cdef Position pos
        cdef vector[MinimizerMetaData_t] map_value
        self._sk.minimizerPosLookupIndex = unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t]()
        for key, value in state["sketch"]["minimizerPosLookupIndex"].items():
            map_value = vector[MinimizerMetaData_t]()
            for pos in value:
                map_value.push_back(pos.to_raw())
            self._sk.minimizerPosLookupIndex.insert(pair[MinimizerMapKeyType_t, MinimizerMapValueType_t](key, map_value))

    # --- Properties ---------------------------------------------------------

    @property
    def lookup_index(self):
        """`MinimizerLookupIndex`: The index of initial minimizer positions.

        This table is used to retrieve at which positions the minimizers
        appear in the reference genomes.

        """
        assert self._sk != nullptr
        cdef MinimizerIndex index = MinimizerIndex.__new__(MinimizerIndex)
        index._map = &self._sk.minimizerPosLookupIndex
        index.owner = self
        return index

    # --- Methods ------------------------------------------------------------

    @staticmethod
    cdef void _do_l1_mappings(
        Map_t& map,
        const int kind,
        const void* data,
        const ssize_t slen,
        QueryMetaData_t[kseq_ptr_t, vector[MinimizerInfo_t]]& query,
        vector[L1_candidateLocus_t]& l1_mappings,
    ) nogil:
        """Compute L1 mappings for the given sequence block.

        Adapted from the `skch::Map::doL1Mapping` in `computeMap.hpp` to avoid
        reading from a `kseq_t`, and support indexing a Python Unicode string.

        """
        cdef vector[MinimizerMetaData_t] seed_hits_l1
        cdef vector[MinimizerInfo_t].iterator uniq_end_iter
        cdef vector[MinimizerInfo_t].iterator it
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].const_iterator seed_find
        cdef MinimizerMapValueType_t hit_position_list

        # compute minimizers
        if map.param.alphabetSize == 4:
            _add_minimizers_nucl(
               query.minimizerTableQuery,
               kind,
               data,
               slen,
               map.param.kmerSize,
               map.param.windowSize,
               0,
            )
        else:
            _add_minimizers_prot(
               query.minimizerTableQuery,
               kind,
               data,
               slen,
               map.param.kmerSize,
               map.param.windowSize,
               0,
            )

        # find the unique minimizers in thos that were just obtained
        sort(query.minimizerTableQuery.begin(), query.minimizerTableQuery.end(), &MinimizerInfo_t.lessByHash)

        # manually implement `unique` as template instantiation has issues on OSX
        it = query.minimizerTableQuery.begin()
        uniq_end_iter = unique(query.minimizerTableQuery.begin(), query.minimizerTableQuery.end(), &MinimizerInfo_t.equalityByHash)

        # early return if no minimizers were found
        query.sketchSize = distance(query.minimizerTableQuery.begin(), uniq_end_iter)
        if query.sketchSize == 0:
            return

        # keep minimizer if it exist in the reference lookup index
        it = query.minimizerTableQuery.begin()
        while it != uniq_end_iter:
            seed_find = map.refSketch.minimizerPosLookupIndex.const_find(dereference(it).hash)
            if seed_find != map.refSketch.minimizerPosLookupIndex.end():
                hit_position_list = dereference(seed_find).second
                if hit_position_list.size() < map.refSketch.getFreqThreshold():
                    seed_hits_l1.insert(seed_hits_l1.end(), hit_position_list.begin(), hit_position_list.end())
            postincrement(it)

        # estimate the number of minimum hits, and compute candidates
        minimum_hits = estimateMinimumHitsRelaxed(query.sketchSize, map.param.kmerSize, map.param.percentageIdentity)
        map.computeL1CandidateRegions(query, seed_hits_l1, minimum_hits, l1_mappings)

    cpdef void _query_fragment(
        self,
        _Map map,
        const int i,
        const seqno_t seq_counter,
        object mem,
        int kind,
        int stride,
        _FinalMappings final_mappings,
    ) except *:
        cdef kseq_t                                         kseq
        cdef QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type] query
        cdef vector[L1_candidateLocus_t]                    l1_mappings
        cdef Py_buffer*                                     buffer
        cdef char*                                          data
        cdef char*                                          fragment

        if __debug__:
            if not PyMemoryView_Check(mem):
                raise TypeError(f"expected memoryview, got {type(mem).__name__!r}")
        buffer = PyMemoryView_GET_BUFFER(mem)

        with nogil:
            data = <char*> buffer.buf
            fragment = &data[i*self._param.minReadLength*stride]

            query.kseq = &kseq
            query.kseq.seq.s = NULL
            query.kseq.seq.l = self._param.minReadLength
            query.seqCounter = seq_counter + i

            Mapper._do_l1_mappings(
                map._map[0],
                # start of the i-th fragment
                kind,
                fragment,
                self._param.minReadLength,
                # outputs
                query,
                l1_mappings,
            )

            map._map.doL2Mapping(
                query,
                l1_mappings,
                final_mappings._vec
            )

    cdef list _query_draft(self, object contigs, int threads=0):
        """Query the sketcher for the given contigs.

        Adapted from the ``skch::Map::mapQuery`` method in ``computeMap.hpp``.

        """
        assert self._sk != nullptr

        # iterators over the contigs
        cdef int                               i               # fragment counter
        cdef object                            contig
        cdef int                               fragment_count
        # bookkeeping
        cdef uint64_t                          total_fragments = 0
        cdef uint64_t                          total_length    = 0
        # compute map and L2 mappings
        cdef _Map                              map
        cdef _FinalMappings                    final_mappings
        # sequence as a unicode object
        cdef const unsigned char[::1]          view
        cdef int                               kind
        cdef int                               stride
        cdef void*                             data
        cdef ssize_t                           slen
        cdef object                            mem
        # parallel computation
        cdef object                            pool
        # core genomic identity results
        cdef vector[CGI_Results]               results
        cdef CGI_Results                       result
        # filtering and reporing results
        cdef uint64_t                          min_length
        cdef uint64_t                          shared_length
        cdef list                              hits            = []

        # run everything in main thread if `threads==1`, otherwise run a
        # `ThreadPool` with the required number of threads, if specified.
        if threads == 0:
            threads = os.cpu_count() or 1
        if threads == 1:
            pool = _DummyPool()
        elif threads > 1:
            pool = multiprocessing.pool.ThreadPool(threads)
        else:
            raise ValueError(f"`threads` must be positive or null, got {threads!r}")

        # create a new mapper and a new l2 mapping vector
        map = _Map.__new__(_Map)
        map._map = new Map_t(self._param, self._sk[0], total_fragments, 0)
        final_mappings = _FinalMappings.__new__(_FinalMappings)

        # spawn a thread pool to map fragments in parallel for all the contigs
        with pool:
            for contig in contigs:
                # check length of contig is enough for computing mapping
                slen = len(contig)
                if slen < min(self._param.windowSize, self._param.kmerSize, self._param.minReadLength):
                    warnings.warn(
                        (
                            "Mapper received a short sequence relative to parameters, "
                            "mapping will not be computed."
                        ),
                        UserWarning,
                    )
                # get a way to read each letter of the contig,
                # independently of it being `str`, `bytes`, `bytearray`, etc.
                if isinstance(contig, str):
                    # make sure the unicode string is in canonical form,
                    # --> won't be needed anymore in Python 3.12
                    IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
                        PyUnicode_READY(contig)
                    # get kind and data for efficient indexing
                    kind = PyUnicode_KIND(contig)
                    data = PyUnicode_DATA(contig)
                    slen = PyUnicode_GET_LENGTH(contig)
                    if kind == PyUnicode_1BYTE_KIND:
                        stride = sizeof(Py_UCS1)
                    elif kind == PyUnicode_2BYTE_KIND:
                        stride = sizeof(Py_UCS2)
                    else:
                        stride = sizeof(Py_UCS4)
                else:
                    # attempt to view the contig as a buffer of contiguous bytes
                    view = contig
                    # pretend the bytes are an ASCII (UCS-1) encoded string
                    kind = PyUnicode_1BYTE_KIND
                    slen = view.shape[0]
                    stride = sizeof(Py_UCS1)
                    if slen != 0:
                        data = <void*> &view[0]
                # turn the void pointer into a Python object so it can
                # be passed to another Python method
                mem = PyMemoryView_FromMemory(<char*> data, slen*stride, PyBUF_READ)
                # compute the expected number of blocks
                fragment_count = slen // self._param.minReadLength
                # map the blocks in parallel
                pool.map(
                    lambda i: self._query_fragment(map, i, total_fragments, mem, kind, stride, final_mappings),
                    range(fragment_count),
                )
                # record the number of fragments
                total_fragments += fragment_count
                total_length += slen

        # compute core genomic identity after successful mapping
        with nogil:
            computeCGI(
                self._param,
                final_mappings._vec,
                map._map[0],
                self._sk[0],
                total_fragments, # total query fragments
                0, # queryFileNo, only used for visualization, ignored
                string(), # fileName, only used for reporting, ignored
                results,
            )

        # build and return the list of hits
        for result in results:
            assert result.refGenomeId < self._lengths.size()
            assert result.refGenomeId < len(self._names)
            min_length = min(total_length, self._lengths[result.refGenomeId])
            shared_length = result.countSeq * self._param.minReadLength
            if shared_length >= min_length * self._param.minFraction:
                hits.append(Hit(
                    name=self._names[result.refGenomeId],
                    identity=result.identity,
                    matches=result.countSeq,
                    fragments=result.totalQueryFragments,
                ))
        return hits

    cpdef list query_draft(self, object contigs, int threads=0):
        """query_draft(self, contigs, threads=0)\n--

        Query the mapper for a complete genome.

        Arguments:
            contigs (iterable or `str` or `bytes`): The genome to query the
                mapper with.
            threads (`int`): The number of threads to use to run the
                fragment mapping in parallel. Pass *0* (the default) to
                auto-detect the number of threads on the local machine.

        Returns:
            `list` of `~pyfastani.Hit`: The hits found for the query.

        Hint:
            Sequence must be larger than the window size, the k-mer size,
            and the fragment length to be mapped, otherwise an empty list
            of hits will be returned.

        Note:
            This method is reentrant and releases the GIL when hashing
            the blocks allowing to query the mapper in parallel for
            several individual genomes.

        .. versionadded:: 0.4.0
           The `threads` argument.

        """
        # delegate to C code
        return self._query_draft(contigs, threads=threads)

    cpdef list query_genome(self, object sequence, int threads=0):
        """query_genome(self, sequence, threads=0)\n--

        Query the mapper for a complete genome.

        Arguments:
            sequence (`str` or `bytes`): The closed genome to query the
                mapper with.
            threads (`int`): The number of threads to use to run the
                fragment mapping in parallel. Pass *0* (the default) to
                auto-detect the number of threads on the local machine.

        Returns:
            `list` of `~pyfastani.Hit`: The hits found for the query.

        Hint:
            Sequence must be larger than the window size, the k-mer size,
            and the fragment length to be mapped, otherwise an empty list
            of hits will be returned.

        Note:
            This method is reentrant and releases the GIL when hashing
            the blocks allowing to query the mapper in parallel for
            several individual genomes.

        .. versionadded:: 0.4.0
           The `threads` argument.

        """
        # delegate to C code
        return self._query_draft((sequence,), threads=threads)


cdef class Minimizers:
    """A read-only view over the minimizers of a `Sketch` or a `Mapper`.
    """

    # --- Attributes ---------------------------------------------------------

    cdef object                   _owner
    cdef vector[MinimizerInfo_t]* _minimizers

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._owner = None
        self._minimizers = NULL

    def __dealloc__(self):
        if self._owner is None:
            del self._minimizers

    def __len__(self):
        return 0 if self._minimizers == NULL else self._minimizers[0].size()

    def __getitem__(self, ssize_t index):
        if self._minimizers == NULL:
            raise IndexError(index)
        cdef ssize_t index_ = index
        cdef ssize_t length = self._minimizers[0].size()
        if index_ < 0:
            index_ += length
        if index_ < 0 or index >= length:
            raise IndexError(index)
        assert self._minimizers != NULL
        return MinimizerInfo.from_raw(self._minimizers[0][index_])

    cpdef dict __getstate__(self):
        cdef MinimizerInfo_t mini
        cdef object          hashes  = array.array("L")
        cdef object          ids     = array.array("l")
        cdef object          offsets = array.array("l")
        cdef size_t          length  = 0
        if self._minimizers != NULL:
            length = self._minimizers.size()
            for mini in self._minimizers[0]:
                hashes.append(mini.hash)
                ids.append(mini.seqId)
                offsets.append(mini.wpos)
        return {
            "hashes": hashes,
            "ids": ids,
            "offsets": offsets,
            "length": length,
        }

    cpdef object __setstate__(self, dict state):
        cdef size_t i
        cdef size_t length  = state["length"]
        cdef object hashes  = state["hashes"]
        cdef object ids     = state["ids"]
        cdef object offsets = state["offsets"]
        if self._minimizers == NULL:
            self._minimizers = new vector[MinimizerInfo_t]()
        self._minimizers.resize(length)
        for i, (hash, seqId, wpos) in enumerate(zip(hashes, ids, offsets)):
            self._minimizers[0][i].hash = hash
            self._minimizers[0][i].seqId = seqId
            self._minimizers[0][i].wpos = wpos


cdef class Hit:
    """A single hit found when querying a `Mapper` with a genome.

    Attributes:
        name (`object`): The name of the genome that produced a hit, as
            given to `Sketch.add_genome` or `Sketch.add_draft`.
        matches (`int`): The number of fragments that matched the target
            genome.
        fragments (`int`): The total number of fragments used to compare
            the query and target genomes.
        identity (`float`): The average nucleotide identity between the
            two genomes, given as a percentage.

    """

    # --- Attributes ---------------------------------------------------------

    cdef readonly object  name
    cdef readonly seqno_t matches
    cdef readonly seqno_t fragments
    cdef readonly float   identity

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object name, float identity, seqno_t matches, seqno_t fragments):
        """__init__(self, name, identity, matches, fragments)\n--

        Create a new `Hit` instance with the given parameters.

        """
        self.name = name
        self.matches = matches
        self.fragments = fragments
        self.identity = identity

    def __repr__(self):
        cdef str ty = type(self).__name__
        return "{}(name={!r}, identity={!r}, matches={!r}, fragments={!r})".format(
            ty, self.name, self.identity, self.matches, self.fragments
        )

    def __eq__(self, Hit other):
        return (
                self.name == other.name
            and self.matches == other.matches
            and self.fragments == other.fragments
            and self.identity == other.identity
        )

    def __reduce__(self):
        return (
            Hit,
            (self.name, self.identity, self.matches, self.fragments)
        )


cdef class MinimizerInfo:
    """The information about a single minimizer.
    """

    # --- Attributes ---------------------------------------------------------

    cdef readonly hash_t   hash
    cdef readonly seqno_t  sequence_id
    cdef readonly offset_t window_position

    # --- Magic methods ------------------------------------------------------

    def __init__(self, hash_t hash, seqno_t sequence_id, offset_t window_position):
        """__init__(self, hash, sequence_id, window_position)\n--

        Create a new `MinimizerInfo` with the given parameters.

        """
        self.hash = hash
        self.sequence_id = sequence_id
        self.window_position = window_position

    def __repr__(self):
        cdef str ty = type(self).__name__
        return "{}(hash={!r}, sequence_id={!r}, window_position={!r})".format(
            ty, self.hash, self.sequence_id, self.window_position
        )

    def __eq__(self, MinimizerInfo other):
        return (
                self.hash == other.hash
            and self.sequence_id == other.sequence_id
            and self.window_position == other.window_position
        )

    def __reduce__(self):
        return (
            MinimizerInfo,
            (self.hash, self.sequence_id, self.window_position)
        )

    # --- Methods ------------------------------------------------------------

    @staticmethod
    cdef MinimizerInfo from_raw(MinimizerInfo_t raw):
        return MinimizerInfo(raw.hash, raw.seqId, raw.wpos)

    cdef MinimizerInfo_t to_raw(self):
        cdef MinimizerInfo_t info
        info.hash = self.hash
        info.seqId = self.sequence_id
        info.wpos = self.window_position
        return info


cdef class Position:

    # --- Attributes ---------------------------------------------------------

    cdef readonly seqno_t  sequence_id
    cdef readonly offset_t window_position

    # --- Magic methods ------------------------------------------------------

    def __init__(self, seqno_t sequence_id, offset_t window_position):
        """__init__(self, sequence_id, window_position)\n--

        Create a new `Position` instance with the given parameters.

        """
        self.sequence_id = sequence_id
        self.window_position = window_position

    def __repr__(self):
        cdef str ty = type(self).__name__
        return "{}(sequence_id={!r}, window_position={!r})".format(
            ty, self.sequence_id, self.window_position
        )

    def __eq__(self, Position other):
        return (
                self.sequence_id == other.sequence_id
            and self.window_position == other.window_position
        )

    def __reduce__(self):
        return (
            Position,
            (self.sequence_id, self.window_position)
        )

    # --- Methods ------------------------------------------------------------

    @staticmethod
    cdef Position from_raw(MinimizerMetaData_t raw):
        return Position(raw.seqId, raw.wpos)

    cdef MinimizerMetaData_t to_raw(self):
        cdef MinimizerMetaData_t m
        m.seqId = self.sequence_id
        m.wpos = self.window_position
        return m


cdef class MinimizerIndex:
    """The index mapping minimizer hash values to their positions.
    """

    # --- Attributes ---------------------------------------------------------

    cdef Sketch_t.MI_Map_t* _map
    cdef object             owner

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        # create a new Sketch with the parameters
        self._map = NULL
        self.owner = None

    def __init__(self):
        self._map = new Sketch_t.MI_Map_t()

    def __dealloc__(self):
        if self.owner is None:
            del self._map

    def __len__(self):
        assert self._map != NULL
        return self._map.size()

    def __iter__(self):
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].iterator it
        it = self._map.begin()
        while it != self._map.end():
            yield dereference(it).first
            preincrement(it)

    def __contains__(self, hash_t item):
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].const_iterator find_it
        find_it = self._map.const_find(item)
        return find_it != self._map.end()

    def __getitem__(self, hash_t item):
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].iterator find_it
        find_it = self._map.find(item)
        if find_it != self._map.end():
            return [Position.from_raw(x) for x in dereference(find_it).second]
        raise KeyError(item)

    def __setitem__(self, hash_t item, object value):
        # convert value to vector
        cdef Position position
        cdef MinimizerMapValueType_t positions = MinimizerMapValueType_t()
        for position in value:
            positions.push_back(position.to_raw())
        # remove previous element if any
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].iterator find_it
        find_it = self._map.find(item)
        if find_it != self._map.end():
            self._map.erase(find_it)
        # add new element
        cdef pair[MinimizerMapKeyType_t, MinimizerMapValueType_t] element
        element.first = item
        element.second = positions
        self._map.insert(element)

    def __delitem__(self, hash_t item):
      # remove previous element if any
      cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].iterator find_it
      find_it = self._map.find(item)
      if find_it != self._map.end():
          self._map.erase(find_it)
      else:
          raise KeyError(item)

    def __reduce__(self):
        return (
            MinimizerIndex,
            (),
            None,
            None,
            self.items(),
        )

    # --- Methods ------------------------------------------------------------

    def items(self):
        cdef MinimizerMetaData_t     pos_raw
        cdef MinimizerMapKeyType_t   key
        cdef MinimizerMapValueType_t value
        cdef size_t                  i
        cdef list                    positions
        cdef Position                position
        cdef unordered_map[MinimizerMapKeyType_t, MinimizerMapValueType_t].iterator it

        it = self._map.begin()
        while it != self._map.end():

            key = dereference(it).first
            value = dereference(it).second

            i = 0
            positions = PyList_New(value.size())

            for pos_raw in value:
                pos = Position.from_raw(pos_raw)
                Py_INCREF(pos)
                PyList_SET_ITEM(positions, i, pos)
                i += 1

            yield key, positions
            preincrement(it)
