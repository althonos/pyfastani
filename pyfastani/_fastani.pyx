# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp
"""Bindings to FastANI, a method for fast whole-genome similarity estimation.
"""

# --- C imports --------------------------------------------------------------

cimport libcpp11.chrono
from cpython cimport PyObject
from libc.limits cimport INT_MAX
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, realloc, free
from libcpp cimport bool, nullptr
from libcpp.functional cimport function
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector
from libcpp11.fstream cimport ofstream

cimport fastani.map.map_stats
from kseq cimport kseq_t, kstring_t
from fastani.cgi.compute_core_identity cimport computeCGI
from fastani.cgi.cgid_types cimport CGI_Results
from fastani.map cimport base_types
from fastani.map.compute_map cimport Map as Map_t
from fastani.map.common_func cimport addMinimizers
from fastani.map.map_parameters cimport Parameters as Parameters_t
from fastani.map.map_stats cimport recommendedWindowSize
from fastani.map.win_sketch cimport Sketch as Sketch_t
from fastani.map.base_types cimport (
    seqno_t,
    ContigInfo as ContigInfo_t,
    MappingResultsVector_t,
    QueryMetaData as QueryMetaData_t,
)

# HACK: we need kseq_t* as a template argument, which is not supported by
#       Cython at the moment, so we just `typedef kseq_t* kseq_ptr_t` in
#       an external C++ header to make Cython happy
from _utils cimport kseq_ptr_t


# --- Python imports ---------------------------------------------------------

import warnings


# --- Cython helpers ---------------------------------------------------------

cdef bool isupper(const unsigned char[::1] seq):
    """Check if a buffer contains only upper case characters
    """
    cdef size_t i
    for i in range(seq.shape[0]):
        if seq[i] >= b'a' and seq[i] <= b'z':
            return False
    return True


cdef void upper(unsigned char[::1] seq):
    """Make the letters in a buffer upper case, in place.
    """
    cdef size_t i
    for i in range(seq.shape[0]):
        if seq[i] >= b'a' and seq[i] <= b'z':
            seq[i] -= 32


# --- Cython classes ---------------------------------------------------------

cdef class Sketch:
    """An index computing minimizers over the reference genomes.
    """

    cdef size_t           _counter
    cdef Sketch_t*        _sk
    cdef vector[uint64_t] _lengths
    cdef list             _names
    cdef Parameters_t     _param

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        # hardcode reporting parameters so that we can control
        # execution flow
        self._param.alphabetSize = 4
        self._param.threads = 1
        self._param.reportAll = True
        self._param.visualize = False
        self._param.matrixOutput = False
        # create a new Sketch with the parameters
        self._sk = new Sketch_t(self._param)

    def __init__(
        self,
        *,
        unsigned int k=16,
        unsigned int fragment_length=3000,
        double minimum_fraction=0.2,
        double p_value=1e-03,
        double percentage_identity=80.0,
        unsigned int reference_size=5_000_000,
    ):
        """__init__(*, k=16, fragment_length=3000, minimum_fraction=0.2, p_value=1e-03, percentage_identity=80, reference_size=5000000)\n--

        Create a new FastANI sequence mapper.

        Keyword Arguments:
            k (`int`): The size of the k-mers. FastANI authors recommend
                a size of at most 16, but any positive number should work.
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

        """
        if minimum_fraction > 1 or minimum_fraction < 0:
            raise ValueError(f"minimum_fraction must be between 0 and 1, got {minimum_fraction!r}")
        if fragment_length <= 0:
            raise ValueError(f"fragment_length must be strictly positive, got {fragment_length!r}")
        if p_value <= 0:
            raise ValueError(f"p_value must be positive, got {p_value!r}")
        if percentage_identity > 100 or percentage_identity < 0:
            raise ValueError(f"percentage_identity must be between 0 and 100, got {percentage_identity!r}")
        if k > 16:
            warnings.warn(
                UserWarning,
                f"Using k-mer size greater than 16 ({k!r}), accuracy of the results is not guaranteed."
            )

        # store parameters
        self._param.kmerSize = k
        self._param.minReadLength = fragment_length
        self._param.minFraction = minimum_fraction
        self._param.p_value = p_value
        self._param.percentageIdentity = percentage_identity
        self._param.referenceSize = reference_size

        # compute the recommended window size
        self._param.windowSize = fastani.map.map_stats.recommendedWindowSize(
            self._param.p_value,
            self._param.kmerSize,
            self._param.alphabetSize,
            self._param.percentageIdentity,
            self._param.minReadLength,
            self._param.referenceSize
        )

        # initialize bookkeeping values
        self._names = []
        self._lengths = vector[uint64_t]()
        self._counter = 0

    def __dealloc__(self):
        del self._sk

    # --- Properties ---------------------------------------------------------

    @property
    def occurences_threshold(self):
        """`int`: The occurence threshold above which minimizers are ignored.
        """
        assert self._sk != nullptr
        return self._sk.freqThreshold

    @property
    def percentage_threshold(self):
        """`float`: The fraction of most frequent minimizers to ignore.
        """
        assert self._sk != nullptr
        return self._sk.percentageThreshold

    @property
    def names(self):
        """`list` of `str`: The names of the sequences currently sketched.
        """
        return self._names[:]

    # --- Methods (adding references) ----------------------------------------

    cdef int _add_draft(self, object name, object contigs) except 1:
        """Add a draft genome to the sketcher.

        Adapted from the ``skch::Sketch::build`` method in ``winSketch.hpp``
        to work without the ``kseq`` I/O.

        """
        assert self._sk != nullptr

        cdef const unsigned char[::1] seq
        cdef kseq_t                   kseq
        cdef ContigInfo_t             info
        cdef size_t                   seql_total = 0
        cdef Parameters_t*            param      = &self._param

        for contig in contigs:
            # encode the sequence if needed
            if isinstance(contig, str):
                contig = contig.encode("ascii")
            # make sure the sequence is uppercase, otherwise a `makeUpperCase`
            # is going to write in our read-only buffer in `addMinimizers`
            if not isupper(contig):
                warnings.warn(
                    UserWarning,
                    "Contig contains lowercase characters, reallocating."
                )
                contig = bytearray(contig)
                upper(contig)

            # get a memory view of the sequence
            seq = contig

            # store the contig name and size
            # (we just use the same name for all contigs)
            # <-- actually useless, only useful for reporting
            # info.name = string(name)
            # info.len = seq.shape[0]
            # self._sk.metadata.push_back(info)

            # check the sequence is large enough to compute minimizers
            if seq.shape[0] >= param.windowSize and seq.shape[0] >= param.kmerSize:
                with nogil:
                    # WARNING: Normally kseq_t owns the buffer, but here we just give
                    #          a view to avoid reallocation if possible. However,
                    #          `addMinimizers` will always attempt to make the
                    #          sequence uppercase, so it should be checked in
                    #          advance that the sequence is already uppercase.
                    kseq.seq.l = kseq.seq.m = seq.shape[0]
                    kseq.seq.s = <char*> <const unsigned char*> &seq[0]
                    # compute the minimizers
                    addMinimizers(
                        self._sk.minimizerIndex,
                        &kseq,
                        param.kmerSize,
                        param.windowSize,
                        param.alphabetSize,
                        self._counter
                    )
            else:
                warnings.warn(UserWarning, (
                    "Sketch received a short contig relative to parameters, "
                    "minimizers will not be added."
                ))

            # compute the genome length with respect to the fragment length
            seql_total += (seq.shape[0] // param.minReadLength) * param.minReadLength

            # record the number of contigs currently stored
            self._counter += 1

        # record the reference genome name and total length
        self._names.append(name)
        self._lengths.push_back(seql_total)

        # record the number of contigs in this genome
        self._sk.sequencesByFileInfo.push_back(self._counter)

    cdef Sketch add_draft(self, object name, object contigs):
        """add_draft(self, name, contigs)\n--

        Add a reference draft genome to the sketcher.

        Using this method is fine even when the genome has a single contig,
        although `Mapper.add_genome` is easier to use in that case.

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

    cdef Sketch add_genome(self, object name, object sequence):
        """add_genome(self, name, sequence)\n--

        Add a reference genome to the sketcher.

        This method is a shortcut for `Mapper.add_draft` when a genome is
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

    # --- Methods (indexing) -------------------------------------------------

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
        mapper._names = self._names
        mapper._lengths.swap(self._lengths)

        # reset the current sketch
        self._names = []
        self._sk = new Sketch_t(self._param)

        # return the new mapper
        return mapper


cdef class Mapper:
    """A genome mapper using Murmur3 hashes and k-mers to compute ANI.
    """

    cdef Sketch_t*        _sk
    cdef vector[uint64_t] _lengths
    cdef list             _names
    cdef Parameters_t     _param

    def __cinit__(self):
        self._sk = NULL

    def __dealloc__(self):
        del self._sk

    @staticmethod
    cdef void _query_fragment(
        int i,
        seqno_t
        seq_counter,
        const unsigned char[::1] seq,
        int min_read_length,
        Map_t* map,
        kseq_t* kseq,
        ofstream* out,
        MappingResultsVector_t* mappings
    ) nogil:
        cdef QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type] query
        query.kseq = kseq
        query.kseq.seq.s = <char*> &seq[i * min_read_length]
        query.kseq.seq.l = query.kseq.seq.m = min_read_length
        query.seqCounter = seq_counter + i
        map.mapSingleQuerySeq[QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type]](query, mappings[0], out[0])

    cdef object _query_draft(self, object contigs):
        """Query the sketcher for the given contigs.

        Adapted from the ``skch::Map::mapQuery`` method in ``computeMap.hpp``.

        """
        assert self._sk != nullptr

        cdef int                      i               # fragment counter
        cdef kseq_t                   kseq            #
        cdef Map_t*                   map
        cdef MappingResultsVector_t   final_mappings
        cdef vector[CGI_Results]      results
        cdef object                   contig
        cdef const unsigned char[::1] seq
        cdef uint64_t                 seql
        cdef int                      fragment_count  = 0
        cdef int                      seq_counter     = 0
        cdef uint64_t                 total_fragments = 0
        cdef Parameters_t             p               = self._param
        cdef list                     hits            = []
        cdef ofstream                 out             = ofstream()

        # create a new mapper with the given mapping result vector
        map = new Map_t(p, self._sk[0], total_fragments, 0)

        # iterate over contigs
        for contig in contigs:
            # make sure the sequence is bytes
            if isinstance(contig, str):
                contig = contig.encode("ascii")
            # make sure the sequence is uppercase, otherwise a `makeUpperCase`
            # is going to write in our read-only buffer in `addMinimizers`
            if not isupper(contig):
                warnings.warn(
                    UserWarning,
                    "Sequence contains lowercase characters, reallocating."
                )
                contig = bytearray(contig)
                upper(contig)
            # get a memory view of the sequence
            seq = contig
            seql = seq.shape[0]

            # query if the sequence is large enough
            if seql >= p.windowSize and seql >= p.kmerSize and seql >= p.minReadLength:
                with nogil:
                    # compute the expected number of blocks
                    fragment_count = seql // p.minReadLength
                    # map the blocks
                    for i in range(fragment_count):
                        Mapper._query_fragment(
                            i,
                            total_fragments,
                            seq,
                            p.minReadLength,
                            map,
                            &kseq,
                            &out,
                            &final_mappings
                        )
                # record the number of fragments
                total_fragments += fragment_count
            else:
                fragmentCount = 0
                warnings.warn(UserWarning, (
                    "Mapper received a short sequence relative to parameters, "
                    "mapping will not be computed."
                ))

        # compute core genomic identity after successful mapping
        computeCGI(
            p,
            final_mappings,
            map[0],
            self._sk[0],
            total_fragments, # total query fragments
            0, # queryFileNo, only used for visualization, ignored
            string(), # fileName, only used for reporting, ignored
            results,
        )
        # free the map
        del map
        # build and return the list of hits
        for res in results:
            assert res.refGenomeId < self._lengths.size()
            assert res.refGenomeId < len(self._names)
            min_length = min(seql, self._lengths[res.refGenomeId])
            shared_length = res.countSeq * p.minReadLength
            if shared_length >= min_length * p.minFraction:
                hits.append(Hit(
                    name=self._names[res.refGenomeId],
                    identity=res.identity,
                    matches=res.countSeq,
                    fragments=res.totalQueryFragments,
                ))
        return hits

    def query_draft(self, object contigs):
        """query_draft(self, contigs)\n--

        Query the mapper for a complete genome.

        Arguments:
            contigs (iterable or `str` or `bytes`): The genome to query the mapper
                with.

        Note:
            Sequence must be larger than the window size, the k-mer size,
            and the fragment length to be mapped, otherwise an empty list
            of hits will be returned.

        Returns:
            `list` of `~pyfastani.Hit`: The hits found for the query.

        """
        # delegate to C code
        return self._query_draft(contigs)

    def query_genome(self, object sequence):
        """query_genome(self, sequence)\n--

        Query the mapper for a complete genome.

        Arguments:
            sequence (`str` or `bytes`): The genome to query the mapper
                with.

        Note:
            Sequence must be larger than the window size, the k-mer size,
            and the fragment length to be mapped, otherwise an empty list
            of hits will be returned.

        Returns:
            `list` of `~pyfastani.Hit`: The hits found for the query.

        """
        # delegate to C code
        return self._query_draft((sequence,))


cdef class Hit:
    """A single hit found when querying the mapper with a genome.
    """
    cdef readonly object  name
    cdef readonly seqno_t matches
    cdef readonly seqno_t fragments
    cdef readonly float   identity

    def __init__(self, object name, float identity, seqno_t matches, seqno_t fragments):
        """__init__(self, name, identity, matches, fragments)\n--

        Create a new `Hit` instance with the given parameters.

        """
        self.name = name
        self.matches = matches
        self.fragments = fragments
        self.identity = identity

    def __repr__(self):
        ty = type(self).__name__
        return "{}(name={!r}, identity={!r} matches={!r}, fragments={!r})".format(
            ty, self.name, self.identity, self.matches, self.fragments
        )

    def __eq__(self, Hit other):
        return (
                self.name == other.name
            and self.matches == other.matches
            and self.fragments == other.fragments
            and self.identity == other.identity
        )
