# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp


# --- C imports --------------------------------------------------------------

cimport libcpp11.chrono
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
    # MappingResult,
    MappingResultsVector_t,
    QueryMetaData as QueryMetaData_t,
)

from _utils cimport kseq_ptr_t


# --- Python imports ---------------------------------------------------------

import os
import warnings


# --- Cython classes ---------------------------------------------------------

cdef class Mapper:

    cdef size_t           _counter
    cdef Sketch_t*        _sk
    cdef bool             _indexed
    cdef vector[uint64_t] _lengths
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

        #
        self._indexed = True
        self._counter = 0

    def __dealloc__(self):
        del self._sk

    # --- Properties ---------------------------------------------------------

    @property
    def occurences_threshold(self):
        """`int`: The occurence threshold above which minimizers are ignored.
        """
        assert self._sk != NULL
        return self._sk.freqThreshold

    @property
    def percentage_threshold(self):
        """`float`: The fraction of most frequent minimizers to ignore.
        """
        assert self._sk != NULL
        return self._sk.percentageThreshold

    @property
    def names(self):
        """`list` of `str`: The names of the sequences currently stored.
        """
        assert self._sk != NULL
        return [
            contig_info.name.decode("utf-8")
            for contig_info in self._sk.metadata
        ]

    # --- Methods ------------------------------------------------------------

    cdef int _add_sequence(self, bytes name, const unsigned char[::1] seq) except +:
        assert self._sk != NULL
        assert self.parameters is not None

        cdef kseq_t        kseq
        cdef ContigInfo_t  info
        cdef size_t        seql  = len(seq)
        cdef Parameters_t* param = &self._param

        # store the sequence name and size
        info.name = string(name)
        info.len = seql
        self._sk.metadata.push_back(info)

        # check the sequence is large enough to compute minimizers
        if seql >= param.windowSize and seql >= param.kmerSize:
            # WARNING: Normally kseq_t owns the buffer, but here we just give
            #          a view to avoid reallocation if possible. However,
            #          `addMinimizers` will always attempt to make the
            #          sequence uppercase, so it should be checked in
            #          advance that the sequence is already uppercase.
            kseq.seq.l = kseq.seq.m = seql
            kseq.seq.s = <char*> <const unsigned char*> &seq[0]
            # compute the minimizers
            addMinimizers(
                self._sk.minimizerIndex,
                &kseq,
                param.kmerSize,
                param.windowSize,
                param.alphabetSize,
                self._counter
            );

        # record the reference sequence name and length
        self._param.refSequences.push_back(<string> name)
        self._lengths.push_back(seql)

        # this genome only contained a single sequence
        self._sk.sequencesByFileInfo.push_back(self._counter)

        # the Sketch will need to be reindexed since we added a new sequence
        self._indexed = False
        self._counter += 1

    def add_sequence(self, str name, str sequence):
        """Add a sequence to the sketcher.

        The sequence will be considered like a complete genome. If
        you would like to create a reference genome from several
        contigs instead, use the `Sketch.add_fragments` method.

        Arguments:
            name (`str`): The name of the sequence to add.
            sequence (`bytes`): The sequence to add.

        Warnings:
            `UserWarning`: When the sequence is too short for minimizers
                to be computed for it.

        """
        assert self.parameters is not None

        cdef Parameters_t* p = &self._param
        if len(sequence) < p.windowSize or len(sequence) < p.kmerSize:
            warnings.warn(UserWarning, (
                "Sketch received a short sequence relative to parameters, "
                "minimizers will not be added."
            ))

        self._add_sequence(name.encode("utf-8"), sequence.upper().encode("utf-8"))

    def add_fragments(self, str name, object fragments):
        raise NotImplementedError("Sketch.add_fragments")

    def index(self):
        """Build the index for fast lookups using minimizer table.
        """
        if not self._indexed:
            # clear the maps in case we are rebuilding over a previous index
            self._sk.minimizerPosLookupIndex.clear()
            self._sk.minimizerFreqHistogram.clear()
            self._sk.freqThreshold = INT_MAX
            # compute the index and the frequency histogram
            self._sk.index()
            self._sk.computeFreqHist()
            # mark the sketch as indexed
            self._indexed = True

    def is_indexed(self):
        """Check whether or not this `Sketch` object has been indexed.
        """
        return self._indexed

    cdef MappingResultsVector_t _query_fragment(self, int i, seqno_t seq_counter, const unsigned char[::1] seq, int min_read_length, Map_t* map) nogil:

        cdef kseq_t                 kseq_buffer
        cdef MappingResultsVector_t l2_mappings = MappingResultsVector_t()
        cdef ofstream               out         = ofstream()

        cdef QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type] query
        query.kseq = &kseq_buffer
        query.kseq.seq.s = <char*> &seq[i * min_read_length]
        query.kseq.seq.l = query.kseq.seq.m = min_read_length
        query.seqCounter = seq_counter + i

        map.mapSingleQuerySeq[QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type]](query, l2_mappings, out)

        return l2_mappings

    cdef object _query_sequence(self, const unsigned char[::1] seq):
        """Query the sketcher for the given sequence.

        Adapted from the ``skch::Map::mapQuery`` method in ``computeMap.hpp``.

        Todo:
            Doctor the memory management so that it's not needed to
            reallocate a new ``Map_t`` and a new result vector at each
            query.

        """
        assert self._sk != NULL
        assert self.parameters is not None

        cdef int                    i
        cdef int                    fragment_count
        cdef Map_t*                 map
        cdef MappingResultsVector_t l2_mappings
        cdef MappingResultsVector_t final_mappings
        cdef vector[CGI_Results]    results
        cdef uint64_t               seql                  = len(seq)
        cdef uint64_t               total_query_fragments = 0
        cdef Parameters_t           p                     = self._param

        if seql >= p.windowSize and seql >= p.kmerSize and seql >= p.minReadLength:
            # create a new mapper with the given mapping result vector
            final_mappings = MappingResultsVector_t()
            map = new Map_t(p, self._sk[0], total_query_fragments, 0)
            # compute the expected number of blocks
            fragment_count = seql // p.minReadLength
            total_query_fragments = fragment_count
            # map the blocks
            for i in range(fragment_count):
                l2_mappings = self._query_fragment(i, 0, seq, p.minReadLength, map)
                for m in l2_mappings:
                    final_mappings.push_back(m)
            # compute core genomic identity after successful mapping
            computeCGI(
                p,
                final_mappings,
                map[0],
                self._sk[0],
                total_query_fragments, # total query fragments, here equal to fragment count
                0, # queryFileNo, only used for visualization, ignored
                string(), # fileName, only used for reporting, ignored
                results,
            )
            # free the map
            del map

        # WIP: return a proper object here
        for res in results:
            min_length = min(seql, self._lengths[res.refGenomeId-1])
            shared_length = res.countSeq * p.minReadLength
            if shared_length >= min_length * p.minFraction:
                print({
                    "name": self._param.refSequences[res.refGenomeId-1].decode("utf-8"),
                    "identity": res.identity,
                    "matches": res.countSeq,
                    "fragments": res.totalQueryFragments,
                })

    def query_sequence(self, str sequence):

        cdef Parameters_t* p = &self._param
        if len(sequence) < p.windowSize or len(sequence) < p.kmerSize or len(sequence) < p.minReadLength:
            warnings.warn(UserWarning, (
                "Sketch received a short sequence relative to parameters, "
                "mapping will not be computed."
            ))

        self._query_sequence(sequence.upper().encode("utf-8"))
