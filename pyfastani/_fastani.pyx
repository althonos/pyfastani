# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp, binding=False
"""Bindings to FastANI, a method for fast whole-genome similarity estimation.
"""

# --- C imports --------------------------------------------------------------

cimport libcpp11.chrono
from cython.operator cimport dereference, preincrement, postincrement
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

cimport fastani.map.map_stats
from kseq cimport kseq_t, kstring_t
from fastani.cgi.compute_core_identity cimport computeCGI
from fastani.cgi.cgid_types cimport CGI_Results
from fastani.map cimport base_types
from fastani.map.compute_map cimport Map as Map_t
from fastani.map.common_func cimport addMinimizers, getHash
from fastani.map.map_parameters cimport Parameters as Parameters_t
from fastani.map.map_stats cimport estimateMinimumHitsRelaxed, recommendedWindowSize
from fastani.map.win_sketch cimport Sketch as Sketch_t
from fastani.map.base_types cimport (
    hash_t,
    seqno_t,
    ContigInfo as ContigInfo_t,
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


# --- Python imports ---------------------------------------------------------

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


cdef int _add_minimizers(
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


# --- Cython classes ---------------------------------------------------------


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


cdef class Sketch(_Parameterized):
    """An index computing minimizers over the reference genomes.
    """

    # --- Attributes ---------------------------------------------------------

    cdef Sketch_t*        _sk       # the internal Sketch_t
    # cdef Parameters_t     _param    # the internal Parameters_t (const)
    cdef size_t           _counter  # the number of contigs (not genomes)
    cdef vector[uint64_t] _lengths  # array mapping each genome to its length
    cdef list             _names    # list mapping each genome to its name

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        # hardcode reporting parameters so that we can control
        # execution flow
        self._param.alphabetSize = 4

        # create a new Sketch with the parameters
        self._sk = new Sketch_t(self._param)
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
    ):
        f"""__init__(*, k=16, fragment_length=3000, minimum_fraction=0.2, p_value=1e-03, percentage_identity=80, reference_size=5000000)\n--

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

        # compute the recommended window size
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
                assert param.alphabetSize == 4 # assumed in `_add_minimizers`
                with nogil:
                    _add_minimizers(
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

    cpdef Sketch add_genome(self, object name, object sequence):
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
        # reset the current sketch
        self._sk = new Sketch_t(self._param)
        self.clear()
        # return the new mapper
        return mapper


cdef class Mapper(_Parameterized):
    """A genome mapper using Murmur3 hashes and k-mers to compute ANI.
    """

    # --- Attributes ---------------------------------------------------------

    cdef Sketch_t*        _sk
    cdef vector[uint64_t] _lengths
    cdef list             _names
    # cdef Parameters_t     _param

    # --- Magic methods ------------------------------------------------------

    def __dealloc__(self):
        del self._sk

    @staticmethod
    cdef void _do_l1_mappings(
        const Parameters_t& param,
        const Sketch_t& ref_sketch,
        Map_t& map,
        const int kind,
        const void* data,
        const ssize_t slen,
        QueryMetaData_t[kseq_ptr_t, vector[MinimizerInfo_t]]& query,
        vector[Map_t.L1_candidateLocus_t]& l1_mappings,
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
        _add_minimizers(
           query.minimizerTableQuery,
           kind,
           data,
           slen,
           param.kmerSize,
           param.windowSize,
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
            seed_find = ref_sketch.minimizerPosLookupIndex.const_find(dereference(it).hash)
            if seed_find != ref_sketch.minimizerPosLookupIndex.end():
                hit_position_list = dereference(seed_find).second
                if hit_position_list.size() < ref_sketch.getFreqThreshold():
                    seed_hits_l1.insert(seed_hits_l1.end(), hit_position_list.begin(), hit_position_list.end())
            postincrement(it)

        # estimate the number of minimum hits, and compute candidates
        minimum_hits = estimateMinimumHitsRelaxed(query.sketchSize, param.kmerSize, param.percentageIdentity)
        map.computeL1CandidateRegions(query, seed_hits_l1, minimum_hits, l1_mappings)

    @staticmethod
    cdef void _query_fragment(
        const Parameters_t& param,
        const Sketch_t& sketch,
        Map_t& map,
        const int i,
        const seqno_t seq_counter,
        const int kind,
        const void* data,
        const ssize_t slen,
        const size_t stride,
        MappingResultsVector_t& l2_mappings
    ) nogil:
        cdef kseq_t                                         kseq
        cdef QueryMetaData_t[kseq_ptr_t, Map_t.MinVec_Type] query
        cdef vector[Map_t.L1_candidateLocus_t]              l1_mappings

        query.kseq = &kseq
        query.kseq.seq.s = NULL
        query.kseq.seq.l = param.minReadLength
        query.seqCounter = seq_counter + i

        Mapper._do_l1_mappings(
            # classes with configuration
            param,
            sketch,
            map,
            # start of the i-th fragment
            kind,
            <const void*> ((<const char*> data) + (i * param.minReadLength) * stride),
            param.minReadLength,
            # outputs
            query,
            l1_mappings,
        )

        map.doL2Mapping(query, l1_mappings, l2_mappings)


    # --- Methods ------------------------------------------------------------

    cdef list _query_draft(self, object contigs):
        """Query the sketcher for the given contigs.

        Adapted from the ``skch::Map::mapQuery`` method in ``computeMap.hpp``.

        """
        assert self._sk != nullptr

        # iterators over the contigs
        cdef int                      i               # fragment counter
        cdef object                   contig
        cdef int                      fragment_count
        # bookkeeping
        cdef uint64_t                 total_fragments = 0
        cdef uint64_t                 total_length    = 0
        # mapping parameters and reference
        cdef Map_t*                   map
        cdef Parameters_t*            param           = &self._param
        cdef Sketch_t*                sketch          = self._sk
        # sequence as a unicode object
        cdef const unsigned char[::1] view
        cdef int                      kind
        cdef void*                    data
        cdef ssize_t                  slen
        cdef size_t                   stride
        # core genomic identity results
        cdef MappingResultsVector_t   final_mappings
        cdef vector[CGI_Results]      results
        cdef CGI_Results              result
        # filtering and reporing results
        cdef uint64_t                 min_length
        cdef uint64_t                 shared_length
        cdef list                     hits            = []

        # create a new mapper with the given mapping result vector
        map = new Map_t(param[0], sketch[0], total_fragments, 0)

        # iterate over contigs
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

            # query if the sequence is large enough
            if slen >= param.windowSize and slen >= param.kmerSize and slen >= param.minReadLength:
                with nogil:
                    # compute the expected number of blocks
                    fragment_count = slen // param.minReadLength
                    # map the blocks
                    for i in range(fragment_count):
                        Mapper._query_fragment(
                            # classes with configuration
                            param[0],
                            sketch[0],
                            map[0],
                            # position in query
                            i,
                            total_fragments,
                            # unicode sequence
                            kind,
                            data,
                            slen,
                            stride,
                            # result vector
                            final_mappings
                        )
                    # record the number of fragments
                    total_fragments += fragment_count
                    total_length += slen
            else:
                warnings.warn(
                    (
                        "Mapper received a short sequence relative to parameters, "
                        "mapping will not be computed."
                    ),
                    UserWarning,
                )

        # compute core genomic identity after successful mapping
        with nogil:
            computeCGI(
                param[0],
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
        for result in results:
            assert result.refGenomeId < self._lengths.size()
            assert result.refGenomeId < len(self._names)
            min_length = min(total_length, self._lengths[result.refGenomeId])
            shared_length = result.countSeq * param.minReadLength
            if shared_length >= min_length * param.minFraction:
                hits.append(Hit(
                    name=self._names[result.refGenomeId],
                    identity=result.identity,
                    matches=result.countSeq,
                    fragments=result.totalQueryFragments,
                ))
        return hits

    cpdef list query_draft(self, object contigs):
        """query_draft(self, contigs)\n--

        Query the mapper for a complete genome.

        Arguments:
            contigs (iterable or `str` or `bytes`): The genome to query the mapper
                with.

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

        """
        # delegate to C code
        return self._query_draft(contigs)

    cpdef list query_genome(self, object sequence):
        """query_genome(self, sequence)\n--

        Query the mapper for a complete genome.

        Arguments:
            sequence (`str` or `bytes`): The genome to query the mapper
                with.

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
        cdef str ty = type(self).__name__
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
