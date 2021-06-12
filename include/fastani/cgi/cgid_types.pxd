from libcpp cimport bool

from fastani.map.base_types cimport seqno_t, offset_t


cdef extern from "cgi/include/cgid_types.hpp" namespace "cgi" nogil:

    cdef cppclass MappingResult_CGI:

        seqno_t refSequenceId
        seqno_t genomeId
        seqno_t querySeqId
        seqno_t refStartPos
        seqno_t queryStartPos
        seqno_t mapRefPosBin
        float   nucIdentity


    cdef cppclass CGI_Results:

        seqno_t refGenomeId
        seqno_t qryGenomeId
        seqno_t countSeq
        seqno_t totalQueryFragments
        float   identity

        bool  operator<(const CGI_Results& x)
