from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp11.chrono cimport high_resolution_clock


cdef extern from "map/include/base_types.hpp" namespace "skch" nogil:

    ctypedef uint32_t hash_t
    ctypedef int      offset_t
    ctypedef int      seqno_t

    # ctypedef high_resolution_clock Time

    cdef cppclass MinimizerInfo:
        hash_t   hash
        seqno_t  seqId
        offset_t wpos

        bool operator<  (const MinimizerInfo &x)
        bool operator== (const MinimizerInfo &x)
        bool operator!= (const MinimizerInfo &x)

        @staticmethod
        bool equalityByHash(const MinimizerInfo &x, const MinimizerInfo &y)
        @staticmethod
        bool lessByHash    (const MinimizerInfo &x, const MinimizerInfo &y)


    cdef cppclass MinimizerMetaData:
        seqno_t  seqId
        offset_t wpos


    ctypedef hash_t                    MinimizerMapKeyType
    ctypedef vector[MinimizerMetaData] MinimizerMapValueType;


    cdef cppclass ContigInfo:
        string   name
        offset_t len


    cdef cppclass QueryMetaData[ K, MV ]:
        K       kseq
        seqno_t seqCounter
        int     sketchSize
        MV      minimizerTableQuery


    cdef cppclass MappingResult:
        offset_t queryLen
        offset_t refStartPos
        offset_t refEndPos
        offset_t queryStartPos
        offset_t queryEndPos
        seqno_t  refSeqId
        seqno_t  querySeqId
        float    nucIdentity
        float    nucIdentityUpperBound
        int      sketchSize
        int      conservedSketches


    ctypedef vector[MappingResult] MappingResultsVector_t;
