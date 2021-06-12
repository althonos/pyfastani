from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.functional cimport function

from fastani.map.base_types cimport seqno_t, offset_t, ContigInfo, MappingResult, MappingResultsVector_t
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch


cdef extern from "map/include/computeMap.hpp" namespace "skch" nogil:


    cdef cppclass Map:

        cppclass L1_candidateLocus_t:
            seqno_t  seqId
            offset_t rangeStartPos
            offset_t rangeEndPos

        cppclass L2_mapLocus_t:
            seqno_t         seqId
            offset_t        meanOptimalPos
            Sketch.MIIter_t optimalStart
            Sketch.MIIter_t optimalEnd
            int             sharedSketchSize

        ctypedef Sketch.MI_Type MinVec_Type
        ctypedef Sketch.MIIter_t MIIter_t

        vector[ContigInfo] metadata

        Map(
            const Parameters &p,
            const Sketch &refsketch,
            uint64_t &totalQueryFragments,
            int queryno,
            function[void(MappingResult&)] f = nullptr
        )

        Map(
            const Parameters &p,
            const Sketch &refsketch,
            uint64_t &totalQueryFragments,
            int queryno,
            MappingResultsVector_t &r
        )