from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.functional cimport function
from libcpp11.fstream cimport ofstream

from fastani.map.base_types cimport seqno_t, offset_t, ContigInfo, MappingResult, MappingResultsVector_t
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch


cdef extern from "map/include/computeMap.hpp" namespace "skch" nogil:

    cdef cppclass Map:

        Parameters &param
        Sketch &refSketch

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

        void mapSingleQuerySeq[Q](Q&, MappingResultsVector_t&, ofstream&) except +
        void doL1Mapping[Q, V](Q&, V&)
        void computeL1CandidateRegions[Q, V1, V2](Q&, V1&, int, V2&) except +
        void doL2Mapping[Q, V1, V2](Q&, V1&, V2&) except +
        void computeL2MappedRegions[Q](Q&, L1_candidateLocus_t&, L2_mapLocus_t&) except +
        void reportL2Mappings(MappingResultsVector_t&, ofstream&) except +


cdef extern from "map/include/computeMap.hpp" namespace "skch::Map" nogil:

    cdef struct L1_candidateLocus_t:
        seqno_t  seqId
        offset_t rangeStartPos
        offset_t rangeEndPos

    cdef struct L2_mapLocus_t:
        seqno_t         seqId
        offset_t        meanOptimalPos
        Sketch.MIIter_t optimalStart
        Sketch.MIIter_t optimalEnd
        int             sharedSketchSize
