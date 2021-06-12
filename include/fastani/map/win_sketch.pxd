cimport libcpp11.iostream
from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map
from libcpp.vector cimport vector

from fastani.map.base_types cimport (
    offset_t,
    seqno_t,
    MinimizerMapKeyType,
    MinimizerMapValueType,
    MinimizerInfo,
    ContigInfo
)
from fastani.map.map_parameters cimport Parameters


cdef extern from "map/include/winSketch.hpp" namespace "skch" nogil:

    IF FASTANI_PRIVATE_ACCESS:

        cdef cppclass Sketch:

            ctypedef vector[MinimizerInfo].const_iterator MIIter_t
            ctypedef unordered_map[MinimizerMapKeyType, MinimizerMapValueType] MI_Map_t
            ctypedef vector[MinimizerInfo] MI_Type

            float              percentageThreshold # private
            int                freqThreshold # private
            Parameters&        param # private
            MI_Type            minimizerIndex # private
            map[int, int]      minimizerFreqHistogram # private
            vector[ContigInfo] metadata
            vector[seqno_t]    sequencesByFileInfo
            MI_Map_t           minimizerPosLookupIndex

            Sketch() # private
            Sketch(const Parameters &p)

            void build() # private
            void index() # private
            void computeFreqHist() # private
            int getFreqThreshold()
            MIIter_t searchIndex(seqno_t seqId, offset_t winpos)
            MIIter_t getMinimizerIndexEnd()

    ELSE:

        cdef cppclass Sketch:

            ctypedef vector[MinimizerInfo].const_iterator MIIter_t
            ctypedef unordered_map[MinimizerMapKeyType, MinimizerMapValueType] MI_Map_t
            ctypedef vector[MinimizerInfo] MI_Type

            vector[ContigInfo] metadata
            vector[seqno_t]    sequencesByFileInfo
            MI_Map_t           minimizerPosLookupIndex

            Sketch(const Parameters &p)

            int getFreqThreshold()
            MIIter_t searchIndex(seqno_t seqId, offset_t winpos)
            MIIter_t getMinimizerIndexEnd()
