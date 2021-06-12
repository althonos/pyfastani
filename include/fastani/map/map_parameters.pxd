from libc.stdint cimport uint64_t
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef extern from "map/include/map_parameters.hpp" namespace "skch" nogil:

    cdef cppclass Parameters:
        int            kmerSize
        int            windowSize
        int            minReadLength
        float          minFraction
        int            threads
        int            alphabetSize
        uint64_t       referenceSize
        float          percentageIdentity
        double         p_value
        vector[string] refSequences
        vector[string] querySequences
        string         outFileName
        bool           reportAll
        bool           visualize
        bool           matrixOutput
