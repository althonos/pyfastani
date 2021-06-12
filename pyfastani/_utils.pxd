from libc.stdint cimport uint64_t

from fastani.map.base_types cimport MappingResultsVector_t
from fastani.map.compute_map cimport Map
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch


cdef extern from "_utils.hpp" nogil:

    Map* new_map_with_result_vector(
        const Parameters &p,
        const Sketch &refsketch,
        uint64_t &totalQueryFragments,
        int queryno,
        MappingResultsVector_t &r
    )
