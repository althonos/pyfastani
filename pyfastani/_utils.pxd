from libc.stdint cimport uint64_t

from kseq cimport kseq_t
from fastani.map.base_types cimport MappingResultsVector_t
from fastani.map.compute_map cimport Map
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch


cdef extern from "_utils.hpp" nogil:

    int omp_get_thread_num()
    int omp_get_num_threads()

    ctypedef kseq_t* kseq_ptr_t
