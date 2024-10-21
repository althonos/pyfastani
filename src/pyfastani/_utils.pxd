from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

from kseq cimport kseq_t
from fastani.map.compute_map cimport Map
from fastani.map.map_parameters cimport Parameters
from fastani.map.win_sketch cimport Sketch
from fastani.map.base_types cimport (
    MappingResultsVector_t,
    MinimizerInfo as MinimizerInfo_t
)


cdef extern from *:
    cppclass auto:
        pass


cdef extern from "<ctype.h>" nogil:
    int toupper(int)


cdef extern from "<iterator>" namespace "std" nogil:
    cdef ssize_t distance[I](I first, I last);


# cdef extern from "<zlib.h>" nogil:
#     cdef struct gzFile_s:
#         pass
#
#     ctypedef gzFile_s* gzFile
#
#     gzFile gzopen(int fd, const char* mode)
#     gzFile gzopen64(const char* path, const char* mode)
#     int gzread(gzFile file, void* buf, unsigned int len)
#     int gzclose(gzFile file)


cdef extern from "_utils.hpp" nogil:

    ctypedef kseq_t* kseq_ptr_t
