cdef extern from "_simd.h" nogil:
    void copy_upper(char*, const char*, size_t)
