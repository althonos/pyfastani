cdef extern from "sequtils.h" nogil:
    void copy_upper(char*, const char*, size_t)
    void reverse_complement(char*, const char*, size_t)
