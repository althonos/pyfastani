
cdef extern from "common/kseq.h":

    ctypedef __kstring_t kstring_t
    cdef struct __kstring_t:
        size_t l
        size_t m
        char*  s

    cdef cppclass kseq_t:
        kstring_t seq
