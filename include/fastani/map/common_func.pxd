from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp.vector cimport vector

from fastani.map.base_types cimport hash_t, seqno_t

cdef extern from "map/include/commonFunc.hpp" namespace "skch::CommonFunc" nogil:

    cdef int seed

    void reverseComplement(const char * src, char * dest, int length)
    hash_t getHash(const char * seq, int length)

    void addMinimizers[T, KSEQ](
        vector[T] &minimizerIndex,
        KSEQ kseq,
        int kmerSize,
        int windowSize,
        int alphabetSize,
        seqno_t seqCounter
    )

    uint64_t getReferenceSize(const vector[string] &refSequences)

    void ltrim(string &s)
    void rtrim(string &s)
    void trim(string &s)
