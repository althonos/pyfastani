from libc.stdint cimport uint64_t


cdef extern from "map/include/map_stats.hpp" namespace "skch::Stat" nogil:

    cdef float  j2md(float j, float k)
    cdef float  md2j(float d, int k)
    cdef float  md_lower_bound(float d, int s, int k, float ci)

    cdef int estimateMinimumHits(int s, int k, float perc_identity)
    cdef int estimateMinimumHitsRelaxed(int s, int k, float perc_identity)

    cdef double estimate_pvalue(
        int s,
        int k,
        int alphabetSize,
        float identity,
        int lengthQuery,
        uint64_t lengthReference
    )

    cdef int recommendedWindowSize(
        double pValue_cutoff,
        int k,
        int alphabetSize,
        float identity,
        int lengthQuery,
        uint64_t lengthReference
    )
