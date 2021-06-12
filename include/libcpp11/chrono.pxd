

cdef extern from "<chrono>" namespace "std::chrono" nogil:

    cdef cppclass high_resolution_clock:
        pass
