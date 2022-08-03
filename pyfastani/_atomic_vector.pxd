from libcpp.vector cimport vector

cdef extern from "_atomic_vector.hpp" nogil:
    cdef cppclass atomic_vector[T](vector[T]):
        pass
