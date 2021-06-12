from libc.stddef cimport wchar_t


cdef extern from "<ostream>" namespace "std" nogil:

    cdef cppclass basic_ostream[CharT]:
        pass

    ctypedef basic_ostream[char]    ostream
    ctypedef basic_ostream[wchar_t] wostream
