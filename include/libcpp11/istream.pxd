from libc.stddef cimport wchar_t


cdef extern from "<istream>" namespace "std" nogil:

    cdef cppclass basic_istream[CharT]:
        pass

    ctypedef basic_istream[char]    istream
    ctypedef basic_istream[wchar_t] wistream

    cdef cppclass basic_iostream[CharT]:
        pass

    ctypedef basic_iostream[char]    iostream
    ctypedef basic_iostream[wchar_t] wiostream
