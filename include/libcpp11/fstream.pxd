from libc.stddef cimport wchar_t


cdef extern from "<ostream>" namespace "std" nogil:

    cdef cppclass basic_ifstream[CharT]:
        basic_ifstream()

    cdef cppclass basic_ofstream[CharT]:
        basic_ofstream()

    cdef cppclass basic_fstream[CharT]:
        basic_fstream()


    ctypedef basic_ifstream[char]    ifstream
    ctypedef basic_ifstream[wchar_t] wifstream

    ctypedef basic_ofstream[char]    ofstream
    ctypedef basic_ofstream[wchar_t] wofstream

    ctypedef basic_fstream[char]    fstream
    ctypedef basic_fstream[wchar_t] wfstream
