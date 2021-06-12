from libcpp11.istream cimport istream
from libcpp11.ostream cimport ostream


cdef extern from "<iostream>" namespace "std" nogil:

    cdef istream cin
    cdef ostream cout
    cdef ostream cerr



#   template<class CharT, class Traits = char_traits<CharT>>
#   class basic_istream;
#
# using istream  = basic_istream<char>;
# using wistream = basic_istream<wchar_t>;
#
# template<class CharT, class Traits = char_traits<CharT>>
#   class basic_iostream;
#
# using iostream  = basic_iostream<char>;
# using wiostream = basic_iostream<wchar_t>;
#
# template<class CharT, class Traits>
#   basic_istream<CharT, Traits>& ws(basic_istream<CharT, Traits>& is);
#
# template<class Istream, class T>
#   Istream&& operator>>(Istream&& is, T&& x);
