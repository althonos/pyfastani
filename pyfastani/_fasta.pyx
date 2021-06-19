# coding: utf-8
# cython: language_level=3, linetrace=True, language=cpp, binding=False

# --- C imports --------------------------------------------------------------

from cpython cimport PyObject
from cpython.bytes cimport PyBytes_FromStringAndSize
from libc.stdlib cimport malloc, realloc, free
from libc.stdio cimport FILE, fopen, fgets, fclose
from libc.errno cimport errno
from libc.string cimport strlen, memcpy, memset

from _unicode cimport PyUnicode_1BYTE_KIND, PyUnicode_FromKindAndData
from _sequtils cimport copy_upper

cdef extern from "<ctype.h>" nogil:
    cdef int toupper(int c)


# --- Python imports ---------------------------------------------------------

import os


# --- Cython classes ---------------------------------------------------------

cdef class Record:
    cdef public str   id
    cdef public bytes seq

    def __init__(self, str id, bytes seq):
        self.id = id
        self.seq = seq


cdef class Parser:
    cdef str    path
    cdef char   line[2048]
    cdef FILE*  file

    cdef char*  seq_buffer
    cdef size_t buffer_length


    def __cinit__(self, str path):
        self.path = path
        self.file = fopen(os.fsencode(path), "r")
        if self.file == NULL:
            raise OSError(errno, path)
        fgets(<char*> &self.line, sizeof(self.line), self.file)

        self.buffer_length = 0x8000
        self.seq_buffer = <char*> malloc(self.buffer_length * sizeof(char))

    def __dealloc__(self):
        fclose(self.file)
        free(self.seq_buffer)

    def __iter__(self):
        return self

    def __next__(self):

        cdef size_t seq_length  = 0
        cdef size_t line_length = 0
        cdef size_t line_pos    = 0
        cdef str    id
        cdef bytes  seq

        if self.line[0] != b">":
            raise StopIteration()

        line_length = strlen(<char*> &self.line)
        if self.line[line_length - 1] != b"\n":
            raise BufferError("FASTA identifier too large for the line buffer")
        id = PyUnicode_FromKindAndData(PyUnicode_1BYTE_KIND, <void*> &self.line[1], line_length - 2)

        with nogil:
            while fgets(<char*> &self.line, sizeof(self.line), self.file) != NULL:
                if self.line[0] == b'>':
                    break

                line_length = strlen(<char*> &self.line)
                if line_length + seq_length >= self.buffer_length:
                    self.buffer_length *= 2
                    self.seq_buffer = <char*> realloc(<void*> self.seq_buffer, self.buffer_length)

                line_pos = 0
                if self.line[line_length - 1] == b"\n":
                    line_length -= 1

                copy_upper(&self.seq_buffer[seq_length], <char*> &self.line, line_length)
                seq_length += line_length

        seq = PyBytes_FromStringAndSize(self.seq_buffer, seq_length)
        return Record(id, seq)
