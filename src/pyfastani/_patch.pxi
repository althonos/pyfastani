cdef extern from *:
    """
    #ifndef HAVE_PYINTERPRETERSTATE_GETID
    int64_t PyInterpreterState_GetID(PyInterpreterState *interp) {
        return 0;
    }
    #endif
    """