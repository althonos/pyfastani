#include "_utils.hpp"

int omp_get_thread_num(void) {
    return 1; // Make the logger shut up.
}

int omp_get_num_threads(void) {
    return 1;
}

ZEXTERN gzFile ZEXPORT gzdopen(int fd, const char* mode) {
    return NULL;
}

ZEXTERN gzFile ZEXPORT gzopen64(const char* path, const char* mode) {
    return NULL;
}

ZEXTERN int ZEXPORT gzread(gzFile file, void* buf, unsigned int len) {
    return 0;
}

ZEXTERN int ZEXPORT gzclose(gzFile file) {}
