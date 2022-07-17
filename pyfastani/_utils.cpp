#include <zlib.h>
#include "_utils.hpp"

ZEXTERN gzFile ZEXPORT gzdopen(int fd, const char* mode) {
    return NULL;
}

ZEXTERN gzFile ZEXPORT gzopen(const char* path, const char* mode) {
    return NULL;
}

ZEXTERN gzFile ZEXPORT gzopen64(const char* path, const char* mode) {
    return NULL;
}

ZEXTERN int ZEXPORT gzread(gzFile file, void* buf, unsigned int len) {
    return 0;
}

ZEXTERN int ZEXPORT gzclose(gzFile file) {
    return 0;
}
