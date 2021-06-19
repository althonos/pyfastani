#ifndef __SEQUTILS_H
#define __SEQUTILS_H

#include "complement.h"

#ifdef __cplusplus
extern "C" {
#endif

void default_copy_upper(char*, const char*, size_t);
#ifdef SSE2_BUILD_SUPPORTED
void sse2_copy_upper(char*, const char*, size_t);
#endif
#ifdef NEON_BUILD_SUPPORTED
void neon_copy_upper(char* dst, const char* src, size_t len);
#endif
void copy_upper(char*, const char*, size_t);

void default_reverse_complement(char*, const char*, size_t);
#ifdef SSSE3_BUILD_SUPPORTED
extern void ssse3_reverse_complement(char* dst, const char* src, size_t len);
#endif
void reverse_complement(char*, const char*, size_t);

#ifdef __cplusplus
}
#endif

#endif
