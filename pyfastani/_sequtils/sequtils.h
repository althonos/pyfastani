#ifndef __SEQUTILS_H
#define __SEQUTILS_H

#ifdef __X86__
#include "cpuinfo_x86.h"
#endif
#ifdef __X86_64__
#include "cpuinfo_x86.h"
#endif
#ifdef __arm__
#include "cpuinfo_arm.h"
#endif
#ifdef __aarch64__
#include "cpuinfo_aarch64.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void default_copy_upper(char*, const char*, size_t);
#ifdef SSE2_BUILD_SUPPORTED
extern void sse2_copy_upper(char*, const char*, size_t);
#endif
#ifdef NEON_BUILD_SUPPORTED
extern void neon_copy_upper(char* dst, const char* src, size_t len);
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
