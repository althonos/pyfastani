#ifndef __SIMD_H
#define __SIMD_H

#ifdef __cplusplus
extern "C" {
#endif

void default_copy_upper(char*, const char*, size_t);
#ifdef __SSE2__
void sse_copy_upper(char*, const char*, size_t);
#endif
#ifdef __ARM_NEON__
void neon_copy_upper(char* dst, const char* src, size_t len);
#endif
void copy_upper(char*, const char*, size_t);

#ifdef __cplusplus
}
#endif

#endif
