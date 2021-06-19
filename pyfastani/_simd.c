#include <ctype.h>
#include <stddef.h>
#ifdef __SSE2__
#include <x86intrin.h>
#endif
#ifdef __ARM_NEON__
#include <arm_neon.h>
#endif

#include "_simd.h"

void default_copy_upper(char* dst, const char* src, size_t len) {
    while (len-- > 0) {
        *dst = toupper(*src);
        src++;
        dst++;
    }
}

#ifdef __SSE2__
void sse2_copy_upper(char* dst, const char* src, size_t len) {

    const __m128i ascii_a = _mm_set1_epi8('a' - 1);
    const __m128i ascii_z = _mm_set1_epi8('z');
    const __m128i offset  = _mm_set1_epi8('a' - 'A');

    while (len >= sizeof(__m128i)) {
        __m128i inp = _mm_loadu_si128((__m128i*) src);
        __m128i greater_than_a = _mm_cmpgt_epi8(inp, ascii_a);
        __m128i less_equal_z = _mm_cmpgt_epi8(ascii_z, inp);
        __m128i mask = _mm_and_si128(greater_than_a, less_equal_z);
        __m128i diff = _mm_and_si128(mask, offset);
        __m128i added = _mm_sub_epi8(inp, diff);
        _mm_storeu_si128((__m128i *) dst, added);
        len -= sizeof(__m128i);
        src += sizeof(__m128i);
        dst += sizeof(__m128i);
    }

    default_copy_upper(dst, src, len);
}
#endif
#ifdef __ARM_NEON__
void neon_copy_upper(char* dst, const char* src, size_t len) {

    const int8x16_t ascii_a = vdupq_n_s8('a' - 1);
    const int8x16_t ascii_z = vdupq_n_s8('z');
    const int8x16_t offset  = vdupq_n_s8('a' - 'A');

    while (len >= sizeof(int8x16_t)) {
        int8x16_t inp = vld1q_u8((int8_t*) src);
        int8x16_t greater_than_a = vcgtq_s8(inp, ascii_a);
        int8x16_t less_equal_z = vcgtq_s8(ascii_z, inp);
        int8x16_t mask = vandq_s8(greater_than_a, less_equal_z);
        int8x16_t diff = vandq_s8(mask, offset);
        int8x16_t added = vsubq_s8(inp, diff);
        vst1q_s8((int8_t*) dst, added);
        len -= sizeof(int8x16_t);
        src += sizeof(int8x16_t);
        dst += sizeof(int8x16_t);
    }

    default_copy_upper(dst, src, len);
}
#endif

void (*resolve_copy_upper (void))(char*, const char*, size_t)
{
  // ifunc resolvers fire before constructors, explicitly call the init
  // function.
#ifdef __SSE2__
  __builtin_cpu_init ();
  if (__builtin_cpu_supports ("sse2"))
    return sse2_copy_upper; // fast copying plus upper.
  else
#endif
#ifdef __ARM_NEON__
__builtin_cpu_init ();
if (__builtin_cpu_supports ("neon"))
  return neon_copy_upper; // fast copying plus upper.
else
#endif
    return default_copy_upper;
}

void copy_upper(char*, const char*, size_t)
     __attribute__ ((ifunc ("resolve_copy_upper")));
