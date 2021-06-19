#include <ctype.h>
#include <x86intrin.h>

extern void sse2_copy_upper(char* dst, const char* src, size_t len) {
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

    while (len-- > 0) {
        *dst = toupper(*src);
        src++;
        dst++;
    }
}
