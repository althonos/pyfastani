#include <ctype.h>
#include <x86intrin.h>

extern void sse2_copy_upper(char* dst, const char* src, size_t len) {
    const __m128i offset  = _mm_set1_epi8('a' - 'A');

    while (len >= sizeof(__m128i)) {
        __m128i inp = _mm_loadu_si128((__m128i*) src);
        __m128i added = _mm_andnot_si128(offset, inp);
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
