#include <x86intrin.h>

#include "complement.h"

// Adapted from: https://github.com/adamkewley/reverse-complement/
extern void ssse3_reverse_complement(char* dst, const char* src, size_t len) {
    const __m128i shuffle_mask = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    const __m128i and_mask = _mm_set1_epi8(0x1f);
    const __m128i sixteen = _mm_set1_epi8(16);
    const __m128i lo_lut = _mm_set_epi8(
        '\0', 'N',  'K', '\0',
        'M',  '\n', '\0', 'D',
        'C',  '\0', '\0', 'H',
        'G',  'V',  'T',  '\0'
    );
    const __m128i hi_lut = _mm_set_epi8(
        '\0', '\0', '\0', '\0',
        '\0', '\0', 'R',  '\0',
        'W',  'B',  'A',  'A',
        'S',  'Y',  '\0', '\0'
    );

    while (len > sizeof(__m128i)) {
        //
        len -= sizeof(__m128i);
        __m128i v = _mm_loadu_si128((__m128i*) &src[len]);

        // reverse elements in the registers
         v =  _mm_shuffle_epi8(v, shuffle_mask);

         // AND all elements with 0x1f, so that a smaller LUT (< 32 bytes)
         // can be used. This is important with SIMD because, unlike
         // single-char complement (above), SIMD uses 16-byte shuffles. The
         // single-char LUT would require four shuffles, this LUT requires
         // two.
         v = _mm_and_si128(v, and_mask);

         // Lookup for all elements <16
         __m128i lo_mask = _mm_cmplt_epi8(v, sixteen);
         __m128i lo_els = _mm_and_si128(v, lo_mask);
         __m128i lo_vals = _mm_shuffle_epi8(lo_lut, lo_els);

         // Lookup for all elements >16
         __m128i hi_els = _mm_sub_epi8(v, sixteen);
         __m128i hi_vals = _mm_shuffle_epi8(hi_lut, hi_els);

         // OR both lookup results
         _mm_storeu_si128((__m128i*) dst, _mm_or_si128(lo_vals, hi_vals));
         dst += sizeof(__m128i);
    }

    while (len > 0) {
        *dst = complement(src[--len]);
        dst++;
    }
}
