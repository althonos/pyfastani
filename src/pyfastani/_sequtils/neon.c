#include <ctype.h>
#include <x86intrin.h>

extern void neon_copy_upper(char* dst, const char* src, size_t len) {
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

    while (len-- > 0) {
        *dst = toupper(*src);
        src++;
        dst++;
    }
}
