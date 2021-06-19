#include <ctype.h>
#include <stddef.h>

#include "sequtils.h"
#include "complement.h"

#ifdef __X86__ || __X86_64__
#include "cpu_features_x86.h"
static const X86Features features = GetX86Info().features;
#endif
#ifdef __arm__
#include "cpu_features_arm.h"
static const ArmFeatures features = GetArmInfo().features;
#endif
#ifdef __aarch64__
#include "cpu_features_aarch64.h"
static const Aarch64Features features = GetAarch64Info().features;
#endif

// --- Fast copy with uppercase ----------------------------------------------

void default_copy_upper(char* dst, const char* src, size_t len) {
    while (len > 16) {
        for (size_t i = 0; i < 16; ++i) {
            *dst = toupper(*src);
            src++;
            dst++;
        }
        len -= 16;
    }
    while (len-- > 0) {
        *dst = toupper(*src);
        src++;
        dst++;
    }
}

void copy_upper(char* dst, const char* src, size_t len) {
    #ifdef __arm__
      if (features.neon)
        return neon_copy_upper(dst, src, len);
      else
    #endif
    #if defined(__X86__) || defined(__X86_64__)
      if (features.sse2)
        return sse2_copy_upper(dst, src, len); // fast copying plus upper.
      else
    #endif
        return default_copy_upper(dst, src, len);
}


// --- Fast reverse complement -----------------------------------------------

void default_reverse_complement(char* dst, const char* src, size_t len) {
    while (len > 16) {
        for (size_t i = 0; i < 16; ++i) {
            *dst = complement(src[--len]);
            dst++;
        }
    }
    while (len > 0) {
        *dst = complement(src[--len]);
        dst++;
    }
}

void reverse_complement(char* dst, const char* src, size_t len) {
    #if defined(__X86__) || defined(__X86_64__)
      if (features.ssse3)
        return ssse3_reverse_complement(dst, src, len); // fast reverse complement.
      else
    #endif
        return default_reverse_complement(dst, src, len);
}
