#include <ctype.h>
#include <stddef.h>

#include "sequtils.h"
#include "complement.h"

#if defined(__x86__) || defined(__x86_64__)
#include "cpuinfo_x86.h"
using namespace cpu_features;
static const X86Features features = GetX86Info().features;
#endif
#ifdef __arm__
#include "cpuinfo_arm.h"
using namespace cpu_features;
static const ArmFeatures features = GetArmInfo().features;
#endif
#ifdef __aarch64__
#include "cpuinfo_aarch64.h"
using namespace cpu_features;
static const Aarch64Features features = GetAarch64Info().features;
#endif

extern "C" {
    // --- Fast copy with uppercase --------------------------------------------

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
        #ifdef NEON_BUILD_SUPPORTED
          if (features.neon)
            return neon_copy_upper(dst, src, len);
          else
        #endif
        #endif
        #ifdef __aarch64__
        #ifdef NEON_BUILD_SUPPORTED
          if (features.neon)
            return neon_copy_upper(dst, src, len);
          else
        #endif
        #endif
        #if defined(__x86__) || defined(__x86_64__)
        #ifdef SSE2_BUILD_SUPPORTED
          if (features.sse2)
            return sse2_copy_upper(dst, src, len); // fast copying plus upper.
          else
        #endif
        #endif
            return default_copy_upper(dst, src, len);
    }


    // --- Fast reverse complement ---------------------------------------------

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
        #if defined(__x86__) || defined(__x86_64__)
        #ifdef SSSE3_BUILD_SUPPORTED
          if (features.ssse3)
            return ssse3_reverse_complement(dst, src, len); // fast reverse complement.
          else
        #endif
        #endif
            return default_reverse_complement(dst, src, len);
    }

}
