#include <ctype.h>
#include <stddef.h>

#include "sequtils.h"
#include "complement.h"


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

void (*resolve_copy_upper (void))(char*, const char*, size_t) {
  // ifunc resolvers fire before constructors, explicitly call the init
  // function.
#ifdef SSE2_BUILD_SUPPORTED
  __builtin_cpu_init ();
  if (__builtin_cpu_supports ("sse2"))
    return sse2_copy_upper; // fast copying plus upper.
  else
#endif
#ifdef NEON_BUILD_SUPPORTED
__builtin_cpu_init ();
if (__builtin_cpu_supports ("neon"))
  return neon_copy_upper; // fast copying plus upper.
else
#endif
    return default_copy_upper;
}

void copy_upper(char*, const char*, size_t)
     __attribute__ ((ifunc ("resolve_copy_upper")));


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

void (*resolve_reverse_complement (void))(char*, const char*, size_t) {
#ifdef SSSE3_BUILD_SUPPORTED
    __builtin_cpu_init ();
    if (__builtin_cpu_supports ("ssse3"))
        return ssse3_reverse_complement; // fast reverse_complement.
    else
#endif
        return default_reverse_complement;
}

void reverse_complement(char*, const char*, size_t)
     __attribute__ ((ifunc ("resolve_reverse_complement")));
