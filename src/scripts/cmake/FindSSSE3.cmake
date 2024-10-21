#.rst:
# FindSSSE3
# --------
#
# Finds SSSE3 support
#
# This module can be used to detect SSSE3 support in a C compiler.  If
# the compiler supports SSSE3, the flags required to compile with
# SSSE3 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support SSSE3.
#
# The following variables are set:
#
# ::
#
#    SSSE3_C_FLAGS - flags to add to the C compiler for SSSE3 support
#    SSSE3_FOUND - true if SSSE3 is detected
#
#=============================================================================

set(_SSSE3_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${SSSE3_FIND_QUIETLY})

# sample SSSE3 source code to test
set(SSSE3_C_TEST_SOURCE
"
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <smmintrin.h>
#endif
int foo() {
    __m128  a = _mm_set1_ps(1.0);
            a = _mm_blend_ps(a, _mm_setzero_ps(), 0b1111);
    float   x = _mm_extract_ps(a, 1);
    return (x == 1) ? 0 : 1;
}
int main(void) { return foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED SSSE3_C_FLAGS) OR (DEFINED HAVE_SSSE3))
else()
  if(WIN32)
    set(SSSE3_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSSE3
      " "
      "/arch:SSSE3")
  else()
    set(SSSE3_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts SSSE3
      " "
      #clang
      "-mssse3"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS SSSE3_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_SSSE3 CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try SSSE3 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${SSSE3_C_TEST_SOURCE}" HAVE_SSSE3)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_SSSE3)
      set(SSSE3_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(SSSE3_C_FLAG_CANDIDATES)

  set(SSSE3_C_FLAGS "${SSSE3_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for SSSE3 intrinsics")
endif()

list(APPEND _SSSE3_REQUIRED_VARS SSSE3_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_SSSE3_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(SSSE3
                                    REQUIRED_VARS ${_SSSE3_REQUIRED_VARS})

  mark_as_advanced(${_SSSE3_REQUIRED_VARS})

  unset(_SSSE3_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindSSSE3 requires C or CXX language to be enabled")
endif()
