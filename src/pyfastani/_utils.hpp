#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <stdint.h>
#include <algorithm>
#include <chrono>
#include <limits>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "common/kseq.h"
#include "map/include/base_types.hpp"
#include "map/include/winSketch.hpp"

extern "C" {
    // compatibility layer for Cython
    typedef kseq_t* kseq_ptr_t;
}

#endif
