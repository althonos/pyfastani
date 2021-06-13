#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <stdint.h>
#include <chrono>
#include <limits>
#include <iostream>
#include <vector>
#include <zlib.h>

#include "common/kseq.h"
#include "map/include/base_types.hpp"
#include "map/include/winSketch.hpp"

// not needed anywhere except in `cgi::correctRefGenomeIds`
// so we can just patch them
extern int omp_get_thread_num(void);
extern int omp_get_num_threads(void);

typedef kseq_t* kseq_ptr_t;

#endif
