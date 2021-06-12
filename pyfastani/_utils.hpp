#ifndef __UTILS_HPP
#define __UTILS_HPP

#include "map/include/computeMap.hpp"

#include <stdint.h>
#include <functional>

#include "map/include/base_types.hpp"
#include "map/include/computeMap.hpp"

skch::Map* new_map_with_result_vector(
    const skch::Parameters &p,
    const skch::Sketch &refsketch,
    uint64_t &totalQueryFragments,
    int queryno,
    skch::MappingResultsVector_t &r
);

// not needed anywhere except in `cgi::correctRefGenomeIds`
// so we can just patch them
extern int omp_get_thread_num(void);
extern int omp_get_num_threads(void);

#endif
