#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>

#include "_utils.hpp"


skch::Map* new_map_with_result_vector(
    const skch::Parameters &p,
    const skch::Sketch &refsketch,
    uint64_t &totalQueryFragments,
    int queryno,
    skch::MappingResultsVector_t &r
) {
  auto fn = std::bind(skch::Map::insertL2ResultsToVec, std::ref(r), std::placeholders::_1);
  return new skch::Map( p, refsketch, totalQueryFragments, queryno, fn );
}

int omp_get_thread_num(void) {
    return 0;
}

int omp_get_num_threads(void) {
    return 1;
}
