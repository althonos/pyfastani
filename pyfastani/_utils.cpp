#include "_utils.hpp"

int omp_get_thread_num(void) {
    return 1; // Make the logger shut up.
}

int omp_get_num_threads(void) {
    return 1;
}
