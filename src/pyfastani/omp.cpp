#include "omp.h"

// not needed anywhere except in `cgi::correctRefGenomeIds` so we can just
// patch these functions to disable logging
int omp_get_thread_num(void) {
    return 1; // Make the logger shut up.
}

int omp_get_num_threads(void) {
    return 1;
}
