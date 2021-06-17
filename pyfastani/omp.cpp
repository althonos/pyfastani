#include "omp.h"

// not needed anywhere except in `cgi::correctRefGenomeIds` so we can just
// patch these functions to disable logging
extern int omp_get_thread_num(void);
extern int omp_get_num_threads(void);
