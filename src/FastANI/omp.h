/* A stub `omp.h` to get compilation to work on OSX */
#ifndef _OMP_H
#define _OMP_H

#ifdef __cplusplus
extern "C" {
#endif

extern int omp_get_thread_num(void);
extern int omp_get_num_threads(void);
extern void omp_set_num_threads(int);

#ifdef __cplusplus
}
#endif
#endif
