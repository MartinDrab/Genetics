
#ifndef KTHREAD_H
#define KTHREAD_H

#ifndef _MSC_VER
#include <pthread.h>
#endif
#include "utils.h"


#ifdef __cplusplus
extern "C" {
#endif


#ifdef _MSC_VER
#define pthread_t	HANDLE
#endif

struct kt_for_t;

typedef struct {
	/** Reference to a structure describing the for loop. */
	struct kt_for_t *t;
	pthread_t ThreadHandle;
	volatile long Terminate;
	volatile long i;
} ktf_worker_t;


typedef void(kt_for_worker_routine)(void* Data, long WOrkIndex, size_t ThreadNo);

/** Parameters for one parallel for loop. */
typedef struct kt_for_t {
	int Allocated;
	/** Number of threads parallelizing the loop. */
	int n_threads;
	/** Number of iterations. */
	long n;
	/** worker structures, one for each thread. */
	ktf_worker_t *w;
	/** Routine representing one iteration. */
	kt_for_worker_routine *func;
	/** Context information passed to the worker routine. */
	void *data;
} kt_for_t;

void kt_for_init(kt_for_t *ForLoop, int n_threads, kt_for_worker_routine *func, void *data, long n);
int kt_for_prepare(kt_for_t *ForLoop);
void kt_for_perform(kt_for_t *ForLoop);
void kt_for_finit(kt_for_t *ForLoop, int Terminate);

void kt_for(int n_threads, kt_for_worker_routine *func, void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

void *kt_forpool_init(int n_threads);
void kt_forpool_destroy(void *_fp);
void kt_forpool(void *_fp, void (*func)(void*,long,int), void *data, long n);

#ifdef __cplusplus
}
#endif

#endif
