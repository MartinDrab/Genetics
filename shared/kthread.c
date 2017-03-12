
#ifdef _MSC_VER
#include <windows.h>
#else
#include <pthread.h>
#endif
#include <stdlib.h>
#include <limits.h>
#include "kthread.h"

/************
 * kt_for() *
 ************/

#ifdef _MSC_VER
#define alloca		_alloca
#endif


static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
#ifndef _MSC_VER
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
#else
	k = InterlockedAdd(&t->w[min_i].i, t->n_threads) - t->n_threads;
#endif
	return k >= t->n? -1 : k;
}

#ifndef _MSC_VER
static void *ktf_worker(void *data)
#else
static DWORD WINAPI ktf_worker(void *data)
#endif
{
	ktf_worker_t *w = (ktf_worker_t*)data;

	if (!w->Terminate) {
		long i;
		for (;;) {
#ifndef _MSC_VER
			i = __sync_fetch_and_add(&w->i, w->t->n_threads);
#else
			i = InterlockedAdd(&w->i, w->t->n_threads) - w->t->n_threads;
#endif
			if (i >= w->t->n) 
				break;
			
			w->t->func(w->t->data, i, w - w->t->w);
		}
		
		while ((i = steal_work(w->t)) >= 0)
			w->t->func(w->t->data, i, w - w->t->w);
	}

#ifndef _MSC_VER
	pthread_exit(0);
#else

	return 0;
#endif
}


static void _kt_for_init_worker_data(kt_for_t *ForLoop, int ThreadNo)
{
	ktf_worker_t *workerData = &ForLoop->w[ThreadNo];

	workerData->t = ForLoop;
	workerData->i = ThreadNo;
	workerData->ThreadHandle = (pthread_t)NULL;
	workerData->Terminate = 0;

	return;
}


/************************************************************************/
/*             PUBLIC FUNCTIONS                                         */
/************************************************************************/


void kt_for_init(kt_for_t *ForLoop, int n_threads, kt_for_worker_routine *func, void *data, long n)
{
	ForLoop->func = func;
	ForLoop->data = data; 
	ForLoop->n_threads = n_threads;
	ForLoop->n = n;
	ForLoop->w = NULL;
	ForLoop->Allocated = 0;

	return;
}


int kt_for_prepare(kt_for_t *ForLoop)
{
	int ret = 0;

	ForLoop->w = (ktf_worker_t*)calloc(ForLoop->n_threads, sizeof(ktf_worker_t));
	if (ForLoop->w != NULL) {
		for (int i = 0; i < ForLoop->n_threads; ++i) {
			_kt_for_init_worker_data(ForLoop, i);
#ifdef _MSC_VER
			DWORD threadId;

			ForLoop->w[i].ThreadHandle = CreateThread(NULL, 0, ktf_worker, &ForLoop->w[i], CREATE_SUSPENDED, &threadId);
			if (ForLoop->w[i].ThreadHandle == NULL) {
				for (int j = 0; j < i; ++j) {
					HANDLE th = ForLoop->w[j].ThreadHandle;

					ForLoop->w[j].Terminate = 1;
					ResumeThread(th);
					WaitForSingleObject(th, INFINITE);
					CloseHandle(th);
				}

				ret = -2;
				break;
			}
#endif
		}

		if (ret == 0)
			ForLoop->Allocated = 1;

		if (ret != 0)
			free(ForLoop->w);
	} else ret = -1;

	return ret;
}


void kt_for_finit(kt_for_t *ForLoop, int Terminate)
{
	ktf_worker_t *w = ForLoop->w;

	for (int i = 0; i < ForLoop->n_threads; ++i) {
		if (w->ThreadHandle != (pthread_t)NULL) {
			w->Terminate = Terminate;
#ifndef _MSC_VER
			pthread_join(w->ThreadHandle, 0);
#else
			ResumeThread(w->ThreadHandle);
			WaitForSingleObject(w->ThreadHandle, INFINITE);
			CloseHandle(w->ThreadHandle);
#endif
			w->ThreadHandle = (pthread_t)NULL;
		}

		++w;
	}

	if (ForLoop->Allocated) {
		free(ForLoop->w);
		ForLoop->w = NULL;
	}

	return;
}


void kt_for_perform(kt_for_t *ForLoop)
{
	ktf_worker_t *w = ForLoop->w;

	for (int i = 0; i < ForLoop->n_threads; ++i) {
#ifndef _MSC_VER
		pthread_create(&w->ThreadHandle, 0, ktf_worker, w);
#else
		if (w->ThreadHandle == NULL) {
			DWORD threadId;

			w->ThreadHandle = CreateThread(NULL, 0, ktf_worker, w, 0, &threadId);
		} else ResumeThread(w->ThreadHandle);
#endif

		++w;
	}

	return;
}


void kt_for(int n_threads, kt_for_worker_routine *func, void *data, long n)
{
	if (n_threads > 1) {
		kt_for_t t;

		kt_for_init(&t, n_threads, func, data, n);
		t.w = (ktf_worker_t*)alloca(n_threads * sizeof(ktf_worker_t));
		for (int i = 0; i < n_threads; ++i)
			_kt_for_init_worker_data(&t, i);

		kt_for_perform(&t);
		kt_for_finit(&t, 0);
	} else {
		for (long j = 0; j < n; ++j) 
			func(data, j, 0);
	}
}

#ifndef _MSC_VER

/***************************
 * kt_for with thread pool *
 ***************************/

struct kt_forpool_t;

typedef struct {
	struct kt_forpool_t *t;
	long i;
	int action;
} kto_worker_t;

typedef struct kt_forpool_t {
	int n_threads, n_pending;
	long n;
	pthread_t *tid;
	kto_worker_t *w;
	void (*func)(void*,long,int);
	void *data;
	pthread_mutex_t mutex;
	pthread_cond_t cv_m, cv_s;
} kt_forpool_t;

static inline long kt_fp_steal_work(kt_forpool_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void *kt_fp_worker(void *data)
{
	kto_worker_t *w = (kto_worker_t*)data;
	kt_forpool_t *fp = w->t;
	for (;;) {
		long i;
		int action;
		pthread_mutex_lock(&fp->mutex);
		if (--fp->n_pending == 0)
			pthread_cond_signal(&fp->cv_m);
		w->action = 0;
		while (w->action == 0) pthread_cond_wait(&fp->cv_s, &fp->mutex);
		action = w->action;
		pthread_mutex_unlock(&fp->mutex);
		if (action < 0) break;
		for (;;) { // process jobs allocated to this worker
			i = __sync_fetch_and_add(&w->i, fp->n_threads);
			if (i >= fp->n) break;
			fp->func(fp->data, i, w - fp->w);
		}
		while ((i = kt_fp_steal_work(fp)) >= 0) // steal jobs allocated to other workers
			fp->func(fp->data, i, w - fp->w);
	}
	pthread_exit(0);
}

void *kt_forpool_init(int n_threads)
{
	kt_forpool_t *fp;
	int i;
	fp = (kt_forpool_t*)calloc(1, sizeof(kt_forpool_t));
	fp->n_threads = fp->n_pending = n_threads;
	fp->tid = (pthread_t*)calloc(fp->n_threads, sizeof(pthread_t));
	fp->w = (kto_worker_t*)calloc(fp->n_threads, sizeof(kto_worker_t));
	for (i = 0; i < fp->n_threads; ++i) fp->w[i].t = fp;
	pthread_mutex_init(&fp->mutex, 0);
	pthread_cond_init(&fp->cv_m, 0);
	pthread_cond_init(&fp->cv_s, 0);
	for (i = 0; i < fp->n_threads; ++i) pthread_create(&fp->tid[i], 0, kt_fp_worker, &fp->w[i]);
	pthread_mutex_lock(&fp->mutex);
	while (fp->n_pending) pthread_cond_wait(&fp->cv_m, &fp->mutex);
	pthread_mutex_unlock(&fp->mutex);
	return fp;
}

void kt_forpool_destroy(void *_fp)
{
	kt_forpool_t *fp = (kt_forpool_t*)_fp;
	int i;
	for (i = 0; i < fp->n_threads; ++i) fp->w[i].action = -1;
	pthread_cond_broadcast(&fp->cv_s);
	for (i = 0; i < fp->n_threads; ++i) pthread_join(fp->tid[i], 0);
	pthread_cond_destroy(&fp->cv_s);
	pthread_cond_destroy(&fp->cv_m);
	pthread_mutex_destroy(&fp->mutex);
	free(fp->w); free(fp->tid); free(fp);
}

void kt_forpool(void *_fp, void (*func)(void*,long,int), void *data, long n)
{
	kt_forpool_t *fp = (kt_forpool_t*)_fp;
	long i;
	if (fp && fp->n_threads > 1) {
		fp->n = n, fp->func = func, fp->data = data, fp->n_pending = fp->n_threads;
		for (i = 0; i < fp->n_threads; ++i) fp->w[i].i = i, fp->w[i].action = 1;
		pthread_mutex_lock(&fp->mutex);
		pthread_cond_broadcast(&fp->cv_s);
		while (fp->n_pending) pthread_cond_wait(&fp->cv_m, &fp->mutex);
		pthread_mutex_unlock(&fp->mutex);
	} else for (i = 0; i < n; ++i) func(data, i, 0);
}

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

static void *ktp_worker(void *data)
{
	ktp_worker_t *w = (ktp_worker_t*)data;
	ktp_t *p = w->pl;
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		w->data = p->func(p->shared, w->step, w->step? w->data : 0); // for the first step, input is NULL

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;
		if (w->step == 0) w->index = p->index++;
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	pthread_t *tid;
	int i;

	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	pthread_mutex_init(&aux.mutex, 0);
	pthread_cond_init(&aux.cv, 0);

	aux.workers = (ktp_worker_t*)alloca(n_threads * sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

	pthread_mutex_destroy(&aux.mutex);
	pthread_cond_destroy(&aux.cv);
}

#endif