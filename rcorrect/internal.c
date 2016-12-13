
#include <omp.h>
#include "fml.h"
#include "internal.h"




void kt_for(int n_threads, void(*func)(void*, long, int), void *data, long n)
{
	int i = 0;

	omp_set_num_threads(n_threads);
#pragma omp parallel for shared(data)
	for (i = 0; i < n; ++i)
		func(data, i, omp_get_thread_num());

	return;
}


void fml_opt_init(fml_opt_t *opt)
{
	opt->n_threads = 1;
	opt->ec_k = 0;
	opt->min_cnt = 4;
	opt->max_cnt = 8;
	opt->min_asm_ovlp = 33;
	opt->min_merge_len = 0;

	opt->mag_opt.flag = MAG_F_NO_SIMPL | MAG_F_POPOPEN;
}


void fml_opt_adjust(fml_opt_t *opt, int n_seqs, const bseq1_t *seqs)
{
	int i, log_len;
	uint64_t tot_len = 0;
	if (opt->n_threads < 1) opt->n_threads = 1;
	for (i = 0; i < n_seqs; ++i) tot_len += seqs[i].l_seq; // compute total length
	for (log_len = 10; log_len < 32; ++log_len) // compute ceil(log2(tot_len))
		if (1ULL << log_len > tot_len) break;
	if (opt->ec_k == 0) opt->ec_k = (log_len + 12) / 2;
	if (opt->ec_k % 2 == 0) ++opt->ec_k;
	opt->mag_opt.min_elen = (int)((double)tot_len / n_seqs * 2.5 + .499);
}
