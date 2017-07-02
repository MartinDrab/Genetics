
#ifndef FML_INTERNAL_H
#define FML_INTERNAL_H

#include "fml.h"
#include "kthread.h"

extern unsigned char seq_nt6_table[256];

#ifdef __cplusplus
extern "C" {
#endif

void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);
struct bfc_ch_s *fml_count(int n, const bseq1_t *seq, int k, int q, int l_pre, int n_threads);

#ifdef __cplusplus
}
#endif

#endif
