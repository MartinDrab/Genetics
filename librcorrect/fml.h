
#ifndef __FML_H__
#define __FML_H__

#define FML_VERSION "r53"

#include <stdint.h>

typedef struct {
	int32_t l_seq;
	char *seq, *qual; // NULL-terminated strings; length expected to match $l_seq
} bseq1_t;

#define MAG_F_AGGRESSIVE 0x20 // pop variant bubbles (not default)
#define MAG_F_POPOPEN    0x40 // aggressive tip trimming (default)
#define MAG_F_NO_SIMPL   0x80 // skip bubble simplification (default)

typedef struct {
	int flag, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bdiff, max_bvtx, min_merge_len, trim_len, trim_depth;
	float min_dratio1, max_bcov, max_bfrac;
} magopt_t;

typedef struct {
	int n_threads;        // number of threads; don't use multi-threading for small data sets
	int ec_k;             // k-mer length for error correction; 0 for auto estimate
	int min_cnt, max_cnt; // both occ threshold in ec and tip threshold in cleaning lie in [min_cnt,max_cnt]
	int min_asm_ovlp;     // min overlap length during assembly
	int min_merge_len;    // during assembly, don't explicitly merge an overlap if shorter than this value
	magopt_t mag_opt;     // graph cleaning options
} fml_opt_t;

struct rld_t;
struct mag_t;

typedef struct {
	uint32_t len:31, from:1; // $from and $to: 0 meaning overlapping 5'-end; 1 overlapping 3'-end
	uint32_t id:31, to:1;    // $id: unitig number
} fml_ovlp_t;

typedef struct {
	int32_t len;      // length of sequence
	int32_t nsr;      // number of supporting reads
	char *seq;        // unitig sequence
	char *cov;        // cov[i]-33 gives per-base coverage at i
	int n_ovlp[2];    // number of 5'-end [0] and 3'-end [1] overlaps
	fml_ovlp_t *ovlp; // overlaps, of size n_ovlp[0]+n_ovlp[1]
} fml_utg_t;

extern int fm_verbose;

#ifdef __cplusplus
extern "C" {
#endif 

void fml_opt_init(fml_opt_t *opt);
void fml_opt_adjust(fml_opt_t *opt, int n_seqs, const bseq1_t *seqs);

float fml_correct(const fml_opt_t *opt, int n, bseq1_t *seq);
float fml_fltuniq(const fml_opt_t *opt, int n, bseq1_t *seq);
bseq1_t *bseq_read(const char *fn, int *n_);



#endif
