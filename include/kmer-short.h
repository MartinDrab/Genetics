
#ifndef __GASSM_KMER_SHORT_H__
#define __GASSM_KMER_SHORT_H__



#include <stdint.h>
#include "utils.h"


typedef struct _KMER_SHORT {
	uint32_t Number;
	uint64_t B[3];
} KMER_SHORT, *PKMER_SHORT;


#define KMER_SHORT_BYTES(aKMerSize)						(sizeof(KMER_SHORT))
#define KMER_SHORT_BYTES_EXTRA(aKMerSize, aExtra)		(KMER_BYTES(aKMerSize) + aExtra)


#define kmer_short_get_size(aKMer)						(0)
#define kmer_sort_set_size(aKMer, aSize)				(0)
#define kmer_short_get_number(aKMer)					((aKMer)->Number)
#define kmer_short_set_number(aKMer, aNumber)			((aKMer)->Number = (aNumber))


void kmer_short_init_by_sequence(PKMER_SHORT KMer, const uint32_t KMerSize, const char *Sequence);;
void kmer_short_advance(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base);
void kmer_short_back(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base);
void kmer_short_set_base(const uint32_t KMerSize, PKMER_SHORT KMer, const uint32_t Pos, const char Base);
char kmer_short_get_base(const uint32_t KMerSize, const KMER_SHORT *KMer, const uint32_t Pos);
uint64_t kmer_short_hash(const uint32_t KMerSize, const KMER_SHORT *KMer);
void kmer_short_print(FILE *Stream, const uint32_t KMerSize, const KMER_SHORT *KMer);
boolean kme_short_seq_equal(const uint32_t KMerSize, const KMER_SHORT *K1, const KMER_SHORT *K2);




#endif
