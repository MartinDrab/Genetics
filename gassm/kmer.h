
#ifndef __GASSM_KMER_H__
#define __GASSM_KMER_H__


#include <stdint.h>
#include "err.h"




typedef struct _KMER {
	uint32_t Size;
	char Bases[1];
} KMER, *PKMER;


#define kmer_get_size(aKMer)					((aKMer)->Size)
#define kmer_get_base(aKMer, aIndex)			((aKMer)->Bases[aIndex])

PKMER kmer_alloc(const uint32_t Size, const char *Sequence);
void kmer_init(PKMER KMer, const char *Sequence);
void kmer_free(PKMER KMer);
PKMER kmer_copy(const PKMER KMer);
void kmer_advance(PKMER KMer, const char Base);
boolean kmer_equal(const PKMER K1, const PKMER K2);



#endif 
