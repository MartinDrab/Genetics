
#ifndef __GASSM_KMER_H__
#define __GASSM_KMER_H__


#include <stdint.h>
#include <malloc.h>
#include "err.h"
#include "utils.h"



typedef struct _KMER {
	uint32_t Size;
	char Bases[1];
} KMER, *PKMER;


#define kmer_get_size(aKMer)					((aKMer)->Size)
#define kmer_get_base(aKMer, aIndex)			((aKMer)->Bases[aIndex])
#define kmer_set_base(aKMer, aIndex, aBase)		((aKmer)->Bases[aIndex] = aBase)

PKMER kmer_alloc(const uint32_t Size, const char *Sequence);
#define KMER_STACK_ALLOC(aVariable, aSize, aSequence)			\
	{															\
		aVariable = (PKMER)alloca(sizeof(KMER) + aSize*sizeof(char));		\
		aVariable->Size = aSize;									\
		memcpy(aVariable->Bases, aSequence, aSize*sizeof(char));	\
	}														\

void kmer_init(PKMER KMer, const char *Sequence);
void kmer_init_from_kmer(PKMER Dest, const PKMER Source);
void kmer_free(PKMER KMer);
PKMER kmer_copy(const PKMER KMer);
void kmer_advance(PKMER KMer, const char Base);
void kmer_back(PKMER KMer, const char Base);
boolean kmer_equal(const PKMER K1, const PKMER K2);
void kmer_print(const PKMER KMer);


#endif 
