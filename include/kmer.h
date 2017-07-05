
#ifndef __GASSM_KMER_H__
#define __GASSM_KMER_H__


#include <stdint.h>
#include <malloc.h>
#include <stdio.h>
#ifndef _MSC_VER
#include <alloca.h>
#endif
#include "err.h"
#include "utils.h"



typedef struct _KMER {
	uint32_t Number;
	uint32_t Size;
	char Bases[1];
} KMER, *PKMER;

#define KMER_BYTES(aKMerSize)					(sizeof(KMER)+aKMerSize*sizeof(char))
#define KMER_BYTES_EXTRA(aKMerSize, aExtra)		(KMER_BYTES(aKMerSize) + aExtra)


#define kmer_get_size(aKMer)					((aKMer)->Size)
#define kmer_get_base(aKMer, aIndex)			((aKMer)->Bases[(aIndex)])
#define kmer_get_last_base(aKMer)				kmer_get_base(aKMer, kmer_get_size(aKMer) - 1)
#define kmer_set_base(aKMer, aIndex, aBase)		((aKMer)->Bases[(aIndex)] = (aBase))
#define kmer_get_number(aKMer)					((aKMer)->Number)
#define kmer_set_number(aKMer, aNumber)			((aKMer)->Number = (aNumber))

ERR_VALUE kmer_alloc(const uint32_t Number, const uint32_t Size, const char *Sequence, PKMER *KMer);
#define KMER_STACK_ALLOC(aVariable, aNumber, aSize, aSequence)				\
	{																		\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));		\
		aVariable->Size = (aSize);											\
		aVariable->Number = (aNumber);										\
		kmer_init_by_sequence(aVariable, aSequence);		\
	}														\

#define KMER_STACK_ALLOC_FROM_KMER(aVariable, aSize, aKMer)				\
	{																	\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));					\
		aVariable->Size = (aSize);										\
		kmer_init_from_kmer(aVariable, aKMer);							\
	}																	\

void kmer_init_by_base(PKMER KMer, char Base);
void kmer_init_by_sequence(PKMER KMer, const char *Sequence);
void kmer_init_from_kmer(PKMER Dest, const KMER *Source);
void kmer_free(PKMER KMer);
ERR_VALUE kmer_copy(PKMER *Dest, const KMER *KMer);
void kmer_advance(PKMER KMer, const char Base);
void kmer_back(PKMER KMer, const char Base);
boolean kmer_equal(const KMER *K1, const KMER *K2);
boolean kmer_seq_equal(const KMER *K1, const KMER *K2);
void kmer_print(FILE *Stream, const KMER *KMer);


INLINE_FUNCTION size_t kmer_hash(const struct _KMER *KMer)
{
	size_t hash = 0;
	const size_t kmerSize = kmer_get_size(KMer);

	for (size_t i = 0 + 1; i < kmerSize; ++i)
		hash = (hash << 5) - hash + kmer_get_base(KMer, i);

	return hash;
}




#endif 
