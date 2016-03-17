
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


#define kmer_get_size(aKMer)					((aKMer)->Size)
#define kmer_get_base(aKMer, aIndex)			((aKMer)->Bases[(aIndex)])
#define kmer_set_base(aKMer, aIndex, aBase)		((aKmer)->Bases[(aIndex)] = (aBase))
#define kmer_get_number(aKMer)					((aKMer)->Number)
#define kmer_set_number(aKMer, aNumber)			((aKMer)->Number = (aNumber))

PKMER kmer_alloc(const uint32_t Number, const uint32_t Size, const char *Sequence);
#define KMER_STACK_ALLOC(aVariable, aNumber, aSize, aSequence)				\
	{																		\
		aVariable = (PKMER)alloca(sizeof(KMER) + aSize*sizeof(char));		\
		aVariable->Size = aSize;											\
		aVariable->Number = aNumber;										\
		memcpy(aVariable->Bases, aSequence, aSize*sizeof(char));			\
	}														\

void kmer_init(PKMER KMer, const char *Sequence);
void kmer_init_from_kmer(PKMER Dest, const KMER *Source);
void kmer_free(PKMER KMer);
PKMER kmer_copy(const KMER *KMer);
void kmer_advance(PKMER KMer, const char Base);
void kmer_back(PKMER KMer, const char Base);
boolean kmer_equal(const KMER *K1, const KMER *K2);
boolean kmer_seq_equal(const KMER *K1, const KMER *K2);
void kmer_print(FILE *Stream, const KMER *KMer);


#endif 
