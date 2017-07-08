
#ifndef __KMER_BASE_H__
#define __KMER_BASE_H__




#include "kmer.h"



INLINE_FUNCTION ERR_VALUE kmer_alloc(const uint32_t Number, const uint32_t Size, const char *Sequence, PKMER *KMer)
{
	PKMER tmpKMer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc(KMER_BYTES(Size), (void **)&tmpKMer);
	if (ret == ERR_SUCCESS) {
		kmer_set_number(tmpKMer, Number);
		kmer_set_size(tmpKMer, Size);
		kmer_seq_init_by_sequence(tmpKMer, Size, Sequence);
		*KMer = tmpKMer;
	}

	return ret;
}


INLINE_FUNCTION void kmer_init_by_base(PKMER KMer, const uint32_t KMerSize, char Base)
{
	for (uint32_t i = 0; i < KMerSize; ++i)
		kmer_set_base(KMerSize, KMer, i, Base);

	return;
}


INLINE_FUNCTION void kmer_init_from_kmer(PKMER Dest, const uint32_t KMerSize,  const KMER *Source)
{
// 	assert(KMerSize == Source->Size);
	kmer_set_size(Dest, kmer_get_size(Source));
	kmer_set_number(Dest, kmer_get_number(Source));
	kmer_seq_init_raw(Dest, KMerSize, Source->Bases);

	return;
}


INLINE_FUNCTION void kmer_free(PKMER KMer)
{
	utils_free(KMer);

	return;
}


INLINE_FUNCTION ERR_VALUE kmer_copy(PKMER *Dest, const uint32_t KMerSize, const KMER *KMer)
{
	PKMER tmpKmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc(KMER_BYTES(KMerSize), (void **)&tmpKmer);
	if (ret == ERR_SUCCESS) {
		kmer_set_size(tmpKmer, kmer_get_size(KMer));
		kmer_set_number(tmpKmer, kmer_get_number(KMer));
		kmer_seq_init_raw(tmpKmer, KMerSize, KMer->Bases);
		*Dest = tmpKmer;
	}

	return ret;
}


INLINE_FUNCTION boolean kmer_equal(uint32_t KMerSize, const KMER *K1, const KMER *K2)
{
//	assert(K1->Size == K2->Size);
//	assert(KMerSize == K1->Size);
	return (kmer_get_number(K1) == kmer_get_number(K2) && kmer_seq_equal(KMerSize, K1, K2));
}


#define KMER_STACK_ALLOC(aVariable, aNumber, aSize, aSequence)				\
	{																		\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));		\
		memset(aVariable, 0, KMER_BYTES(aSize));	\
		kmer_set_size(aVariable, (aSize));										\
		kmer_set_number(aVariable, (aNumber));	\
		if ((aSequence)!= NULL)	\
			kmer_seq_init_by_sequence(aVariable, aSize, aSequence);		\
	}														\

#define KMER_STACK_ALLOC_FROM_KMER(aVariable, aSize, aKMer)				\
	{																	\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));					\
		memset(aVariable, 0, KMER_BYTES(aSize));	\
		kmer_set_size(aVariable, (aSize));										\
		kmer_init_from_kmer(aVariable, aSize, aKMer);							\
	}																	\


#endif
