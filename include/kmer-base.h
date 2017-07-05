
#ifndef __KMER_BASE_H__
#define __KMER_BASE_H__




#include "kmer.h"



INLINE_FUNCTION ERR_VALUE kmer_alloc(const uint32_t Number, const uint32_t Size, const char *Sequence, PKMER *KMer)
{
	PKMER tmpKMer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc(KMER_BYTES(Size), &tmpKMer);
	if (ret == ERR_SUCCESS) {
		kmer_set_number(tmpKMer, Number);
		tmpKMer->Size = Size;
		kmer_init_by_sequence(tmpKMer, Sequence);
		*KMer = tmpKMer;
	}

	return ret;
}


INLINE_FUNCTION void kmer_init_by_base(PKMER KMer, char Base)
{
	for (size_t i = 0; i < KMer->Size; ++i)
		kmer_set_base(KMer, i, Base);
}


INLINE_FUNCTION void kmer_init_from_kmer(PKMER Dest, const KMER *Source)
{
	Dest->Size = Source->Size;
	Dest->Number = Source->Number;
	kmer_init_by_sequence(Dest, Source->Bases);

	return;
}


INLINE_FUNCTION void kmer_free(PKMER KMer)
{
	utils_free(KMer);

	return;
}


INLINE_FUNCTION ERR_VALUE kmer_copy(PKMER *Dest, const KMER *KMer)
{
	return kmer_alloc(KMer->Number, KMer->Size, KMer->Bases, Dest);
}


INLINE_FUNCTION boolean kmer_equal(const KMER *K1, const KMER *K2)
{
	assert(K1->Size == K2->Size);

	return (kmer_get_number(K1) == kmer_get_number(K2) && kmer_seq_equal(K1, K2));
}


#define KMER_STACK_ALLOC(aVariable, aNumber, aSize, aSequence)				\
	{																		\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));		\
		kmer_set_size(aVariable, (aSize));										\
		kmer_set_number(aVariable, (aNumber));											\
		kmer_init_by_sequence(aVariable, aSequence);		\
	}														\

#define KMER_STACK_ALLOC_FROM_KMER(aVariable, aSize, aKMer)				\
	{																	\
		aVariable = (PKMER)alloca(KMER_BYTES(aSize));					\
		kmer_set_size(aVariable, (aSize));										\
		kmer_init_from_kmer(aVariable, aKMer);							\
	}																	\


#endif
