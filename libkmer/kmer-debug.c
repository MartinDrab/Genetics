
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer-debug.h"




/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/

ERR_VALUE kmer_alloc(const uint32_t Number, const uint32_t Size, const char *Sequence, PKMER *KMer)
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


void kmer_init_by_sequence(PKMER KMer, const char *Sequence)
{
	if (Sequence != NULL)
		memcpy(KMer->Bases, Sequence, KMer->Size*sizeof(char));

	return;
}


void kmer_init_by_base(PKMER KMer, char Base)
{
	for (size_t i = 0; i < KMer->Size; ++i)
		kmer_set_base(KMer, i, Base);
}


void kmer_init_from_kmer(PKMER Dest, const KMER *Source)
{
	Dest->Size = Source->Size;
	Dest->Number = Source->Number;
	kmer_init_by_sequence(Dest, Source->Bases);

	return;
}


void kmer_free(PKMER KMer)
{
	utils_free(KMer);

	return;
}

ERR_VALUE kmer_copy(PKMER *Dest, const KMER *KMer)
{
	return kmer_alloc(KMer->Number, KMer->Size, KMer->Bases, Dest);
}


void kmer_advance(PKMER KMer, const char Base)
{
	memmove(KMer->Bases, KMer->Bases + 1, (KMer->Size - 1)*sizeof(char));
	KMer->Bases[KMer->Size - 1] = Base;

	return;
}

void kmer_back(PKMER KMer, const char Base)
{
	memmove(KMer->Bases + 1, KMer->Bases, (KMer->Size - 1)*sizeof(char));
	KMer->Bases[0] = Base;

	return;
}

boolean kmer_seq_equal(const KMER *K1, const KMER *K2)
{
	assert(K1->Size == K2->Size);

	return (memcmp(K1->Bases, K2->Bases, K1->Size*sizeof(char)) == 0);
}


boolean kmer_equal(const KMER *K1, const KMER *K2)
{
	assert(K1->Size == K2->Size);

	return (K1->Number == K2->Number && kmer_seq_equal(K1, K2));
}

void kmer_print(FILE *Stream, const KMER *KMer)
{
	for (size_t i = 0; i < KMer->Size; ++i)
		fputc(kmer_get_base(KMer, i), Stream);

	fprintf(Stream, "_%u", KMer->Number);

	return;
}
