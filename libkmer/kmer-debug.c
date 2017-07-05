
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer-debug.h"




/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


void kmer_init_by_sequence(PKMER KMer, const uint32_t KMerSize, const char *Sequence)
{
	if (Sequence != NULL)
		memcpy(KMer->Bases, Sequence, KMerSize*sizeof(char));

	return;
}


void kmer_advance(const uint32_t KMerSize, PKMER KMer, const char Base)
{
	assert(KMer->Size == KMerSize);

	memmove(KMer->Bases, KMer->Bases + 1, (KMerSize - 1)*sizeof(char));
	KMer->Bases[KMer->Size - 1] = Base;

	return;
}

void kmer_back(const uint32_t KMerSize, PKMER KMer, const char Base)
{
	assert(KMerSize == KMer->Size);

	memmove(KMer->Bases + 1, KMer->Bases, (KMerSize - 1)*sizeof(char));
	KMer->Bases[0] = Base;

	return;
}

boolean kmer_seq_equal(const KMER *K1, const KMER *K2)
{
	assert(K1->Size == K2->Size);

	return (memcmp(K1->Bases, K2->Bases, K1->Size*sizeof(char)) == 0);
}


void kmer_print(FILE *Stream, const uint32_t KMerSize, const KMER *KMer)
{
	for (size_t i = 0; i < KMerSize; ++i)
		fputc(kmer_get_base(KMer, i), Stream);

	fprintf(Stream, "_%u", KMer->Number);

	return;
}
