
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer-debug.h"




/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


void kmer_debug_seq_init_by_sequence(PKMER_DEBUG KMer, const uint32_t KMerSize, const char *Sequence)
{
	if (Sequence != NULL)
		memcpy(KMer->Bases, Sequence, KMerSize*sizeof(char));

	return;
}


void kmer_debug_advance(const uint32_t KMerSize, PKMER_DEBUG KMer, const char Base)
{
	assert(KMer->Size == KMerSize);
	assert(Base != '\0');

	memmove(KMer->Bases, KMer->Bases + 1, (KMerSize - 1)*sizeof(char));
	KMer->Bases[KMer->Size - 1] = Base;

	return;
}

void kmer_debug_back(const uint32_t KMerSize, PKMER_DEBUG KMer, const char Base)
{
	assert(KMerSize == KMer->Size);

	memmove(KMer->Bases + 1, KMer->Bases, (KMerSize - 1)*sizeof(char));
	KMer->Bases[0] = Base;

	return;
}

boolean kmer_debug_seq_equal(const uint32_t KMerSize, const KMER_DEBUG *K1, const KMER_DEBUG *K2)
{
	assert(K1->Size == K2->Size);
	assert(KMerSize == K1->Size);

	return (memcmp(K1->Bases, K2->Bases, KMerSize*sizeof(char)) == 0);
}


void kmer_debug_print(FILE *Stream, const uint32_t KMerSize, const KMER_DEBUG *KMer)
{
	for (size_t i = 0; i < KMerSize; ++i)
		fputc(kmer_debug_get_base(KMerSize, KMer, i), Stream);

	fprintf(Stream, "_%u", KMer->Number);

	return;
}


size_t kmer_debug_hash(const uint32_t Context, const KMER_DEBUG *KMer)
{
	size_t hash = 0;
	const uint32_t kmerSize = (uint32_t)Context;

	assert(kmerSize == kmer_debug_get_size(KMer));
	for (size_t i = 0 + 1; i < kmerSize; ++i)
		hash = (hash << 5) - hash + kmer_debug_get_base(Context, KMer, i);

	return hash;
}
