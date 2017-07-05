
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer-debug.h"




/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


void kmer_init_by_sequence(PKMER KMer, const char *Sequence)
{
	if (Sequence != NULL)
		memcpy(KMer->Bases, Sequence, KMer->Size*sizeof(char));

	return;
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


void kmer_print(FILE *Stream, const KMER *KMer)
{
	for (size_t i = 0; i < KMer->Size; ++i)
		fputc(kmer_get_base(KMer, i), Stream);

	fprintf(Stream, "_%u", KMer->Number);

	return;
}
