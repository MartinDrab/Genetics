
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"




/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/

PKMER kmer_alloc(const uint32_t Size, const char *Sequence)
{
	PKMER ret = NULL;

	ret = (PKMER)malloc(sizeof(KMER) + Size*sizeof(char));
	if (ret != NULL) {
		ret->Size = Size;
		kmer_init(ret, Sequence);
	}

	return ret;
}

void kmer_init(PKMER KMer, const char *Sequence)
{
	if (Sequence != NULL)
		memcpy(KMer->Bases, Sequence, KMer->Size*sizeof(char));

	return;
}

void kmer_init_from_kmer(PKMER Dest, const PKMER Source)
{
	assert(Dest->Size == Source->Size);
	memcpy(Dest->Bases, Source->Bases, Dest->Size * sizeof(char));

	return;
}


void kmer_free(PKMER KMer)
{
	free(KMer);

	return;
}

PKMER kmer_copy(const PKMER KMer)
{
	return kmer_alloc(KMer->Size, KMer->Bases);
}


void kmer_advance(PKMER KMer, const char Base)
{
	memmove(KMer->Bases, KMer->Bases + 1, (KMer->Size - 1)*sizeof(char));
	KMer->Bases[KMer->Size - 1] = Base;

	return;
}

boolean kmer_equal(const PKMER K1, const PKMER K2)
{
	assert(K1->Size == K2->Size);

	return (memcmp(K1->Bases, K2->Bases, K1->Size*sizeof(char)) == 0);
}

void kmer_print(const PKMER KMer)
{
	for (size_t i = 0; i < KMer->Size; ++i)
		putchar(kmer_get_base(KMer, i));

	return;
}
