
#include <stdint.h>
#include <malloc.h>
#include <stdio.h>
#ifndef _MSC_VER
#include <alloca.h>
#endif
#include "err.h"
#include "utils.h"
#include "kmer-short.h"



static unsigned char _baseToShortBaseTable[256] = {
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 0, 4, 1, 6, 8, 8, 2, 5, 8, 8, 8, 8, 8, 7, 8,
	8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 0, 4, 1, 6, 8, 8, 2, 5, 8, 8, 8, 8, 8, 7, 8,
	8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};

static char _shortBaseToBaseTable[8] = {
	'A', 'C', 'G', 'T', 'B', 'H', 'D', 'N'
};


#ifdef KMER_SHORT_CHECKS
static void _kmer_short_test(const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	const uint64_t *x = KMer->Bases;
	uint64_t mask = (0ULL - 1ULL) - ((1ULL << KMerSize) - 1);

	assert((x[0] & mask) == 0);
	assert((x[1] & mask) == 0);
	assert((x[2] & mask) == 0);

	return;
}
#endif

static void _kmer_short_to_buffer(const uint32_t KMerSize, const KMER_SHORT *KMer, char *Buffer)
{
#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
#endif
	for (uint32_t i = 0; i < KMerSize; ++i)
		Buffer[i] = kmer_short_get_base(KMerSize, KMer, i);

	Buffer[KMerSize] = '\0';

	return;
}



void kmer_short_seq_init_by_sequence(PKMER_SHORT KMer, const uint32_t KMerSize, const char *Sequence)
{
	memset(KMer->Bases, 0, sizeof(KMer->Bases));
	for (uint32_t i = 0; i < KMerSize; ++i)
		kmer_short_set_base(KMerSize, KMer, i, Sequence[i]);

	return;
}


void kmer_short_seq_init_raw(PKMER_SHORT KMer, const uint32_t KMerSize, const uint64_t *Data)
{
	memcpy(&KMer->Bases, Data, sizeof(KMer->Bases));
#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
#endif

	return;
}


void kmer_short_advance(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	unsigned int c = 0;
	uint64_t *x = KMer->Bases;
	uint64_t mask = (1ULL << KMerSize) - 1;

#ifdef KMER_SHORT_CHECKS
	char buf1[KMER_SHORT_MAXIMUM_SIZE + 1];
	char buf2[KMER_SHORT_MAXIMUM_SIZE + 1];
	_kmer_short_to_buffer(KMerSize, KMer, buf1);
	_kmer_short_test(KMerSize, KMer);
#endif
	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (x[0] << 1 | (c & 1))  & mask;
	x[1] = (x[1] << 1 | ((c >> 1) & 1)) & mask;
	x[2] = (x[2] << 1 | (c >> 2)) & mask;

#ifdef KMER_SHORT_CHECKS
	_kmer_short_to_buffer(KMerSize, KMer, buf2);
	if (buf2[KMerSize - 1] != Base)
		__debugbreak();

	for (uint32_t i = 0; i < KMerSize - 2; ++i) {
		if (buf1[i + 1] != buf2[i])
			__debugbreak();
	}
#endif

	return;
}


void kmer_short_back(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	int c = 0;
	uint64_t *x = KMer->Bases;
	const uint32_t shift = KMerSize - 1;

#ifdef KMER_SHORT_CHECKS
	char buf1[KMER_SHORT_MAXIMUM_SIZE + 1];
	char buf2[KMER_SHORT_MAXIMUM_SIZE + 1];
	_kmer_short_to_buffer(KMerSize, KMer, buf1);
	_kmer_short_test(KMerSize, KMer);
#endif
	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = ((x[0] >> 1) | ((uint64_t)(c & 1) << shift));
	x[1] = ((x[1] >> 1) | ((uint64_t)((c >> 1) & 1) << shift));
	x[2] = ((x[2] >> 1) | ((uint64_t)(c >> 2) << shift));

#ifdef KMER_SHORT_CHECKS
	_kmer_short_to_buffer(KMerSize, KMer, buf2);
	if (buf2[0] != Base)
		__debugbreak();

	for (uint32_t i = 1; i < KMerSize - 1; ++i) {
		if (buf2[i + 1] != buf1[i])
			__debugbreak();
	}
#endif

	return;
}


char kmer_short_get_base(const uint32_t KMerSize, const KMER_SHORT *KMer, const uint32_t Pos)
{
	int c = 0;
	const uint64_t *data = KMer->Bases;
	const uint32_t d = KMerSize - Pos - 1;
	const uint64_t bit = (1ULL << d);

#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
#endif
	c = ((data[0] & bit) >> d) |
		(((data[1] & bit) >> d) << 1) |
		(((data[2] & bit) >> d) << 2);

	assert(c < 8);

	return _shortBaseToBaseTable[c];
}


char kmer_short_get_last_base(const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	int c = 0;
	const uint64_t *data = &KMer->Bases;

#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
#endif
	c = (data[0] & 1) |
		((data[1] & 1) << 1) |
		((data[2] & 1) << 2);

	assert(c < 8);
#ifdef KMER_SHORT_CHECKS
	char buf1[KMER_SHORT_MAXIMUM_SIZE + 1];
	_kmer_short_to_buffer(KMerSize, KMer, buf1);
	assert(_shortBaseToBaseTable[c] == buf1[KMerSize - 1]);
#endif

	return _shortBaseToBaseTable[c];
}


void kmer_short_set_base(const uint32_t KMerSize, PKMER_SHORT KMer, const uint32_t Pos, const char Base)
// d-bp from the 3'-end of k-mer; 0<=d<k
{ // IMPORTANT: 0 <= c < 4
	int c = 0;
	const uint32_t d = KMerSize - Pos - 1;
	const uint64_t t = ~(1ULL << d);
	uint64_t *x = KMer->Bases;

#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
	
	char buf1[KMER_SHORT_MAXIMUM_SIZE + 1];
	char buf2[KMER_SHORT_MAXIMUM_SIZE + 1];
	_kmer_short_to_buffer(KMerSize, KMer, buf1);
#endif

	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = ((uint64_t)(c & 1) << d) | (x[0] & t);
	x[1] = ((uint64_t)(((c >> 1) & 1)) << d) | (x[1] & t);
	x[2] = ((uint64_t)(c >> 2) << d) | (x[2] & t);

#ifdef KMER_SHORT_CHECKS
	_kmer_short_to_buffer(KMerSize, KMer, buf2);
	assert(buf2[Pos] == Base);
	for (uint32_t i = 0; i < KMerSize; ++i) {
		if (i == Pos)
			continue;

		if (buf1[i] != buf2[i])
			__debugbreak();
	}
#endif

	return;
}


size_t kmer_short_hash(const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	size_t ret = 0;
	const uint64_t *h = KMer->Bases;

#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, KMer);
#endif
	ret = (size_t)(
		(h[0] + h[1] + h[2]) |
		((h[0] << KMerSize) ^ (h[1] << KMerSize) ^ (h[2] << KMerSize)));

	return ret;
}


void kmer_short_print(FILE *Stream, const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	char *buf = alloca((KMerSize + 1) * sizeof(char));

	_kmer_short_to_buffer(KMerSize, KMer, buf);
	fputs(buf, Stream);

	return;
}


boolean kmer_short_seq_equal(const uint32_t KMerSize, const KMER_SHORT *K1, const KMER_SHORT *K2)
{
	const uint64_t *B1 = K1->Bases;
	const uint64_t *B2 = K2->Bases;

#ifdef KMER_SHORT_CHECKS
	_kmer_short_test(KMerSize, K1);
	_kmer_short_test(KMerSize, K2);
#endif
	return (
		B1[0] == B2[0] &&
		B1[1] == B2[1] &&
		B1[2] == B2[2]);
}
