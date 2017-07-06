
#include "err.h"
#include "utils.h"
#include "kmer-short.h"



static unsigned char _baseToShortBaseTable[256] = {
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 0, 4, 1, 6, 5, 8, 2, 8, 8, 8, 8, 8, 8, 7, 8,
	8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 0, 4, 1, 6, 5, 8, 2, 8, 8, 8, 8, 8, 8, 7, 8,
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
	'A', 'C', 'G', 'T', 'B', 'E', 'D', 'N'
};



void kmer_short_seq_init_by_sequence(PKMER_SHORT KMer, const uint32_t KMerSize, const char *Sequence)
{
	for (uint32_t i = 0; i < KMerSize; ++i)
		kmer_short_set_base(KMerSize, KMer, i, Sequence[i]);

	return;
}


void kmer_short_seq_init_raw(PKMER_SHORT KMer, const uint32_t KMerSize, const uint64_t *Data)
{
	memcpy(&KMer->B, Data, sizeof(KMer->B));

	return;
}


void kmer_short_advance(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	unsigned int c = 0;
	uint64_t *x = KMer->B;
	uint64_t mask = (1ULL << KMerSize) - 1;

	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (x[0] << 1 | (c & 1))  & mask;
	x[1] = (x[1] << 1 | ((c >> 1) & 1)) & mask;
	x[2] = (x[2] << 1 | (c >> 2)) & mask;

	return;
}


void kmer_short_back(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	int c = 0;
	uint64_t *x = KMer->B;
	const uint64_t mask = (1ULL << KMerSize) - 1;

	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (x[0] >> 1 | (c & 1) << (KMerSize - 1))  & mask;
	x[1] = (x[1] >> 1 | ((c >> 1) & 1) << (KMerSize - 1)) & mask;
	x[2] = (x[2] >> 1 | (c >> 2) << (KMerSize - 1)) & mask;

	return;
}


char kmer_short_get_base(const uint32_t KMerSize, const KMER_SHORT *KMer, const uint32_t Pos)
{
	char b = 0;
	int c = 0;
	const uint64_t *data = &KMer->B;
	const uint64_t bit = (1 << Pos);

	c = ((data[0] & bit) >> Pos) |
		((data[0] & bit) >> (Pos - 1)) |
		((data[0] & bit) >> (Pos - 2) );

	assert(c < 8);

	return _shortBaseToBaseTable[c];
}

void kmer_short_set_base(const uint32_t KMerSize, PKMER_SHORT KMer, const uint32_t Pos, const char Base)
// d-bp from the 3'-end of k-mer; 0<=d<k
{ // IMPORTANT: 0 <= c < 4
	int c = 0;
	const uint64_t t = ~(1ULL << Pos);
	uint64_t *x = KMer->B;

	c = _baseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (uint64_t)(c & 1) << Pos | (x[0] & t);
	x[1] = (uint64_t)(((c >> 1) & 1)) << Pos | (x[1] & t);
	x[2] = (uint64_t)(c >> 2) << Pos | (x[2] & t);

	return;
}


// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
static uint64_t kmer_short_hash_64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;

	return key;
}


size_t kmer_short_hash(const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	int t = KMerSize >> 1, u = ((KMer->B[1] >> t & 1) > (KMer->B[0] >> t & 1)); // the middle base is always different
	uint64_t mask = (1ULL << KMerSize) - 1, ret;
	uint64_t h[3];

	h[0] = kmer_short_hash_64((KMer->B[u << 1 | 0] + KMer->B[u << 1 | 1]) & mask, mask);
	h[1] = kmer_short_hash_64(h[0] ^ KMer->B[u << 1 | 1], mask);
	h[2] = kmer_short_hash_64(h[1] ^ KMer->B[u << 2 | 1], mask);

	ret = (h[2] ^ h[0] ^ h[1]) << KMerSize | (h[2] + (h[0] + h[1]) & mask);

	return (size_t)ret;
}


void kmer_short_print(FILE *Stream, const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	char *buf = alloca((KMerSize + 1) * sizeof(char));

	for (size_t l = 0; l < KMerSize; ++l)
		buf[KMerSize - 1 - l] = "ACGT"[(((KMer->B[2] >> 2) & 1) << 2) | (KMer->B[1] >> l & 1) << 1 | (KMer->B[0] >> l & 1)];

	buf[KMerSize] = 0;
	fputs(buf, Stream);

	return;
}


boolean kme_short_seq_equal(const uint32_t KMerSize, const KMER_SHORT *K1, const KMER_SHORT *K2)
{
	const uint64_t *B1 = K1->B;
	const uint64_t *B2 = K2->B;

	return (
		B1[0] == B2[0] &&
		B1[1] == B2[1] &&
		B1[2] == B2[2]);
}
