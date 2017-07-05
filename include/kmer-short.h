
#ifndef __GASSM_KMER_SHORT_H__
#define __GASSM_KMER_SHORT_H__



#include <stdint.h>
#include "utils.h"


typedef struct _KMER_SHORT {
	uint32_t Number;
	uint64_t B[3];
} KMER_SHORT, *PKMER_SHORT;


static unsigned char BaseToShortBaseTable[256] = {
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 
};


INLINE_FUNCTION void kmer_short_advance(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	unsigned int c = 0;
	uint64_t *x = KMer->B;
	uint64_t mask = (1ULL<<KMerSize) - 1;

	c = BaseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (x[0] << 1 | (c&1))  & mask;
	x[1] = (x[1] << 1 | ((c>>1)&1)) & mask;
	x[2] = (x[2] << 1 | (c >> 2)) & mask;

	return;
}


INLINE_FUNCTION void kmer_short_back(const uint32_t KMerSize, PKMER_SHORT KMer, const char Base)
{
	int c = 0;
	uint64_t *x = KMer->B;
	uint64_t mask = (1ULL << KMerSize) - 1;

	c = BaseToShortBaseTable[Base];
	assert(c < 8);
	x[0] = (x[0] >> 1 | (c & 1) << (KMerSize - 1))  & mask;
	x[1] = (x[1] >> 1 | ((c >> 1) & 1) << (KMerSize - 1)) & mask;
	x[2] = (x[2] >> 1 | (c >> 2) << (KMerSize - 1)  ) & mask;

	return;
}


INLINE_FUNCTION void kmer_short_set_base(const uint32_t KMerSize, PKMER_SHORT KMer, const uint32_t Pos, const char Base)
// d-bp from the 3'-end of k-mer; 0<=d<k
{ // IMPORTANT: 0 <= c < 4
	uint64_t t = ~(1ULL<<Pos);
	int c = 0;
	uint64_t *x = KMer->B;

	// TODO: Translate the Base into c
	x[0] = (uint64_t) (c&1)<<Pos | (x[0]&t);
	x[1] = (uint64_t)((c>>1))<<Pos | (x[1]&t);
	x[2] = (uint64_t)(c >> 2) << Pos | (x[2] & t);

	return;
}


// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
INLINE_FUNCTION uint64_t kmer_short_hash_64(uint64_t key, uint64_t mask)
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


INLINE_FUNCTION uint64_t kmer_short_hash(const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	int t = KMerSize>>1, u = ((KMer->B[1]>>t&1) > (KMer->B[0]>>t&1)); // the middle base is always different
	uint64_t mask = (1ULL<<KMerSize) - 1, ret;
	uint64_t h[3];

	h[0] = kmer_short_hash_64((KMer->B[u<<1|0] + KMer->B[u<<1|1]) & mask, mask);
	h[1] = kmer_short_hash_64(h[0] ^ KMer->B[u<<1|1], mask);
	h[2] = kmer_short_hash_64(h[1] ^ KMer->B[u << 2 | 1], mask);
	
	ret = (h[2] ^ h[0] ^ h[1]) << KMerSize | ((h[0] + h[1]) & mask);
	
	return ret;
}


INLINE_FUNCTION void kmer_short_print(FILE *Stream, const uint32_t KMerSize, const KMER_SHORT *KMer)
{
	char *buf = alloca((KMerSize + 1)*sizeof(char));

	for (size_t l = 0; l < KMerSize; ++l)
		buf[KMerSize - 1 - l] = "ACGT"[(((KMer->B[2]>>2)&1)<<2) | (KMer->B[1]>>l&1)<<1 | (KMer->B[0]>>l&1)];
	
	buf[KMerSize] = 0;
	
	return;
}


INLINE_FUNCTION boolean kme_short_seq_equal(const uint32_t KMerSize, const KMER_SHORT *K1, const KMER_SHORT *K2)
{
	return (memcmp(K1->B, K2->B, sizeof(K1->B)));
}



#endif
