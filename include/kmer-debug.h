
#ifndef __GASSM_KMER_DEBUG_H__
#define __GASSM_KMER_DEBUG_H__



typedef struct _KMER {
	uint32_t Number;
	uint32_t Size;
	char Bases[1];
} KMER, *PKMER;

#define KMER_BYTES(aKMerSize)					(sizeof(KMER)+aKMerSize*sizeof(char))
#define KMER_BYTES_EXTRA(aKMerSize, aExtra)		(KMER_BYTES(aKMerSize) + aExtra)


#define kmer_get_size(aKMer)					((aKMer)->Size)
#define kmer_set_size(aKMer, aSize)				((aKMer)->Size = (aSize))
#define kmer_get_base(aKMer, aIndex)			((aKMer)->Bases[(aIndex)])
#define kmer_get_last_base(aKMer)				kmer_get_base(aKMer, kmer_get_size(aKMer) - 1)
#define kmer_set_base(aKMer, aIndex, aBase)		((aKMer)->Bases[(aIndex)] = (aBase))
#define kmer_get_number(aKMer)					((aKMer)->Number)
#define kmer_set_number(aKMer, aNumber)			((aKMer)->Number = (aNumber))

void kmer_init_by_sequence(PKMER KMer, const uint32_t KMerSize, const char *Sequence);
void kmer_advance(const uint32_t KMerSize, PKMER KMer, const char Base);
void kmer_back(const uint32_t KMerSize, PKMER KMer, const char Base);
boolean kmer_seq_equal(const uint32_t KMerSize, const KMER *K1, const KMER *K2);
void kmer_print(FILE *Stream, const uint32_t KMerSize, const KMER *KMer);
size_t kmer_hash(const uint32_t Context, const KMER *KMer);






#endif 
