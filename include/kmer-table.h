
#ifndef __GASSM_KMER_TABLE_H__
#define __GASSM_KMER_TABLE_H__


#include "err.h"
#include "kmer.h"


typedef struct _KMER_TABLE_ENTRY_ADVANCED_INFO {
	union {
		struct {
			PKMER Input;
			PKMER Output;
		} PassThroughVertex;
	};
} KMER_TABLE_ENTRY_ADVANCED_INFO, *PKMER_TABLE_ENTRY_ADVANCED_INFO;

typedef struct _KMER_TABLE_ENTRY {
	unsigned char Deleted;
	PKMER KMer;
	uint32_t DegreeIn;
	uint32_t degreeOut;
	uint32_t Order;
	KMER_TABLE_ENTRY_ADVANCED_INFO AdvancedInfo;
} KMER_TABLE_ENTRY, *PKMER_TABLE_ENTRY;

typedef struct _KMER_TABLE {
	size_t Size;
	size_t Inverse;
	size_t X;
	size_t PowX;
	size_t KMerSize;
	uint32_t LastOrder;
	PKMER_TABLE_ENTRY Entries;
} KMER_TABLE, *PKMER_TABLE;


ERR_VALUE kmer_table_create(const size_t KMerSize, const size_t X, const size_t Size, PKMER_TABLE *Table);
void kmer_table_destroy(PKMER_TABLE Table);
ERR_VALUE kmer_table_extend(PKMER_TABLE Table);
ERR_VALUE kmer_table_copy(const PKMER_TABLE Source, PKMER_TABLE * Copied);
void kmer_table_print(FILE *Stream, const PKMER_TABLE Table);

ERR_VALUE kmer_table_delete(PKMER_TABLE Table, const PKMER KMer);
ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const PKMER KMer);
ERR_VALUE kmer_table_insert_hint(PKMER_TABLE Table, const PKMER KMer, const size_t Hash);
PKMER_TABLE_ENTRY kmer_table_get(const PKMER_TABLE Table, const PKMER KMer);

ERR_VALUE kmer_table_first(const PKMER_TABLE Table, PKMER_TABLE_ENTRY *Slot);
ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const PKMER_TABLE_ENTRY Current, PKMER_TABLE_ENTRY *Next);

size_t kmer_hash(const PKMER_TABLE Table, const PKMER KMer);
size_t kmer_hash_advance(const PKMER_TABLE Table, const PKMER KMer, const size_t Hash, const char NewBase);



#endif 
