
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
	void *Data;
} KMER_TABLE_ENTRY, *PKMER_TABLE_ENTRY;

typedef struct _KMER_TABLE;

typedef void(KMER_TABLE_ON_INSERT_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order);
typedef void(KMER_TABLE_ON_DELETE_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData);
typedef ERR_VALUE (KMER_TABLE_ON_COPY_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, void **Copy);
typedef void(KMER_TABLE_ON_PRINT_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream);

typedef struct _KMER_TABLE_CALLBACKS {
	KMER_TABLE_ON_INSERT_CALLBACK *OnInsert;
	KMER_TABLE_ON_DELETE_CALLBACK *OnDelete;
	KMER_TABLE_ON_COPY_CALLBACK *OnCopy;
	KMER_TABLE_ON_PRINT_CALLBACK *OnPrint;
} KMER_TABLE_CALLBACKS, *PKMER_TABLE_CALLBACKS;

typedef struct _KMER_TABLE {
	size_t Size;
	size_t Inverse;
	size_t X;
	size_t PowX;
	size_t KMerSize;
	uint32_t LastOrder;
	KMER_TABLE_CALLBACKS Callbacks;
	PKMER_TABLE_ENTRY Entries;
} KMER_TABLE, *PKMER_TABLE;


ERR_VALUE kmer_table_create(const size_t KMerSize, const size_t X, const size_t Size, const KMER_TABLE_CALLBACKS *Callbacks, PKMER_TABLE *Table);
void kmer_table_destroy(PKMER_TABLE Table);
ERR_VALUE kmer_table_extend(PKMER_TABLE Table);
ERR_VALUE kmer_table_copy(const PKMER_TABLE Source, PKMER_TABLE * Copied);
void kmer_table_print(FILE *Stream, const PKMER_TABLE Table);

ERR_VALUE kmer_table_delete(PKMER_TABLE Table, const PKMER KMer);
ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const PKMER KMer, void *Data);
ERR_VALUE kmer_table_insert_hint(PKMER_TABLE Table, const PKMER KMer, const size_t Hash, void *Data);
void *kmer_table_get(const PKMER_TABLE Table, const PKMER KMer);

ERR_VALUE kmer_table_first(const PKMER_TABLE Table, PKMER_TABLE_ENTRY *Slot);
ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const PKMER_TABLE_ENTRY Current, PKMER_TABLE_ENTRY *Next);

size_t kmer_hash(const PKMER_TABLE Table, const PKMER KMer);
size_t kmer_hash_advance(const PKMER_TABLE Table, const PKMER KMer, const size_t Hash, const char NewBase);



#endif 
