
#ifndef __GASSM_KMER_TABLE_H__
#define __GASSM_KMER_TABLE_H__


#include "khash.h"
#include "err.h"
#include "kmer.h"


__KHASH_TYPE(vertexTable, const KMER *, void *)

struct _KMER_TABLE;

typedef void(KMER_TABLE_ON_INSERT_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order);
typedef void(KMER_TABLE_ON_DELETE_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, void *Context);
typedef ERR_VALUE (KMER_TABLE_ON_COPY_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, void **Copy, void *Context);
typedef void(KMER_TABLE_ON_PRINT_CALLBACK)(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream);

typedef struct _KMER_TABLE_CALLBACKS {
	void *Context;
	KMER_TABLE_ON_INSERT_CALLBACK *OnInsert;
	KMER_TABLE_ON_DELETE_CALLBACK *OnDelete;
	KMER_TABLE_ON_COPY_CALLBACK *OnCopy;
	KMER_TABLE_ON_PRINT_CALLBACK *OnPrint;
} KMER_TABLE_CALLBACKS, *PKMER_TABLE_CALLBACKS;


typedef struct _KMER_TABLE {
	uint32_t KMerSize;
	uint32_t LastOrder;
	KMER_TABLE_CALLBACKS Callbacks;
	khash_t(vertexTable) *KHashTable;
} KMER_TABLE, *PKMER_TABLE;


ERR_VALUE kmer_table_create(const uint32_t KMerSize, const size_t Size, const KMER_TABLE_CALLBACKS *Callbacks, PKMER_TABLE *Table);
void kmer_table_destroy(PKMER_TABLE Table);
void kmer_table_print(FILE *Stream, const PKMER_TABLE Table);

#define kmer_table_size(aTable)	\
	(kh_size((aTable)->KHashTable))

ERR_VALUE kmer_table_delete(PKMER_TABLE Table, const KMER *KMer);
ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const KMER *KMer, void *Data);
void *kmer_table_get(const struct _KMER_TABLE *Table, const struct _KMER *KMer);

ERR_VALUE kmer_table_first(const PKMER_TABLE Table, void **Slot, void **Data);
ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const void *Current, void **Next, void **Data);



#endif 
