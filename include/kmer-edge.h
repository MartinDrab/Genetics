
#ifndef __KMER_EDGE_H__
#define __KMER_EDGE_H__


#include <stdint.h>
#include "err.h"
#include "khash.h"
#include "utils.h"
#include "kmer.h"


typedef struct _KMER_EDGE_TABLE_KEY {
	const KMER *Source;
	const KMER *Dest;
} KMER_EDGE_TABLE_KEY, *PKMER_EDGE_TABLE_KEY;

struct _KMER_EDGE_TABLE;

typedef void(KMER_EDGE_TABLE_ON_INSERT_CALLBACK)(struct _KMER_EDGE_TABLE *Table, void *ItemData, const uint32_t Order);
typedef void(KMER_EDGE_TABLE_ON_DELETE_CALLBACK)(struct _KMER_EDGE_TABLE *Table, void *ItemData, void *Context);
typedef ERR_VALUE(KMER_EDGE_TABLE_ON_COPY_CALLBACK)(struct _KMER_EDGE_TABLE *Table, void *ItemData, void **Copy, void *Context);
typedef void(KMER_EDGE_TABLE_ON_PRINT_CALLBACK)(struct _KMER_EDGE_TABLE *Table, void *ItemData, void *Context, FILE *Stream);

typedef struct _KMER_EDGE_TABLE_CALLBACKS {
	void *Context;
	KMER_EDGE_TABLE_ON_INSERT_CALLBACK *OnInsert;
	KMER_EDGE_TABLE_ON_DELETE_CALLBACK *OnDelete;
	KMER_EDGE_TABLE_ON_COPY_CALLBACK *OnCopy;
	KMER_EDGE_TABLE_ON_PRINT_CALLBACK *OnPrint;
} KMER_EDGE_TABLE_CALLBACKS, *PKMER_EDGE_TABLE_CALLBACKS;


__KHASH_TYPE(edgeTable, KMER_EDGE_TABLE_KEY, void *)


typedef struct _KMER_EDGE_TABLE {
	uint32_t KMerSize;
	unsigned int LastOrder;
	khash_t(edgeTable) *KHashTable;
	KMER_EDGE_TABLE_CALLBACKS Callbacks;
} KMER_EDGE_TABLE, *PKMER_EDGE_TABLE;




ERR_VALUE kmer_edge_table_create(const uint32_t KMerSize, const size_t Size, const PKMER_EDGE_TABLE_CALLBACKS Callbacks, PKMER_EDGE_TABLE *Table);
void kmer_edge_table_destroy(PKMER_EDGE_TABLE Table);

#define kmer_edge_table_size(aTable)	\
	(kh_size((aTable)->KHashTable))

ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const KMER *Source, const KMER *Dest, void *Data);
ERR_VALUE kmer_edge_table_delete(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
void *kmer_edge_table_get(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest);
ERR_VALUE kmer_edge_table_first(const PKMER_EDGE_TABLE Table, void **Slot, void **Data);
ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const void *Current, void **Next, void **Data);

void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table, void *Context);




#endif 
