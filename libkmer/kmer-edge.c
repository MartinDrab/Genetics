
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "khash.h"
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"




/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE_TABLE)


static INLINE_FUNCTION boolean _kmer_edge_key_equal(const KMER_EDGE_TABLE_KEY Key1, const KMER_EDGE_TABLE_KEY Key2)
{
	return (
		kmer_get_number(Key1.Source) == kmer_get_number(Key2.Source) &&
		kmer_get_number(Key1.Dest) == kmer_get_number(Key2.Dest) &&
		kmer_seq_equal(Key1.Source, Key2.Source) &&
		kmer_seq_equal(Key1.Dest, Key2.Dest)
		);
}

static INLINE_FUNCTION size_t _kmer_edge_key_hash(const KMER_EDGE_TABLE_KEY Key)
{
	size_t hash1 = Key.Source->Number;
	size_t hash2 = Key.Dest->Number;
	const size_t kmerSize = Key.Source->Size;

	for (size_t i = 0 + 1; i < kmerSize; ++i) {
		hash1 = (hash1 << 5) - hash1 + kmer_get_base(Key.Source, i);
		hash2 = (hash2 << 5) - hash2 + kmer_get_base(Key.Dest, i);
	}

	return (((Key.Source->Number + 1)*hash1 << 1) + (Key.Dest->Number + 1)*hash2);
}


__KHASH_IMPL(edgeTable, INLINE_FUNCTION, KMER_EDGE_TABLE_KEY, void *, TRUE, _kmer_edge_key_hash, _kmer_edge_key_equal)


static void _on_insert_dummy_callback(struct _KMER_EDGE_TABLE *Table, void *ItemData, const uint32_t Order)
{
	return;
}


static void _on_delete_dummy_callback(struct _KMER_EDGE_TABLE *Table, void *ItemData)
{
	return;
}


static ERR_VALUE _on_copy_dummy_callback(struct _KMER_EDGE_TABLE *Table, void *ItemData, void **Copy)
{
	*Copy = ItemData;

	return ERR_SUCCESS;
}


static void _on_print_dummy_callback(struct _KMER_EDGE_TABLE *Table, void *ItemData, FILE *Stream)
{
	return;
}



/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t Size, const PKMER_EDGE_TABLE_CALLBACKS Callbacks, PKMER_EDGE_TABLE *Table)
{
	PKMER_EDGE_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE_TABLE(&tmpTable);
	if (ret == ERR_SUCCESS) {
		tmpTable->KHashTable = kh_init(edgeTable);
		tmpTable->KMerSize = KMerSize;
		if (Callbacks == NULL) {
			tmpTable->Callbacks.OnInsert = _on_insert_dummy_callback;
			tmpTable->Callbacks.OnDelete = _on_delete_dummy_callback;
			tmpTable->Callbacks.OnCopy = _on_copy_dummy_callback;
			tmpTable->Callbacks.OnPrint = _on_print_dummy_callback;
		} else tmpTable->Callbacks = *Callbacks;
			
		tmpTable->LastOrder = 0;
		*Table = tmpTable;
		ret = ERR_SUCCESS;
	}

	return ret;
}


void kmer_edge_table_destroy(PKMER_EDGE_TABLE Table)
{	
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it))
			Table->Callbacks.OnDelete(Table, kh_val(Table->KHashTable, it));
	}

	kh_destroy(edgeTable, Table->KHashTable);
	utils_free(Table);

	return;
}


void *kmer_edge_table_get(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest)
{
	khiter_t it;
	void *ret = NULL;
	KMER_EDGE_TABLE_KEY key;

	key.Source = Source;
	key.Dest = Dest;
	it = kh_get(edgeTable, Table->KHashTable, key);
	if (it != kh_end(Table->KHashTable))
		ret = kh_val(Table->KHashTable, it);

	return ret;
}


ERR_VALUE kmer_edge_table_delete(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	khiter_t it;
	KMER_EDGE_TABLE_KEY key;

	key.Source = Source;
	key.Dest = Dest;
	it = kh_get(edgeTable, Table->KHashTable, key);
	if (it != kh_end(Table->KHashTable)) {
		Table->Callbacks.OnDelete(Table, kh_val(Table->KHashTable, it));
		kh_del(edgeTable, Table->KHashTable, it);
		ret = ERR_SUCCESS;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const KMER *Source, const KMER *Dest, void *Data)
{
	int r;
	khiter_t it;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	KMER_EDGE_TABLE_KEY key;

	key.Source = Source;
	key.Dest = Dest;
	it = kh_put(edgeTable, Table->KHashTable, key, &r);
	switch (r) {
		case 0:
			ret = ERR_ALREADY_EXISTS;
			break;
		case 1:
		case 2:
			kh_val(Table->KHashTable, it) = Data;
			Table->Callbacks.OnInsert(Table, Data, Table->LastOrder);
			++Table->LastOrder;
			ret = ERR_SUCCESS;
			break;
		case -1:
			ret = ERR_OUT_OF_MEMORY;
			break;
		default:
			__debugbreak();
			ret = ERR_INTERNAL_ERROR;
			break;
	}

	return ret;
}


void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table)
{
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it))
			Table->Callbacks.OnPrint(Table, kh_val(Table->KHashTable, it), Stream);
	}

	return;
}


ERR_VALUE kmer_edge_table_first(const PKMER_EDGE_TABLE Table, void **Slot, void **Data)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it)) {
			*Slot = (void *)it;
			*Data = kh_val(Table->KHashTable, it);
			ret = ERR_SUCCESS;
		}
	}

	return ret;
}


ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const void *Current, void **Next, void **Data)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = ((khiter_t)Current) + 1; it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it)) {
			*Next = (void *)it;
			*Data = kh_val(Table->KHashTable, it);
			ret = ERR_SUCCESS;
		}
	}

	return ret;
}
