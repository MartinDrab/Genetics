
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "khash.h"
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"


/************************************************************************/
/*                      HELPER MACROS AND TYPES                         */
/************************************************************************/

/************************************************************************/
/*                   HELPER FUNCTIONS                                   */
/************************************************************************/


UTILS_TYPED_MALLOC_FUNCTION(KMER_TABLE)


__KHASH_IMPL(vertexTable, INLINE_FUNCTION, const KMER *, void *, TRUE, kmer_hash, kmer_equal)


static void _on_insert_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order)
{
	return;
}


static void _on_delete_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, void *Context)
{
	return;
}


static ERR_VALUE _on_copy_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, void **Copy, void *Context)
{
	*Copy = ItemData;

	return ERR_SUCCESS;
}


static void _on_print_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream)
{
	return;
}


/************************************************************************/
/*                  PUBLIC FUNCTIONS                                    */
/************************************************************************/

ERR_VALUE kmer_table_create(const size_t KMerSize, const size_t Size, const KMER_TABLE_CALLBACKS *Callbacks, PKMER_TABLE *Table)
{
	PKMER_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_TABLE(&tmpTable);
	if (ret == ERR_SUCCESS) {
		tmpTable->KHashTable = kh_init(vertexTable);
		if (tmpTable != NULL) {
			tmpTable->KHashTable->Context = KMerSize;
			tmpTable->LastOrder = 0;
			tmpTable->KMerSize = KMerSize;
			if (Callbacks == NULL) {
				tmpTable->Callbacks.Context = NULL;
				tmpTable->Callbacks.OnInsert = _on_insert_dummy_callback;
				tmpTable->Callbacks.OnDelete = _on_delete_dummy_callback;
				tmpTable->Callbacks.OnCopy = _on_copy_dummy_callback;
				tmpTable->Callbacks.OnPrint = _on_print_dummy_callback;
			} else tmpTable->Callbacks = *Callbacks;

			*Table = tmpTable;
			ret = ERR_SUCCESS;
		} else ret = ERR_OUT_OF_MEMORY;
		
		if (ret != ERR_SUCCESS)
			utils_free(tmpTable);
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}


void kmer_table_destroy(PKMER_TABLE Table)
{
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it))
			Table->Callbacks.OnDelete(Table, kh_val(Table->KHashTable, it), Table->Callbacks.Context);
	}

	kh_destroy(vertexTable, Table->KHashTable);
	utils_free(Table);

	return;
}


void kmer_table_print(FILE *Stream, const PKMER_TABLE Table)
{
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it))
			Table->Callbacks.OnPrint(Table, kh_val(Table->KHashTable, it), Stream);
	}

	return;
}


ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const KMER *KMer, void *Data)
{
	int r = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	khiter_t it;

	it = kh_put(vertexTable, Table->KHashTable, KMer, &r);
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
			assert(FALSE);
			break;
	}

	return ret;
}


ERR_VALUE kmer_table_delete(PKMER_TABLE Table, const KMER *KMer)
{
	khiter_t it;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	it = kh_get(vertexTable, Table->KHashTable, KMer);
	if (it != kh_end(Table->KHashTable)) {
		Table->Callbacks.OnDelete(Table, kh_val(Table->KHashTable, it), Table->Callbacks.Context);
		kh_del(vertexTable, Table->KHashTable, it);
		ret = ERR_SUCCESS;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


void *kmer_table_get(const struct _KMER_TABLE *Table, const struct _KMER *KMer)
{
	void *ret = NULL;
	khiter_t it;

	it = kh_get(vertexTable, Table->KHashTable, KMer);
	if (it != kh_end(Table->KHashTable))
		ret = kh_val(Table->KHashTable, it);

	return ret;
}


ERR_VALUE kmer_table_first(const PKMER_TABLE Table, void **Slot, void **Data)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = kh_begin(Table->KHashTable); it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it)) {
			*Slot = (void *)it;
			*Data = kh_val(Table->KHashTable, it);
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}


ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const void *Current, void **Next, void **Data)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NO_MORE_ENTRIES;
	for (khiter_t it = (khiter_t)Current + 1; it != kh_end(Table->KHashTable); ++it) {
		if (kh_exist(Table->KHashTable, it)) {
			*Next = (void *)it;
			*Data = kh_val(Table->KHashTable, it);
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}
