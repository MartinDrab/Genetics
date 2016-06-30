
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"


/************************************************************************/
/*                         HELPER MACROS                                */
/************************************************************************/

typedef enum _ETableOpType {
	totSearch,
	totInsert,
	totDelete,
} ETableOpType, *PETableOpType;


#define _kmer_edge_table_entry_empty(aEntry)						((aEntry)->Source == NULL)
#define _next_hash_attempt(aHash, aAttempt, aModulus)				((aHash + 2 * aAttempt + 1) % aModulus)



/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE_TABLE)
static UTILS_TYPED_CALLOC_FUNCTION(KMER_EDGE_TABLE_ENTRY)


static PKMER_EDGE_TABLE_ENTRY _kmer_edge_table_get_slot_insert_hint(const KMER_EDGE_TABLE *Table, size_t Hash, const KMER *Source, const KMER *Dest)
{
	PKMER_EDGE_TABLE_ENTRY firstDeleted = NULL;
	PKMER_EDGE_TABLE_ENTRY ret = NULL;

	ret = Table->Entries + Hash;
	if ((!_kmer_edge_table_entry_empty(ret) || ret->Deleted) && (ret->Deleted || !kmer_equal(ret->Source, Source) || !kmer_equal(ret->Dest, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE_TABLE_ENTRY first = ret;

		if (ret->Deleted)
			firstDeleted = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
			if (firstDeleted == NULL && ret->Deleted)
				firstDeleted = ret;
			
			if (!ret->Deleted && (_kmer_edge_table_entry_empty(ret) || (kmer_equal(ret->Source, Source) && kmer_equal(ret->Dest, Dest))))
				break;

			if (first == ret) {
				ret = NULL;
				break;
			}

			++attempt;
		} while (TRUE);

		if (firstDeleted != NULL && ret != NULL && _kmer_edge_table_entry_empty(ret))
			ret = firstDeleted;
	}

	return ret;
}


static PKMER_EDGE_TABLE_ENTRY _kmer_edge_table_get_slot_delsearch_hint(const KMER_EDGE_TABLE *Table, size_t Hash, const KMER *Source, const KMER *Dest)
{
	PKMER_EDGE_TABLE_ENTRY ret = NULL;

	ret = Table->Entries + Hash;
	if ((!_kmer_edge_table_entry_empty(ret) || ret->Deleted) && (ret->Deleted || !kmer_equal(ret->Source, Source) || !kmer_equal(ret->Dest, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
			if (!ret->Deleted && (_kmer_edge_table_entry_empty(ret) || (kmer_equal(ret->Source, Source) && kmer_equal(ret->Dest, Dest))))
				break;

			if (first == ret) {
				ret = NULL;
				break;
			}

			++attempt;
		} while (TRUE);
	}

	return ret;
}


static PKMER_EDGE_TABLE_ENTRY _kmer_edge_table_get_slot(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest, const ETableOpType OpType)
{
	size_t hash = 0;
	PKMER_EDGE_TABLE_ENTRY ret = NULL;

	hash = kmer_edge_hash(Table, Source, Dest);
	assert(hash < Table->Size);
	switch (OpType) {
		case totInsert:
			ret = _kmer_edge_table_get_slot_insert_hint(Table, hash, Source, Dest);
			break;
		case totDelete:
		case totSearch:
			ret = _kmer_edge_table_get_slot_delsearch_hint(Table, hash, Source, Dest);
			break;
	}

	return ret;
}



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

ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t X, const size_t Size, const PKMER_EDGE_TABLE_CALLBACKS Callbacks, PKMER_EDGE_TABLE *Table)
{
	PKMER_EDGE_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (utils_is_prime(Size)) {
		ret = utils_malloc_KMER_EDGE_TABLE(&tmpTable);
		if (ret == ERR_SUCCESS) {
			tmpTable->Size = Size;
			tmpTable->X = X;
			tmpTable->KMerSize = KMerSize;
			tmpTable->PowX = utils_pow_mod(tmpTable->X, KMerSize - 1, tmpTable->Size);
			if (Callbacks == NULL) {
				tmpTable->Callbacks.OnInsert = _on_insert_dummy_callback;
				tmpTable->Callbacks.OnDelete = _on_delete_dummy_callback;
				tmpTable->Callbacks.OnCopy = _on_copy_dummy_callback;
				tmpTable->Callbacks.OnPrint = _on_print_dummy_callback;
			} else tmpTable->Callbacks = *Callbacks;
			
			tmpTable->LastOrder = 0;
			ret = utils_mul_inverse(X, Size, &tmpTable->Inverse);
			if (ret == ERR_SUCCESS) {
				assert((X*tmpTable->Inverse % Size) == 1);
				ret = utils_calloc_KMER_EDGE_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
				if (ret == ERR_SUCCESS) {
					memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_EDGE_TABLE_ENTRY));
					*Table = tmpTable;
				}
			}

			if (ret != ERR_SUCCESS)
				utils_free(tmpTable);
		}
	} else ret = ERR_NOT_A_PRIME;

	return ret;
}


void kmer_edge_table_destroy(PKMER_EDGE_TABLE Table)
{
	PKMER_EDGE_TABLE_ENTRY entry = Table->Entries;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_edge_table_entry_empty(entry)) {
			Table->Callbacks.OnDelete(Table, entry->Data);
			memset(entry, 0, sizeof(KMER_EDGE_TABLE_ENTRY));
		}

		++entry;
	}
	
	utils_free(Table->Entries);
	utils_free(Table);

	return;
}


ERR_VALUE kmer_edge_table_extend(PKMER_EDGE_TABLE Table)
{
	size_t newSize = 0;
	PKMER_EDGE_TABLE newTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	newSize = Table->Size * 2;
	newSize = utils_next_prime(newSize);
	ret = kmer_edge_table_create(Table->KMerSize, Table->X, newSize, &Table->Callbacks, &newTable);
	if (ret == ERR_SUCCESS) {
		const KMER_EDGE_TABLE_ENTRY *entry = Table->Entries;
		PKMER_EDGE_TABLE_ENTRY newSlot = NULL;

		for (size_t i = 0; i < Table->Size; ++i) {
			if (!_kmer_edge_table_entry_empty(entry)) {
				newSlot = _kmer_edge_table_get_slot(newTable, entry->Source, entry->Dest, totInsert);
				if (newSlot != NULL)
					memcpy(newSlot, entry, sizeof(KMER_EDGE_TABLE_ENTRY));
				else ret = ERR_TABLE_FULL;
			}

			if (ret != ERR_SUCCESS)
				break;

			++entry;
		}

		if (ret == ERR_SUCCESS) {
			utils_free(Table->Entries);
			memcpy(Table, newTable, sizeof(KMER_EDGE_TABLE));
		}

		if (ret != ERR_SUCCESS)
			kmer_edge_table_destroy(newTable);
	}

	return ret;
}


PKMER_EDGE_TABLE_ENTRY kmer_edge_table_get(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest)
{
	PKMER_EDGE_TABLE_ENTRY ret = NULL;

	ret = _kmer_edge_table_get_slot(Table, Source, Dest, totSearch);
	if (ret != NULL && _kmer_edge_table_entry_empty(ret))
		ret = NULL;

	return ret;
}


void *kmer_edge_table_get_data(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest)
{
	void *ret = NULL;
	PKMER_EDGE_TABLE_ENTRY entry = kmer_edge_table_get(Table, Source, Dest);

	if (entry != NULL)
		ret = entry->Data;

	return ret;
}


void kmer_edge_table_delete_by_entry(PKMER_EDGE_TABLE Table, PKMER_EDGE_TABLE_ENTRY Entry)
{
	Table->Callbacks.OnDelete(Table, Entry->Data);
	memset(Entry, 0, sizeof(KMER_EDGE_TABLE_ENTRY));
	Entry->Deleted = TRUE;

	return;
}


ERR_VALUE kmer_edge_table_delete(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_edge_table_get_slot(Table, Source, Dest, totDelete);
	if (entry != NULL && !_kmer_edge_table_entry_empty(entry)) {
		kmer_edge_table_delete_by_entry(Table, entry);
		ret = ERR_SUCCESS;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const KMER *Source, const KMER *Dest, void *Data)
{
	PKMER_EDGE_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_edge_table_get_slot(Table, Source, Dest, totInsert);
	if (entry != NULL) {
		if (_kmer_edge_table_entry_empty(entry) || entry->Deleted) {
			entry->Deleted = FALSE;
			entry->Source = Source;
			entry->Dest = Dest;
			entry->Data = Data;
			Table->Callbacks.OnInsert(Table, Data, Table->LastOrder);
			Table->LastOrder++;
			ret = ERR_SUCCESS;
		} else ret = ERR_ALREADY_EXISTS;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


size_t kmer_edge_hash(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest)
{
	size_t hash = 0;

	for (size_t i = 0; i < kmer_get_size(Source); ++i) {
		hash = (hash*Table->X + kmer_get_base(Source, i));
		hash %= Table->Size;
	}

	for (size_t i = 0; i < kmer_get_size(Dest); ++i) {
		hash = (hash*Table->X + kmer_get_base(Dest, i));
		hash %= Table->Size;
	}

	return hash;
}


void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table)
{
	PKMER_EDGE_TABLE_ENTRY edge = NULL;

	edge = Table->Entries;
	for (size_t j = 0; j < Table->Size; ++j) {
		if (!_kmer_edge_table_entry_empty(edge))
			Table->Callbacks.OnPrint(Table, edge->Data, Stream);

		++edge;
	}

	return;
}


ERR_VALUE kmer_edge_table_first(const PKMER_EDGE_TABLE Table, PKMER_EDGE_TABLE_ENTRY *Slot)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE_ENTRY entry = Table->Entries;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_edge_table_entry_empty(entry)) {
			*Slot = entry;
			ret = ERR_SUCCESS;
			break;
		}

		++entry;
	}

	return ret;
}


ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const PKMER_EDGE_TABLE_ENTRY Current, PKMER_EDGE_TABLE_ENTRY *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE_ENTRY entry = Current + 1;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = (Current - Table->Entries) + 1; i < Table->Size; ++i) {
		if (!_kmer_edge_table_entry_empty(entry)) {
			*Next = entry;
			ret = ERR_SUCCESS;
			break;
		}

		++entry;
	}

	return ret;
}
