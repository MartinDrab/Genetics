
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"


/************************************************************************/
/*                      HELPER MACROS AND TYPES                         */
/************************************************************************/

#define _kmer_table_entry_empty(aEntry)						((aEntry)->KMer == NULL)
#define _kmer_table_next_entry(aEntry)						((aEntry) + 1)
#define _kmer_table_nth_entry(aTable, aIndex)				((aTable)->Entries + (aIndex))
#define _next_hash_attempt(aHash, aAttempt)		(aHash + 2 * aAttempt + 1)




typedef enum _ETableOpType {
	totSearch,
	totInsert,
	totDelete,
} ETableOpType, *PETableOpType;

/************************************************************************/
/*                   HELPER FUNCTIONS                                   */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_TABLE)
static UTILS_TYPED_CALLOC_FUNCTION(KMER_TABLE_ENTRY)


static PKMER_TABLE_ENTRY _kmer_table_get_slot_insert_hint(const KMER_TABLE *Table, size_t Hash, const KMER *KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;
	PKMER_TABLE_ENTRY firstDeleted = NULL;

	ret = _kmer_table_nth_entry(Table, Hash);
	if ((!_kmer_table_entry_empty(ret) || ret->Deleted) && (_kmer_table_entry_empty(ret) || !kmer_equal(ret->KMer, KMer))) {
		size_t attempt = 1;
		PKMER_TABLE_ENTRY first = ret;

		if (ret->Deleted)
			firstDeleted = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt) % Table->Size;
			ret = _kmer_table_nth_entry(Table,  Hash);
			if (firstDeleted == NULL && ret->Deleted)
				firstDeleted = ret;
			
			if (!ret->Deleted && (_kmer_table_entry_empty(ret) || kmer_equal(ret->KMer, KMer)))
				break;

			if (first == ret) {
				ret = NULL;
				break;
			}

			++attempt;
		} while (TRUE);

		if (firstDeleted != NULL && ret != NULL && _kmer_table_entry_empty(ret))
			ret = firstDeleted;
	}

	return ret;
}


static PKMER_TABLE_ENTRY _kmer_table_get_slot_delsearch_hint(const KMER_TABLE *Table, size_t Hash, const KMER *KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = _kmer_table_nth_entry(Table, Hash);
	if ((!_kmer_table_entry_empty(ret) || ret->Deleted) && (_kmer_table_entry_empty(ret) || !kmer_equal(ret->KMer, KMer))) {
		size_t attempt = 1;
		PKMER_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt) % Table->Size;
			ret = _kmer_table_nth_entry(Table, Hash);
			if (!ret->Deleted && (_kmer_table_entry_empty(ret) || kmer_equal(ret->KMer, KMer)))
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


static PKMER_TABLE_ENTRY _kmer_table_get_slot(const KMER_TABLE *Table, const KMER *KMer, const ETableOpType OpType)
{
	size_t hash = 0;
	PKMER_TABLE_ENTRY ret = NULL;

	hash = kmer_hash(KMer) % Table->Size;
	assert(hash < Table->Size);
	switch (OpType) {
		case totInsert:
			ret = _kmer_table_get_slot_insert_hint(Table, hash, KMer);
			break;
		case totSearch:
		case totDelete:
			ret = _kmer_table_get_slot_delsearch_hint(Table, hash, KMer);
			break;
	}

	return ret;
}


static void _on_insert_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order)
{
	return;
}


static void _on_delete_dummy_callback(struct _KMER_TABLE *Table, void *ItemData)
{
	return;
}


static ERR_VALUE _on_copy_dummy_callback(struct _KMER_TABLE *Table, void *ItemData, void **Copy)
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

	if (utils_is_prime(Size)) {
		ret = utils_malloc_KMER_TABLE(&tmpTable);
		if (ret == ERR_SUCCESS) {
			tmpTable->NumberOfItems = 0;
			tmpTable->LastOrder = 0;
			tmpTable->Size = Size;
			tmpTable->KMerSize = KMerSize;
			if (Callbacks == NULL) {
				tmpTable->Callbacks.OnInsert = _on_insert_dummy_callback;
				tmpTable->Callbacks.OnDelete = _on_delete_dummy_callback;
				tmpTable->Callbacks.OnCopy = _on_copy_dummy_callback;
				tmpTable->Callbacks.OnPrint = _on_print_dummy_callback;
			} else tmpTable->Callbacks = *Callbacks;

			ret = utils_calloc_KMER_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
			if (ret == ERR_SUCCESS) {
				memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_TABLE_ENTRY));
				*Table = tmpTable;
			}

			if (ret != ERR_SUCCESS)
				utils_free(tmpTable);
		}
	} else ret = ERR_NOT_A_PRIME;

	return ret;
}


void kmer_table_destroy(PKMER_TABLE Table)
{
	PKMER_TABLE_ENTRY entry = Table->Entries;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry))
			Table->Callbacks.OnDelete(Table, entry->Data);

		entry = _kmer_table_next_entry(entry);
	}

	utils_free(Table->Entries);
	utils_free(Table);

	return;
}


ERR_VALUE kmer_table_extend(PKMER_TABLE Table)
{
	size_t newSize = 0;
	PKMER_TABLE newTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	newSize = Table->Size * 2;
	newSize = utils_next_prime(newSize);
	ret = kmer_table_create(Table->KMerSize, newSize, &Table->Callbacks, &newTable);
	if (ret == ERR_SUCCESS) {
		PKMER_TABLE_ENTRY entry = Table->Entries;
		PKMER_TABLE_ENTRY newSlot = NULL;

		newTable->LastOrder = Table->LastOrder;
		for (size_t i = 0; i < Table->Size; ++i) {
			if (!_kmer_table_entry_empty(entry)) {				
				newSlot = _kmer_table_get_slot(newTable, entry->KMer, totInsert);
				if (newSlot != NULL)
					memcpy(newSlot, entry, sizeof(KMER_TABLE_ENTRY));
				else ret = ERR_TABLE_FULL;
			}

			if (ret != ERR_SUCCESS)
				break;

			entry = _kmer_table_next_entry(entry);
		}

		if (ret == ERR_SUCCESS) {
			utils_free(Table->Entries);
			newTable->NumberOfItems = Table->NumberOfItems;
			newTable->LastOrder = Table->LastOrder;
			memcpy(Table, newTable, sizeof(KMER_TABLE));
		}

		if (ret != ERR_SUCCESS)
			kmer_table_destroy(newTable);
	}

	return ret;
}


void kmer_table_print(FILE *Stream, const PKMER_TABLE Table)
{
	const KMER_TABLE_ENTRY *entry = Table->Entries;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry))
			Table->Callbacks.OnPrint(Table, entry->Data, Stream);

		entry = _kmer_table_next_entry(entry);
	}

	return;
}


ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const KMER *KMer, void *Data)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Table->NumberOfItems * 100 / Table->Size < 50) {
		entry = _kmer_table_get_slot(Table, KMer, totInsert);
		if (entry != NULL) {
			if (_kmer_table_entry_empty(entry)) {
				entry->Deleted = FALSE;
				entry->Data = Data;
				entry->KMer = KMer;
				Table->Callbacks.OnInsert(Table, entry->Data, Table->LastOrder);
				++Table->LastOrder;
				++Table->NumberOfItems;
				ret = ERR_SUCCESS;
			} else ret = ERR_ALREADY_EXISTS;
		} else ret = ERR_TABLE_FULL;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


ERR_VALUE kmer_table_delete(PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot(Table, KMer, totDelete);
	if (entry != NULL) {
		if (!_kmer_table_entry_empty(entry)) {
			Table->Callbacks.OnDelete(Table, entry->Data);
			memset(entry, 0, sizeof(KMER_TABLE_ENTRY));
			entry->Deleted = TRUE;
			--Table->NumberOfItems;
			ret = ERR_SUCCESS;
		} else ret = ERR_NOT_FOUND;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


void *kmer_table_get(const struct _KMER_TABLE *Table, const struct _KMER *KMer)
{
	void *ret = NULL;
	PKMER_TABLE_ENTRY entry = NULL;

	entry = _kmer_table_get_slot(Table, KMer, totSearch);
	if (entry != NULL && _kmer_table_entry_empty(entry))
		entry = NULL;

	if (entry != NULL)
		ret = entry->Data;

	return ret;
}


ERR_VALUE kmer_table_first(const PKMER_TABLE Table, PKMER_TABLE_ENTRY *Slot)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_TABLE_ENTRY entry = Table->Entries;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry)) {
			*Slot = entry;
			ret = ERR_SUCCESS;
			break;
		}

		entry = _kmer_table_next_entry(entry);
	}

	return ret;
}


ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const PKMER_TABLE_ENTRY Current, PKMER_TABLE_ENTRY *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_TABLE_ENTRY entry = _kmer_table_next_entry(Current);

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = (Current - Table->Entries) + 1; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry)) {
			*Next = entry;
			ret = ERR_SUCCESS;
			break;
		}

		entry = _kmer_table_next_entry(entry);
	}

	return ret;
}


size_t kmer_hash(const struct _KMER *KMer)
{
	size_t hash = 0;
	const size_t kmerSize = kmer_get_size(KMer);

	for (size_t i = 0; i < kmerSize; ++i) {
		hash <<= 1;
		hash += (kmer_get_base(KMer, i));
	}

	return hash;
}


ERR_VALUE kmer_table_get_multiple(const KMER_TABLE *Table, const KMER *KMer, PDYM_ARRAY DataArray)
{
	size_t hash = 0;
	const KMER_TABLE_ENTRY *entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	hash = kmer_hash(KMer) % Table->Size;
	entry = _kmer_table_nth_entry(Table, hash);
	if (!_kmer_table_entry_empty(entry) && kmer_seq_equal(KMer, entry->KMer))
		ret = dym_array_push_back(DataArray, entry->Data);

	if (ret == ERR_SUCCESS && (!_kmer_table_entry_empty(entry) || entry->Deleted)) {
		size_t attempt = 1;
		const KMER_TABLE_ENTRY *first = entry;

		do {
			hash = _next_hash_attempt(hash, attempt) % Table->Size;
			entry = _kmer_table_nth_entry(Table, hash);
			if (_kmer_table_entry_empty(entry) && !entry->Deleted)
				break;
				
			if (first == entry) {
				entry = NULL;
				break;
			}

			if (!_kmer_table_entry_empty(entry) && kmer_seq_equal(KMer, entry->KMer))
				ret = dym_array_push_back(DataArray, entry->Data);

			++attempt;
		} while (ret == ERR_SUCCESS);
	}

	return ret;
}
