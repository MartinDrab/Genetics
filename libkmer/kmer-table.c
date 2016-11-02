
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

#define _kmer_table_slot_used(aFlags, aIndex)				(((aFlags)[(aIndex)]) & KMER_TABLE_ENTRY_FLAG_USED)
#define _kmer_table_slot_deleted(aFlags, aIndex)			(((aFlags)[(aIndex)]) & KMER_TABLE_ENTRY_FLAG_DELETED)
#define _next_hash_attempt(aHash, aAttempt)		(aHash + 2 * aAttempt + 1)
#define _kmer_table_nth_key(aKeys, aIndex)		((aKeys)[(aIndex)])


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


static size_t _kmer_table_get_slot_insert_hint(const KMER_TABLE *Table, size_t Hash, const KMER *KMer)
{
	size_t ret = (size_t)-1;
	size_t firstDeleted = (size_t)-1;
	const uint8_t *flags = Table->Flags;
	const KMER **keys = Table->Keys;
	const size_t tableSize = Table->Size;

	if (flags[Hash] != 0 && (!_kmer_table_slot_used(flags, Hash) || !kmer_equal(KMer, _kmer_table_nth_key(keys, Hash)))) {
		if (_kmer_table_slot_deleted(flags, Hash))
			firstDeleted = Hash;

		for (size_t attempt = 1; attempt <= tableSize; ++attempt) {
			Hash = _next_hash_attempt(Hash, attempt) % tableSize;
			if (firstDeleted == (size_t)-1 && _kmer_table_slot_deleted(flags, Hash))
				firstDeleted = Hash;
			
			if (!_kmer_table_slot_deleted(flags, Hash) && (!_kmer_table_slot_used(flags, Hash) || kmer_equal(_kmer_table_nth_key(keys, Hash), KMer))) {
				ret = Hash;
				break;
			}
		}

		if (firstDeleted != (size_t)-1 && ret != (size_t)-1 && !_kmer_table_slot_used(flags, Hash))
			ret = firstDeleted;
	} else ret = Hash;

	return ret;
}


static size_t _kmer_table_get_slot_delsearch_hint(const KMER_TABLE *Table, size_t Hash, const KMER *KMer)
{
	size_t ret = (size_t)-1;
	const uint8_t *flags = Table->Flags;
	const KMER **keys = Table->Keys;
	const size_t tableSize = Table->Size;

	if ((flags[Hash] != 0) && (!_kmer_table_slot_used(flags, Hash) || !kmer_equal(_kmer_table_nth_key(keys, Hash), KMer))) {
		for (size_t attempt = 1; attempt <= tableSize; ++attempt) {
			Hash = _next_hash_attempt(Hash, attempt) % tableSize;
			if (!_kmer_table_slot_deleted(flags, Hash) && (!_kmer_table_slot_used(flags, Hash) || kmer_equal(_kmer_table_nth_key(keys , Hash), KMer))) {
				ret = Hash;
				break;
			}
		}
	} else ret = Hash;

	return ret;
}


static size_t _kmer_table_get_slot(const KMER_TABLE *Table, const KMER *KMer, const ETableOpType OpType)
{
	size_t hash = 0;

	hash = kmer_hash(KMer) % Table->Size;
	assert(hash < Table->Size);
	switch (OpType) {
		case totInsert:
			hash = _kmer_table_get_slot_insert_hint(Table, hash, KMer);
			break;
		case totSearch:
		case totDelete:
			hash = _kmer_table_get_slot_delsearch_hint(Table, hash, KMer);
			break;
	}

	return hash;
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

			tmpTable->Entries = (PKMER_TABLE_ENTRY)malloc(tmpTable->Size*(sizeof(KMER_TABLE_ENTRY) + sizeof(PKMER) + sizeof(uint8_t)));
			tmpTable->Keys = (PKMER *)(tmpTable->Entries + tmpTable->Size);
			tmpTable->Flags = (uint8_t *)(tmpTable->Keys + tmpTable->Size);
			if (tmpTable->Entries != NULL) {
				memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_TABLE_ENTRY));
				memset(tmpTable->Flags, 0, sizeof(uint8_t)*tmpTable->Size);
				*Table = tmpTable;
			} else ret = ERR_OUT_OF_MEMORY;

			if (ret != ERR_SUCCESS)
				utils_free(tmpTable);
		}
	} else ret = ERR_NOT_A_PRIME;

	return ret;
}


void kmer_table_destroy(PKMER_TABLE Table)
{
	const uint8_t *flags = Table->Flags;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (_kmer_table_slot_used(flags, i))
			Table->Callbacks.OnDelete(Table, Table->Entries[i].Data);
	}

	free(Table->Entries);
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
		const uint8_t *flags = Table->Flags;
		const KMER **keys = Table->Keys;
		const KMER_TABLE_ENTRY *entry = Table->Entries;
		KMER_TABLE_ENTRY *newEntries = newTable->Entries;
		size_t slotIndex = 0;
		
		for (size_t i = 0; i < Table->Size; ++i) {
			if (_kmer_table_slot_used(flags, i)) {				
				slotIndex = _kmer_table_get_slot(newTable, _kmer_table_nth_key(keys, i), totInsert);
				if (slotIndex != (size_t)-1) {
					memcpy(newEntries + slotIndex, entry, sizeof(KMER_TABLE_ENTRY));
					newTable->Flags[slotIndex] = KMER_TABLE_ENTRY_FLAG_USED;
					newTable->Keys[slotIndex] = keys[i];
				} else ret = ERR_TABLE_FULL;
			}

			if (ret != ERR_SUCCESS)
				break;

			++entry;
		}

		if (ret == ERR_SUCCESS) {
			free(Table->Entries);
			newTable->NumberOfItems = Table->NumberOfItems;
			newTable->LastOrder = Table->LastOrder;
			memcpy(Table, newTable, sizeof(KMER_TABLE));
			utils_free(newTable);
		}

		if (ret != ERR_SUCCESS)
			kmer_table_destroy(newTable);
	}

	return ret;
}


void kmer_table_print(FILE *Stream, const PKMER_TABLE Table)
{
	const uint8_t *flags = Table->Flags;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (_kmer_table_slot_used(flags, i))
			Table->Callbacks.OnPrint(Table, Table->Entries[i].Data, Stream);
	}

	return;
}


ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const KMER *KMer, void *Data)
{
	size_t slotIndex = (size_t)-1;
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Table->NumberOfItems * 100 / Table->Size < 50) {
		slotIndex = _kmer_table_get_slot(Table, KMer, totInsert);
		if (slotIndex != (size_t)-1) {
			entry = Table->Entries + slotIndex;
			if (!_kmer_table_slot_used(Table->Flags, slotIndex)) {
				entry->Data = Data;
				Table->Keys[slotIndex] = KMer;
				Table->Flags[slotIndex] = KMER_TABLE_ENTRY_FLAG_USED;
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
	size_t slot = (size_t)-1;
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	slot = _kmer_table_get_slot(Table, KMer, totDelete);
	if (slot != (size_t)-1) {
		entry = Table->Entries + slot;
		if (_kmer_table_slot_used(Table->Flags, slot)) {
			Table->Callbacks.OnDelete(Table, entry->Data);
			memset(entry, 0, sizeof(KMER_TABLE_ENTRY));
			Table->Flags[slot] = KMER_TABLE_ENTRY_FLAG_DELETED;
			--Table->NumberOfItems;
			ret = ERR_SUCCESS;
		} else ret = ERR_NOT_FOUND;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


void *kmer_table_get(const struct _KMER_TABLE *Table, const struct _KMER *KMer)
{
	void *ret = NULL;
	size_t slot = (size_t)-1;

	slot = _kmer_table_get_slot(Table, KMer, totSearch);
	if (slot != (size_t)-1 && !_kmer_table_slot_used(Table->Flags, slot))
		slot = (size_t)-1;

	if (slot != (size_t)-1)
		ret = Table->Entries[slot].Data;

	return ret;
}


ERR_VALUE kmer_table_first(const PKMER_TABLE Table, PKMER_TABLE_ENTRY *Slot)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const uint8_t *flags = Table->Flags;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = 0; i < Table->Size; ++i) {
		if (_kmer_table_slot_used(flags, i)) {
			*Slot = Table->Entries + i;
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}


ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const PKMER_TABLE_ENTRY Current, PKMER_TABLE_ENTRY *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const uint8_t *flags = Table->Flags;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = (Current - Table->Entries) + 1; i < Table->Size; ++i) {
		if (_kmer_table_slot_used(flags, i)) {
			*Next = Table->Entries + i;
			ret = ERR_SUCCESS;
			break;
		}
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
	const uint8_t *flags = Table->Flags;
	const KMER **keys = Table->Keys;
	const size_t tableSize = Table->Size;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	hash = kmer_hash(KMer) % tableSize;
	if (_kmer_table_slot_used(flags, hash)) {
		if (kmer_seq_equal(KMer, _kmer_table_nth_key(keys, hash) ))
			ret = dym_array_push_back(DataArray, (Table->Entries + hash)->Data);
	}

	if (ret == ERR_SUCCESS && flags[hash] != 0) {
		for (size_t attempt = 1; attempt <= tableSize; ++attempt) {
			hash = _next_hash_attempt(hash, attempt) % tableSize;
			if (flags[hash] == 0)
				break;

			if (_kmer_table_slot_used(flags, hash)) {
				if (kmer_seq_equal(KMer, _kmer_table_nth_key(keys, hash))) {
					ret = dym_array_push_back(DataArray, (Table->Entries + hash)->Data);
					if (ret != ERR_SUCCESS)
						break;
				}
			}
		}
	}

	return ret;
}
