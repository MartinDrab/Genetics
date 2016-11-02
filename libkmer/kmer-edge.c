
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
/*                         HELPER MACROS                                */
/************************************************************************/

typedef enum _ETableOpType {
	totSearch,
	totInsert,
	totDelete,
} ETableOpType, *PETableOpType;


#define _kmer_edge_table_entry_empty(aEntry)						((aEntry)->Source == NULL)
#define _kmer_edge_table_next_entry(aEntry)							((aEntry) + 1)
#define _kmer_edge_table_nth_entry(aTable, aIndex)				((aTable)->Entries + (aIndex))
#define _next_hash_attempt(aHash, aAttempt, aModulus)				((aHash + 2 * aAttempt + 1) % aModulus)



/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE_TABLE)
static UTILS_TYPED_CALLOC_FUNCTION(KMER_EDGE_TABLE_ENTRY)



static boolean _kmer_edge_entry_equal(const KMER_EDGE_TABLE_ENTRY *Entry, const KMER *Source, const KMER *Dest)
{
	return (
		kmer_get_number(Source) == kmer_get_number(Entry->Source) &&
		kmer_get_number(Dest) == kmer_get_number(Entry->Dest) &&
		kmer_seq_equal(Source, Entry->Source) &&
		kmer_seq_equal(Dest, Entry->Dest)
		);
}

static boolean _kmer_edge_key_equal(const KMER_EDGE_TABLE_KEY Key1, const KMER_EDGE_TABLE_KEY Key2)
{
	return (
		kmer_get_number(Key1.Source) == kmer_get_number(Key2.Source) &&
		kmer_get_number(Key1.Dest) == kmer_get_number(Key2.Dest) &&
		kmer_seq_equal(Key1.Source, Key2.Source) &&
		kmer_seq_equal(Key1.Dest, Key2.Dest)
		);
}

static size_t _kmer_edge_key_hash(const KMER_EDGE_TABLE_KEY Key)
{
	size_t hash1 = 0;
	size_t hash2 = 0;
	const size_t kmerSize = Key.Source->Size;

	for (size_t i = 0; i < kmerSize; ++i) {
		hash1 <<= 1;
		hash2 <<= 1;
		hash1 += (kmer_get_base(Key.Source, i));
		hash2 += (kmer_get_base(Key.Dest, i));
	}

	return (((Key.Source->Number + 1)*hash1 << 1) + (Key.Dest->Number + 1)*hash2);
}

__KHASH_TYPE(edgeTable, KMER_EDGE_TABLE_KEY, void *)

__KHASH_IMPL(edgeTable, INLINE_FUNCTION, KMER_EDGE_TABLE_KEY, void *, TRUE, _kmer_edge_key_hash, _kmer_edge_key_equal)


static PKMER_EDGE_TABLE_ENTRY _kmer_edge_table_get_slot_insert_hint(const KMER_EDGE_TABLE *Table, size_t Hash, const KMER *Source, const KMER *Dest)
{
	PKMER_EDGE_TABLE_ENTRY firstDeleted = NULL;
	PKMER_EDGE_TABLE_ENTRY ret = NULL;

	ret = _kmer_edge_table_nth_entry(Table, Hash);
	if ((!_kmer_edge_table_entry_empty(ret) || ret->Deleted) && (ret->Deleted || !_kmer_edge_entry_equal(ret, Source, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE_TABLE_ENTRY first = ret;

		if (ret->Deleted)
			firstDeleted = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = _kmer_edge_table_nth_entry(Table, Hash);
			if (firstDeleted == NULL && ret->Deleted)
				firstDeleted = ret;
			
			if (!ret->Deleted && (_kmer_edge_table_entry_empty(ret) || _kmer_edge_entry_equal(ret, Source, Dest)))
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

	ret = _kmer_edge_table_nth_entry(Table, Hash);
	if ((!_kmer_edge_table_entry_empty(ret) || ret->Deleted) && (ret->Deleted || !_kmer_edge_entry_equal(ret, Source, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = _kmer_edge_table_nth_entry(Table, Hash);
			if (!ret->Deleted && (_kmer_edge_table_entry_empty(ret) || _kmer_edge_entry_equal(ret, Source, Dest)))
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

	hash = kmer_edge_hash(Source, Dest) % Table->Size;
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

ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t Size, const PKMER_EDGE_TABLE_CALLBACKS Callbacks, PKMER_EDGE_TABLE *Table)
{
	PKMER_EDGE_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (utils_is_prime(Size)) {
		ret = utils_malloc_KMER_EDGE_TABLE(&tmpTable);
		if (ret == ERR_SUCCESS) {
			tmpTable->NumberOfItems = 0;
			tmpTable->Size = Size;
			tmpTable->KMerSize = KMerSize;
			if (Callbacks == NULL) {
				tmpTable->Callbacks.OnInsert = _on_insert_dummy_callback;
				tmpTable->Callbacks.OnDelete = _on_delete_dummy_callback;
				tmpTable->Callbacks.OnCopy = _on_copy_dummy_callback;
				tmpTable->Callbacks.OnPrint = _on_print_dummy_callback;
			} else tmpTable->Callbacks = *Callbacks;
			
			tmpTable->LastOrder = 0;
			ret = utils_calloc_KMER_EDGE_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
			if (ret == ERR_SUCCESS) {
				memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_EDGE_TABLE_ENTRY));
				*Table = tmpTable;
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

		entry = _kmer_edge_table_next_entry(entry);
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
	ret = kmer_edge_table_create(Table->KMerSize, newSize, &Table->Callbacks, &newTable);
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

			entry = _kmer_edge_table_next_entry(entry);
		}

		if (ret == ERR_SUCCESS) {
			utils_free(Table->Entries);
			newTable->LastOrder = Table->LastOrder;
			newTable->NumberOfItems = Table->NumberOfItems;
			memcpy(Table, newTable, sizeof(KMER_EDGE_TABLE));
			utils_free(newTable);
		}

		if (ret != ERR_SUCCESS)
			kmer_edge_table_destroy(newTable);
	}

	return ret;
}


void *kmer_edge_table_get(const struct _KMER_EDGE_TABLE *Table, const struct _KMER *Source, const struct _KMER *Dest)
{
	void *ret = NULL;
	PKMER_EDGE_TABLE_ENTRY entry = NULL;

	entry = _kmer_edge_table_get_slot(Table, Source, Dest, totSearch);
	if (entry != NULL && !_kmer_edge_table_entry_empty(entry))
		ret = entry->Data;

	return ret;
}


void kmer_edge_table_delete_by_entry(PKMER_EDGE_TABLE Table, PKMER_EDGE_TABLE_ENTRY Entry)
{
	Table->Callbacks.OnDelete(Table, Entry->Data);
	memset(Entry, 0, sizeof(KMER_EDGE_TABLE_ENTRY));
	Entry->Deleted = TRUE;
	--Table->NumberOfItems;

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

	if (Table->NumberOfItems * 100 / Table->Size < 50) {
		entry = _kmer_edge_table_get_slot(Table, Source, Dest, totInsert);
		if (entry != NULL) {
			if (_kmer_edge_table_entry_empty(entry) || entry->Deleted) {
				entry->Deleted = FALSE;
				entry->Source = Source;
				entry->Dest = Dest;
				entry->Data = Data;
				Table->Callbacks.OnInsert(Table, Data, Table->LastOrder);
				++Table->LastOrder;
				++Table->NumberOfItems;
				ret = ERR_SUCCESS;
			} else ret = ERR_ALREADY_EXISTS;
		} else ret = ERR_TABLE_FULL;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


size_t kmer_edge_hash(const struct _KMER *Source, const struct _KMER *Dest)
{
	size_t hash1 = 0;
	size_t hash2 = 0;
	const size_t kmerSize = Source->Size;

	assert(Source->Size == Dest->Size);
	for (size_t i = 0; i < kmerSize; ++i) {
		hash1 <<= 1;
		hash2 <<= 1;
		hash1 += (kmer_get_base(Source, i));
		hash2 += (kmer_get_base(Dest, i));
	}

	return (((Source->Number + 1)*hash1 << 1) + (Dest->Number + 1)*hash2);
}


void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table)
{
	const KMER_EDGE_TABLE_ENTRY *edge = NULL;

	edge = Table->Entries;
	for (size_t j = 0; j < Table->Size; ++j) {
		if (!_kmer_edge_table_entry_empty(edge))
			Table->Callbacks.OnPrint(Table, edge->Data, Stream);

		edge = _kmer_edge_table_next_entry(edge);
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

		entry = _kmer_edge_table_next_entry(entry);
	}

	return ret;
}


ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const PKMER_EDGE_TABLE_ENTRY Current, PKMER_EDGE_TABLE_ENTRY *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE_ENTRY entry = _kmer_edge_table_next_entry(Current);

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = (Current - Table->Entries) + 1; i < Table->Size; ++i) {
		if (!_kmer_edge_table_entry_empty(entry)) {
			*Next = entry;
			ret = ERR_SUCCESS;
			break;
		}

		entry = _kmer_edge_table_next_entry(entry);
	}

	return ret;
}
