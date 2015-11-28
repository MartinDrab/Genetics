
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"


/************************************************************************/
/*                      HELPER MACROS AND TYPES                         */
/************************************************************************/

#define _kmer_table_entry_empty(aEntry)						((aEntry)->KMer == NULL)
#define _next_hash_attempt(aHash, aAttempt, aModulus)		((aHash + 2 * aAttempt + 1) % aModulus)

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


static PKMER_TABLE_ENTRY _kmer_table_get_slot_insert_hint(const PKMER_TABLE Table, size_t Hash, const PKMER KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = Table->Entries + Hash;
	if (!_kmer_table_entry_empty(ret) && !kmer_equal(ret->KMer, KMer)) {
		size_t attempt = 1;
		PKMER_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
			if (_kmer_table_entry_empty(ret) || kmer_equal(ret->KMer, KMer))
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


static PKMER_TABLE_ENTRY _kmer_table_get_slot_delsearch_hint(const PKMER_TABLE Table, size_t Hash, const PKMER KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = Table->Entries + Hash;
	if ((!_kmer_table_entry_empty(ret) || ret->Deleted) && (_kmer_table_entry_empty(ret) || !kmer_equal(ret->KMer, KMer))) {
		size_t attempt = 1;
		PKMER_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
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


static PKMER_TABLE_ENTRY _kmer_table_get_slot(const PKMER_TABLE Table, const PKMER KMer, const ETableOpType OpType)
{
	size_t hash = 0;
	PKMER_TABLE_ENTRY ret = NULL;

	hash = kmer_hash(Table, KMer);
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


/************************************************************************/
/*                  PUBLIC FUNCTIONS                                    */
/************************************************************************/

ERR_VALUE kmer_table_create(const size_t KMerSize, const size_t X, const size_t Size, PKMER_TABLE *Table)
{
	PKMER_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (utils_is_prime(Size)) {
		ret = utils_malloc_KMER_TABLE(&tmpTable);
		if (ret == ERR_SUCCESS) {
			tmpTable->LastOrder = 0;
			tmpTable->Size = Size;
			tmpTable->X = X;
			tmpTable->KMerSize = KMerSize;
			tmpTable->PowX = utils_pow_mod(tmpTable->X, KMerSize - 1, tmpTable->Size);
			ret = utils_mul_inverse(X, Size, &tmpTable->Inverse);
			if (ret == ERR_SUCCESS) {
				assert((X*tmpTable->Inverse % Size) == 1);
				ret = utils_calloc_KMER_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
				if (ret == ERR_SUCCESS) {
					memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_TABLE_ENTRY));
					*Table = tmpTable;
				}
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
			kmer_free(entry->KMer);

		++entry;
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
	ret = kmer_table_create(Table->KMerSize, Table->X, newSize, &newTable);
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

			++entry;
		}

		if (ret == ERR_SUCCESS) {
			utils_free(Table->Entries);
			memcpy(Table, newTable, sizeof(KMER_TABLE));
		}

		if (ret != ERR_SUCCESS)
			kmer_table_destroy(newTable);
	}

	return ret;
}


ERR_VALUE kmer_table_copy(const PKMER_TABLE Source, PKMER_TABLE * Copied)
{
	PKMER_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_TABLE(&tmpTable);
	if (ret == ERR_SUCCESS) {
		tmpTable->Inverse = Source->Inverse;
		tmpTable->KMerSize = Source->KMerSize;
		tmpTable->PowX = Source->PowX;
		tmpTable->Size = Source->Size;
		tmpTable->LastOrder = Source->LastOrder;
		ret = utils_calloc_KMER_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
		if (ret == ERR_SUCCESS) {
			PKMER_TABLE_ENTRY sourceEntry = Source->Entries;
			PKMER_TABLE_ENTRY destEntry = tmpTable->Entries;

			memset(tmpTable->Entries, 0, sizeof(KMER_TABLE_ENTRY)*tmpTable->Size);
			for (size_t i = 0; i < tmpTable->Size; ++i) {
				if (!_kmer_table_entry_empty(sourceEntry) || sourceEntry->Deleted) {
					memcpy(destEntry, sourceEntry, sizeof(KMER_TABLE_ENTRY));
					if (sourceEntry->KMer != NULL) {
						destEntry->KMer = kmer_copy(sourceEntry->KMer);
						if (destEntry->KMer == NULL)
							ret = ERR_OUT_OF_MEMORY;
					}
				}

				if (ret != ERR_SUCCESS) {
					--destEntry;
					for (size_t j = 0; j < i; ++j) {
						if (!_kmer_table_entry_empty(destEntry))
							kmer_free(destEntry->KMer);

						--destEntry;
					}
				}

				++sourceEntry;
				++destEntry;
			}

			if (ret == ERR_SUCCESS)
				*Copied = tmpTable;
				
			if (ret != ERR_SUCCESS)
				utils_free(tmpTable->Entries);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpTable);
	}

	return ret;
}


void kmer_table_print(FILE *Stream, const PKMER_TABLE Table)
{
	PKMER_TABLE_ENTRY entry = NULL;

	entry = Table->Entries;
	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry)) {
			fprintf(Stream, "\t");
			kmer_print(Stream, entry->KMer);
			fprintf(Stream, "[label=\"");
			kmer_print(Stream, entry->KMer);
			fprintf(Stream, "\\nIN=%u; OUT=%u\"]", entry->DegreeIn, entry->degreeOut);
			fprintf(Stream, ";\n");
		}

		++entry;
	}

	return;
}


ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot(Table, KMer, totInsert);
	if (entry != NULL) {
		if (_kmer_table_entry_empty(entry)) {
			entry->Deleted = FALSE;
			entry->DegreeIn = 0;
			entry->degreeOut = 0;
			entry->KMer = kmer_copy(KMer);
			entry->Order = Table->LastOrder;
			++Table->LastOrder;
			memset(&entry->AdvancedInfo, 0, sizeof(entry->AdvancedInfo));
			ret = (!_kmer_table_entry_empty(entry)) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		} else ret = ERR_ALREADY_EXISTS;
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
			kmer_free(entry->KMer);
			memset(entry, 0, sizeof(KMER_TABLE_ENTRY));
			entry->Deleted = TRUE;
			ret = ERR_SUCCESS;
		} else ret = ERR_NOT_FOUND;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_table_insert_hint(PKMER_TABLE Table, const PKMER KMer, const size_t Hash)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot_insert_hint(Table, Hash, KMer);
	if (entry != NULL) {
		if (_kmer_table_entry_empty(entry)) {
			entry->Deleted = FALSE;
			entry->DegreeIn = 0;
			entry->degreeOut = 0;
			entry->Order = Table->LastOrder;
			++Table->LastOrder;
			entry->KMer = kmer_copy(KMer);
			ret = (_kmer_table_entry_empty(entry)) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		} else ret = ERR_ALREADY_EXISTS;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


PKMER_TABLE_ENTRY kmer_table_get(const PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = _kmer_table_get_slot(Table, KMer, totSearch);
	if (ret != NULL && _kmer_table_entry_empty(ret))
		ret = NULL;

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

		++entry;
	}

	return ret;
}


ERR_VALUE kmer_table_next(const PKMER_TABLE Table, const PKMER_TABLE_ENTRY Current, PKMER_TABLE_ENTRY *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_TABLE_ENTRY entry = Current + 1;

	ret = ERR_NO_MORE_ENTRIES;
	for (size_t i = (Current - Table->Entries) + 1; i < Table->Size; ++i) {
		if (!_kmer_table_entry_empty(entry)) {
			*Next = entry;
			ret = ERR_SUCCESS;
			break;
		}

		++entry;
	}

	return ret;
}


size_t kmer_hash(const PKMER_TABLE Table, const PKMER KMer)
{
	size_t hash = 0;

	for (size_t i = 0; i < kmer_get_size(KMer); ++i) {
		hash = (hash*Table->X + kmer_get_base(KMer, i));
		hash %= Table->Size;
	}

	return hash;
}


size_t kmer_hash_advance(const PKMER_TABLE Table, const PKMER KMer, const size_t Hash, const char NewBase)
{
	return ((Hash + Table->Size - kmer_get_base(KMer, 0))*Table->Inverse + NewBase*Table->PowX);
}
