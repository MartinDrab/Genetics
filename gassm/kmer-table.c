
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"


/************************************************************************/
/*                      HELPER MACROS                                   */
/************************************************************************/

#define _kmer_table_entry_empty(aEntry)						((aEntry)->KMer == NULL)
#define _next_hash_attempt(aHash, aAttempt, aModulus)		((aHash + 2 * aAttempt + 1) % aModulus)

/************************************************************************/
/*                   HELPER FUNCTIONS                                   */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_TABLE)
static UTILS_TYPED_CALLOC_FUNCTION(KMER_TABLE_ENTRY)


static PKMER_TABLE_ENTRY _kmer_table_get_slot_hint(const PKMER_TABLE Table, size_t Hash, const PKMER KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = Table->Entries + Hash;
	if (!_kmer_table_entry_empty(ret) && (KMer == NULL || !kmer_equal(ret->KMer, KMer))) {
		size_t attempt = 1;
		PKMER_TABLE_ENTRY first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
			if (_kmer_table_entry_empty(ret) || (KMer != NULL && kmer_equal(ret->KMer, KMer)))
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

static PKMER_TABLE_ENTRY _kmer_table_get_slot(const PKMER_TABLE Table, const PKMER KMer)
{
	size_t hash = 0;
	PKMER_TABLE_ENTRY ret = NULL;

	hash = kmer_hash(Table, KMer);
	assert(hash < Table->Size);
	ret = _kmer_table_get_slot_hint(Table, hash, KMer);

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

		for (size_t i = 0; i < Table->Size; ++i) {
			if (!_kmer_table_entry_empty(entry)) {				
				newSlot = _kmer_table_get_slot(newTable, entry->KMer);
				if (newSlot != NULL)
					memcpy(newSlot, entry, sizeof(KMER_TABLE_ENTRY));
				else ret = ERR_TABLE_FULL;
			}

			if (ret != ERR_SUCCESS)
				break;

			++entry;
		}

		if (ret == ERR_SUCCESS) {
			free(Table->Entries);
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
		ret = utils_calloc_KMER_TABLE_ENTRY(tmpTable->Size, &tmpTable->Entries);
		if (ret == ERR_SUCCESS) {
			PKMER_TABLE_ENTRY sourceEntry = Source->Entries;
			PKMER_TABLE_ENTRY destEntry = tmpTable->Entries;

			memset(tmpTable->Entries, 0, sizeof(KMER_TABLE_ENTRY)*tmpTable->Size);
			for (size_t i = 0; i < tmpTable->Size; ++i) {
				if (!_kmer_table_entry_empty(sourceEntry)) {
					memcpy(destEntry, sourceEntry, sizeof(KMER_TABLE_ENTRY));
					destEntry->KMer = kmer_copy(sourceEntry->KMer);
					if (destEntry->KMer == NULL)
						ret = ERR_OUT_OF_MEMORY;
				}

				if (ret != ERR_SUCCESS) {
					--destEntry;
					for (size_t j = 0; j < i; ++j) {
						if (_kmer_table_entry_empty(destEntry))
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


void kmer_table_print(const PKMER_TABLE Table)
{
	PKMER_TABLE_ENTRY entry = NULL;

	entry = Table->Entries;
	for (size_t i = 0; i < Table->Size; ++i) {
		if (entry->KMer != NULL) {
			printf("\t");
			kmer_print(entry->KMer);
			printf("[label=\"");
			kmer_print(entry->KMer);
			printf("\\nIN=%u; OUT=%u\"]", entry->DegreeIn, entry->degreeOut);
			printf(";\n");
		}

		++entry;
	}

	return;
}


ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot(Table, KMer);
	if (entry != NULL) {
		if (_kmer_table_entry_empty(entry)) {
			entry->DegreeIn = 0;
			entry->degreeOut = 0;
			entry->KMer = kmer_copy(KMer);
			ret = (!_kmer_table_entry_empty(entry)) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		} else ret = ERR_ALREADY_EXISTS;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


ERR_VALUE kmer_table_insert_hint(PKMER_TABLE Table, const PKMER KMer, const size_t Hash)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot_hint(Table, Hash, KMer);
	if (entry != NULL) {
		if (_kmer_table_entry_empty(entry)) {
			entry->DegreeIn = 0;
			entry->degreeOut = 0;
			entry->KMer = kmer_copy(KMer);
			ret = (_kmer_table_entry_empty(entry)) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		} else ret = ERR_ALREADY_EXISTS;
	} else ret = ERR_TABLE_FULL;

	return ret;
}


PKMER_TABLE_ENTRY kmer_table_get(const PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY ret = NULL;

	ret = _kmer_table_get_slot(Table, KMer);
	if (ret != NULL && _kmer_table_entry_empty(ret))
		ret = NULL;

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
