
#include <stdlib.h>
#include <assert.h>
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
		tmpTable = (PKMER_TABLE)malloc(sizeof(KMER_TABLE));
		if (tmpTable != NULL) {
			tmpTable->Size = Size;
			tmpTable->X = X;
			tmpTable->KMerSize = KMerSize;
			tmpTable->PowX = utils_pow_mod(tmpTable->X, KMerSize - 1, tmpTable->Size);
			ret = utils_mul_inverse(X, Size, &tmpTable->Inverse);
			if (ret == ERR_SUCCESS) {
				assert((X*tmpTable->Inverse % Size) == 1);
				tmpTable->Entries = (PKMER_TABLE_ENTRY)calloc(tmpTable->Size, sizeof(KMER_TABLE_ENTRY));
				if (tmpTable->Entries != NULL) {
					memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_TABLE_ENTRY));
					*Table = tmpTable;
					ret = ERR_SUCCESS;
				} else ret = ERR_OUT_OF_MEMORY;
			}

			if (ret != ERR_SUCCESS)
				free(tmpTable);
		} else ret = ERR_OUT_OF_MEMORY;
	} else ret = ERR_NOT_A_PRIME;

	return ret;
}

void kmer_table_destroy(PKMER_TABLE Table)
{
	free(Table->Entries);
	free(Table);

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
					newSlot->KMer = entry->KMer;
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

ERR_VALUE kmer_table_insert(PKMER_TABLE Table, const PKMER KMer)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_table_get_slot(Table, KMer);
	if (entry != NULL) {
		if (_kmer_table_entry_empty(entry)) {
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
