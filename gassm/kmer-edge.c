
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"


/************************************************************************/
/*                         HELPER MACROS                                */
/************************************************************************/

#define _kmer_edge_table_entry_empty(aEntry)						((aEntry)->Source == NULL)
#define _next_hash_attempt(aHash, aAttempt, aModulus)				((aHash + 2 * aAttempt + 1) % aModulus)


/************************************************************************/
/*                          GLOBAL VARIABLES                            */
/************************************************************************/

static unsigned int _lastOrder = 0;

/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

static PKMER_EDGE _kmer_edge_table_get_slot_hint(const PKMER_EDGE_TABLE Table, size_t Hash, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE ret = NULL;

	ret = Table->Entries + Hash;
	if (!_kmer_edge_table_entry_empty(ret) && (Source == NULL || !kmer_equal(ret->Source, Source) || !kmer_equal(ret->Dest, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE first = ret;

		do {
			Hash = _next_hash_attempt(Hash, attempt, Table->Size);
			ret = Table->Entries + Hash;
			if (_kmer_edge_table_entry_empty(ret) || (Source != NULL && kmer_equal(ret->Source, Source) && kmer_equal(ret->Dest, Dest)))
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

static PKMER_EDGE _kmer_edge_table_get_slot(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	size_t hash = 0;
	PKMER_EDGE ret = NULL;

	hash = kmer_edge_hash(Table, Source, Dest);
	assert(hash < Table->Size);
	ret = _kmer_edge_table_get_slot_hint(Table, hash, Source, Dest);

	return ret;
}


/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t X, const size_t Size, PKMER_EDGE_TABLE *Table)
{
	PKMER_EDGE_TABLE tmpTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (utils_is_prime(Size)) {
		tmpTable = (PKMER_EDGE_TABLE)malloc(sizeof(KMER_EDGE_TABLE));
		if (tmpTable != NULL) {
			tmpTable->Size = Size;
			tmpTable->X = X;
			tmpTable->KMerSize = KMerSize;
			tmpTable->PowX = utils_pow_mod(tmpTable->X, KMerSize - 1, tmpTable->Size);
			ret = utils_mul_inverse(X, Size, &tmpTable->Inverse);
			if (ret == ERR_SUCCESS) {
				assert((X*tmpTable->Inverse % Size) == 1);
				tmpTable->Entries = (PKMER_EDGE)calloc(tmpTable->Size, sizeof(KMER_EDGE));
				if (tmpTable->Entries != NULL) {
					memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_EDGE));
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

void kmer_edge_table_destroy(PKMER_EDGE_TABLE Table)
{
	free(Table->Entries);
	free(Table);

	return;
}


ERR_VALUE kmer_edge_table_extend(PKMER_EDGE_TABLE Table)
{
	size_t newSize = 0;
	PKMER_EDGE_TABLE newTable = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	newSize = Table->Size * 2;
	newSize = utils_next_prime(newSize);
	ret = kmer_edge_table_create(Table->KMerSize, Table->X, newSize, &newTable);
	if (ret == ERR_SUCCESS) {
		PKMER_EDGE entry = Table->Entries;
		PKMER_EDGE newSlot = NULL;

		for (size_t i = 0; i < Table->Size; ++i) {
			if (!_kmer_edge_table_entry_empty(entry)) {
				newSlot = _kmer_edge_table_get_slot(newTable, entry->Source, entry->Dest);
				if (newSlot != NULL)
					memcpy(newSlot, entry, sizeof(KMER_EDGE));
				else ret = ERR_TABLE_FULL;
			}

			if (ret != ERR_SUCCESS)
				break;

			++entry;
		}

		if (ret == ERR_SUCCESS) {
			free(Table->Entries);
			memcpy(Table, newTable, sizeof(KMER_EDGE_TABLE));
		}

		if (ret != ERR_SUCCESS)
			kmer_edge_table_destroy(newTable);
	}

	return ret;
}

PKMER_EDGE kmer_edge_table_get(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE ret = NULL;

	ret = _kmer_edge_table_get_slot(Table, Source, Dest);
	if (ret != NULL && _kmer_edge_table_entry_empty(ret))
		ret = NULL;

	return ret;
}


ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_edge_table_get_slot(Table, Source, Dest);
	if (entry != NULL) {
		if (_kmer_edge_table_entry_empty(entry)) {
			entry->Source = kmer_copy(Source);
			if (entry->Source != NULL) {
				entry->Dest = kmer_copy(Dest);
				ret = (entry->Dest != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
				if (ret == ERR_SUCCESS) {
					entry->Order = _lastOrder;
					++_lastOrder;
				}
			} else ret = ERR_OUT_OF_MEMORY;
		} else ret = ERR_ALREADY_EXISTS;
	} else ret = ERR_TABLE_FULL;

	return ret;
}

size_t kmer_edge_hash(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
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

void kmer_edge_table_print(const PKMER_EDGE_TABLE Table)
{
	PKMER_EDGE edge = NULL;

	for (unsigned int i = 0; i < _lastOrder; ++i) {
		edge = Table->Entries;
		for (size_t j = 0; j < Table->Size; ++j) {
			if (!_kmer_edge_table_entry_empty(edge) && edge->Order == i) {
				printf("\t");
				kmer_print(edge->Source);
				printf(" -> ");
				kmer_print(edge->Dest);
				printf(" [weight=%u];\n", edge->Weight);
			}

			++edge;
		}
	}


	return;
}
