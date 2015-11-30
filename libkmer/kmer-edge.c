
#include <stdint.h>
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
/*                          GLOBAL VARIABLES                            */
/************************************************************************/

static unsigned int _lastOrder = 0;

/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE_TABLE)
static UTILS_TYPED_CALLOC_FUNCTION(KMER_EDGE)


static PKMER_EDGE _kmer_edge_table_get_slot_insert_hint(const PKMER_EDGE_TABLE Table, size_t Hash, const PKMER Source, const PKMER Dest)
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


static PKMER_EDGE _kmer_edge_table_get_slot_delsearch_hint(const PKMER_EDGE_TABLE Table, size_t Hash, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE ret = NULL;

	ret = Table->Entries + Hash;
	if ((!_kmer_edge_table_entry_empty(ret) || ret->Deleted) && (ret->Deleted || !kmer_equal(ret->Source, Source) || !kmer_equal(ret->Dest, Dest))) {
		size_t attempt = 1;
		PKMER_EDGE first = ret;

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


static PKMER_EDGE _kmer_edge_table_get_slot(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest, const ETableOpType OpType)
{
	size_t hash = 0;
	PKMER_EDGE ret = NULL;

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


/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t X, const size_t Size, PKMER_EDGE_TABLE *Table)
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
			ret = utils_mul_inverse(X, Size, &tmpTable->Inverse);
			if (ret == ERR_SUCCESS) {
				assert((X*tmpTable->Inverse % Size) == 1);
				ret = utils_calloc_KMER_EDGE(tmpTable->Size, &tmpTable->Entries);
				if (ret == ERR_SUCCESS) {
					memset(tmpTable->Entries, 0, tmpTable->Size*sizeof(KMER_EDGE));
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
	PKMER_EDGE entry = Table->Entries;

	for (size_t i = 0; i < Table->Size; ++i) {
		if (!_kmer_edge_table_entry_empty(entry)) {
			kmer_free(entry->Source);
			kmer_free(entry->Dest);
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
	ret = kmer_edge_table_create(Table->KMerSize, Table->X, newSize, &newTable);
	if (ret == ERR_SUCCESS) {
		PKMER_EDGE entry = Table->Entries;
		PKMER_EDGE newSlot = NULL;

		for (size_t i = 0; i < Table->Size; ++i) {
			if (!_kmer_edge_table_entry_empty(entry)) {
				newSlot = _kmer_edge_table_get_slot(newTable, entry->Source, entry->Dest, totInsert);
				if (newSlot != NULL)
					memcpy(newSlot, entry, sizeof(KMER_EDGE));
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


ERR_VALUE kmer_edge_table_copy(const PKMER_EDGE_TABLE Source, PKMER_EDGE_TABLE * Copied)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE tmpTable = NULL;

	ret = utils_malloc_KMER_EDGE_TABLE(&tmpTable);
	if (ret == ERR_SUCCESS) {
		tmpTable->Inverse = Source->Inverse;
		tmpTable->KMerSize = Source->KMerSize;
		tmpTable->PowX = Source->PowX;
		tmpTable->Size = Source->Size;
		tmpTable->X = Source->X;
		ret = utils_calloc_KMER_EDGE(tmpTable->Size, &tmpTable->Entries);
		if (ret == ERR_SUCCESS) {
			PKMER_EDGE sourceEntry = Source->Entries;
			PKMER_EDGE destEntry = tmpTable->Entries;

			memset(tmpTable->Entries, 0, sizeof(KMER_EDGE)*tmpTable->Size);
			for (size_t i = 0; i < tmpTable->Size; ++i) {
				if (!_kmer_edge_table_entry_empty(sourceEntry) || sourceEntry->Deleted) {
					memcpy(destEntry, sourceEntry, sizeof(KMER_EDGE));
					if (!sourceEntry->Deleted) {
						destEntry->Source = kmer_copy(sourceEntry->Source);
						if (destEntry->Source != NULL) {
							destEntry->Dest = kmer_copy(sourceEntry->Dest);
							if (destEntry->Dest == NULL) {
								kmer_free(destEntry->Source);
								ret = ERR_OUT_OF_MEMORY;
							}
						} else ret = ERR_OUT_OF_MEMORY;
					}
				}

				if (ret != ERR_SUCCESS) {
					--destEntry;
					for (size_t j = 0; j < i; ++j) {
						if (!_kmer_edge_table_entry_empty(destEntry)) {
							kmer_free(destEntry->Dest);
							kmer_free(destEntry->Source);
						}

						--destEntry;
					}
					
					break;
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


PKMER_EDGE kmer_edge_table_get(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE ret = NULL;

	ret = _kmer_edge_table_get_slot(Table, Source, Dest, totSearch);
	if (ret != NULL && _kmer_edge_table_entry_empty(ret))
		ret = NULL;

	return ret;
}


void kmer_edge_table_delete_by_entry(PKMER_EDGE_TABLE Table, PKMER_EDGE Entry)
{
	kmer_free(Entry->Source);
	kmer_free(Entry->Dest);
	memset(Entry, 0, sizeof(KMER_EDGE));
	Entry->Deleted = TRUE;

	return;
}


ERR_VALUE kmer_edge_table_delete(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	edge = _kmer_edge_table_get_slot(Table, Source, Dest, totDelete);
	if (edge != NULL && !_kmer_edge_table_entry_empty(edge))
		kmer_edge_table_delete_by_entry(Table, edge);
	else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_edge_table_insert_ex(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest, PKMER_EDGE *Edge)
{
	PKMER_EDGE entry = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	entry = _kmer_edge_table_get_slot(Table, Source, Dest, totInsert);
	if (entry != NULL) {
		if (_kmer_edge_table_entry_empty(entry)) {
			entry->Deleted = FALSE;
			entry->Source = kmer_copy(Source);
			if (entry->Source != NULL) {
				entry->Dest = kmer_copy(Dest);
				entry->Type = kmetReference;
				entry->Probability = 0;
				ret = (entry->Dest != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
				if (ret == ERR_SUCCESS) {
					entry->Order = _lastOrder;
					*Edge = entry;
					++_lastOrder;
				}
			} else ret = ERR_OUT_OF_MEMORY;
		} else {
			*Edge = entry;
			ret = ERR_ALREADY_EXISTS;
		}
	} else ret = ERR_TABLE_FULL;

	return ret;
}


ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest)
{
	PKMER_EDGE dummy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_edge_table_insert_ex(Table, Source, Dest, &dummy);

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

void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table)
{
	PKMER_EDGE edge = NULL;

	edge = Table->Entries;
	for (size_t j = 0; j < Table->Size; ++j) {
		if (!_kmer_edge_table_entry_empty(edge)) {
			fprintf(Stream, "\t");
			kmer_print(Stream, edge->Source);
			fprintf(Stream, " -> ");
			kmer_print(Stream, edge->Dest);
			fprintf(Stream, " [");
			fprintf(Stream, "weight=%lu", edge->Weight);
			fprintf(Stream, ",label=\"L: %u; W: %li; P: %u\"", edge->Length, edge->Weight, (uint8_t)(edge->Probability * 100));
			fprintf(Stream, "];\n");
		}

		++edge;
	}

	return;
}


ERR_VALUE kmer_edge_table_first(const PKMER_EDGE_TABLE Table, PKMER_EDGE *Slot)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE entry = Table->Entries;

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


ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const PKMER_EDGE Current, PKMER_EDGE *Next)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE entry = Current + 1;

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
