
#include <assert.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "dym-array.h"


/************************************************************************/
/*                       HELPER FUNCTIONS                               */
/************************************************************************/


/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

void dym_array_create(PDYM_ARRAY Array, const uint8_t PreallocRatio)
{
	memset(Array, 0, sizeof(DYM_ARRAY));
	Array->PreallocRatio = (PreallocRatio == 0 || PreallocRatio < 100) ? 140 : PreallocRatio;

	return;
}


void dym_array_destroy(PDYM_ARRAY Array)
{
	if (Array->AllocatedLength > 0)
		utils_free(Array->Data);

	memset(Array, 0, sizeof(DYM_ARRAY));

	return;
}


ERR_VALUE dym_array_reserve(PDYM_ARRAY Array, const size_t Items)
{
	void *tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(Items, sizeof(void *), &tmp);
	if (ret == ERR_SUCCESS) {
		Array->AllocatedLength = Items;
		if (Array->ValidLength > 0) {
			Array->ValidLength = min(Array->ValidLength, Items);;
			memcpy(tmp, Array->Data, Array->ValidLength*sizeof(void *));
		}

		utils_free(Array->Data);
		Array->Data = (void **)tmp;
	}

	return ret;
}


ERR_VALUE dym_array_push_back(PDYM_ARRAY Array, const void *Value)
{
	size_t newSize = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (Array->ValidLength == Array->AllocatedLength) {
		// WARNING: Be aware of a possible overflow
		newSize = (Array->AllocatedLength * Array->PreallocRatio / 100) + 5;
		ret = dym_array_reserve(Array, newSize);
	}
	
	if (ret == ERR_SUCCESS)
		dym_array_push_back_no_alloc(Array, Value);

	return ret;
}


void dym_array_push_back_no_alloc(PDYM_ARRAY Array, const void *Value)
{
	assert(Array->AllocatedLength > Array->ValidLength);
	Array->Data[Array->ValidLength] = (void *)Value;
	++Array->ValidLength;

	return;
}


void *dym_array_get(const PDYM_ARRAY Array, const size_t Index)
{
	assert(Index < Array->ValidLength);
	return Array->Data[Index];
}


void dym_array_set(PDYM_ARRAY Array, const size_t Index, const void *Value)
{
	assert(Index < Array->ValidLength);
	Array->Data[Index] = (void *)Value;

	return;
}

void *dym_array_remove(PDYM_ARRAY Array, const size_t Index)
{
	void *ret = NULL;

	assert(Index < Array->ValidLength);
	ret = Array->Data[Index];
	memmove(Array->Data + Index, Array->Data + Index + 1, (Array->ValidLength - Index - 1)*sizeof(void *));

	return ret;
}



void *dym_array_remove_back(PDYM_ARRAY Array)
{
	assert(Array->ValidLength > 0);
	--Array->ValidLength;

	return Array->Data[Array->ValidLength];
}


ERR_VALUE dym_array_to_array(const PDYM_ARRAY Array, void ***NormalArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	void ** tmpArray = NULL;

	ret = utils_calloc(Array->ValidLength, sizeof(void *), (void **)&tmpArray);
	if (ret == ERR_SUCCESS) {
		memcpy(tmpArray, Array->Data, Array->ValidLength*sizeof(void *));
		*NormalArray = tmpArray;
	}

	return ret;
}


ERR_VALUE dym_array_compact(PDYM_ARRAY Array)
{
	return dym_array_reserve(Array, Array->ValidLength);
}


ERR_VALUE dym_array_copy(PDYM_ARRAY Dest, const DYM_ARRAY *Source)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (Dest->AllocatedLength < Source->AllocatedLength) {
		void **tmp = NULL;

		ret = utils_calloc(Source->ValidLength, sizeof(void *), (void **)&tmp);
		if (ret == ERR_SUCCESS) {
			utils_free(Dest->Data);
			Dest->Data = tmp;
			Dest->AllocatedLength = Source->ValidLength;
		}
	}

	if (ret == ERR_SUCCESS) {
		memcpy(Dest->Data, Source->Data, Source->ValidLength*sizeof(void *));
		Dest->ValidLength = Source->ValidLength;
	}

	return ret;
}


ERR_VALUE dym_array_prepare_for_insert(PDYM_ARRAY Array, const size_t NumberOfItems)
{
	ERR_VALUE ret = ERR_SUCCESS;

	ret = ERR_SUCCESS;
	if (Array->AllocatedLength < Array->ValidLength + NumberOfItems)
		ret = dym_array_reserve(Array, Array->ValidLength + NumberOfItems);

	return ret;
}


void dym_array_replace(PDYM_ARRAY Array, const void *Item, const void *New)
{
	boolean found = FALSE;

	for (size_t i = 0; i < Array->ValidLength; ++i) {
		found = Array->Data[i] == Item;
		if (found) {
			Array->Data[i] = (void *)New;
			break;
		}
	}

	if (!found) {
		assert(FALSE);
	}

	return;
}


void dym_array_remove_by_item_fast(PDYM_ARRAY Array, const void *Item)
{
	boolean found = FALSE;

	for (size_t i = 0; i < Array->ValidLength; ++i) {
		found = Array->Data[i] == Item;
		if (found) {
			--Array->ValidLength;
			Array->Data[i] = Array->Data[Array->ValidLength];
			break;
		}
	}

	if (!found) {
		assert(FALSE);
	}

	return;
}

