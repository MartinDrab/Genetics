
#ifndef __DYM_ARRAY_H__
#define __DYM_ARRAY_H__




typedef struct _DYM_ARRAY {
	void **Data;
	size_t ValidLength;
	size_t AllocatedLength;
	uint8_t PreallocRatio;
} DYM_ARRAY, *PDYM_ARRAY;

#define dym_array_size(aArray)				((aArray)->ValidLength)
#define dym_array_alloc_size(aArray)		((aArray)->AllocatedLength)
#define dym_array_data(aArray)				((aArray)->Data)


void dym_array_create(PDYM_ARRAY Array, const uint8_t PreallocRatio);
void dym_array_destroy(PDYM_ARRAY Array);
ERR_VALUE dym_array_reserve(PDYM_ARRAY Array, const size_t Items);
ERR_VALUE dym_array_push_back(PDYM_ARRAY Array, const void *Value);
void dym_array_push_back_no_alloc(PDYM_ARRAY Array, const void *Value);
void *dym_array_get(const struct _DYM_ARRAY *Array, const size_t Index);
void dym_array_set(PDYM_ARRAY Array, const size_t Index, const void *Value);
void *dym_array_remove(PDYM_ARRAY Array, const size_t Index);
void *dym_array_remove_back(PDYM_ARRAY Array);
ERR_VALUE dym_array_to_array(const PDYM_ARRAY Array, void ***NormalArray);
ERR_VALUE dym_array_compact(PDYM_ARRAY Array);
ERR_VALUE dym_array_copy(PDYM_ARRAY Dest, const DYM_ARRAY *Source);
ERR_VALUE dym_array_prepare_for_insert(PDYM_ARRAY Array, const size_t NumberOfItems);
void dym_array_replace(PDYM_ARRAY Array, const void *Item, const void *New);
void dym_array_remove_by_item_fast(PDYM_ARRAY Array, const void *Item);



#endif
