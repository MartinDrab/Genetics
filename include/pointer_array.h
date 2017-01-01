
#ifndef __POINTER_ARRAY_H__
#define __POINTER_ARRAY_H__


#include <assert.h>
#include "err.h"
#include "utils.h"

#define POINTER_ARRAY_STATIC_ALLOC							4

#define POINTER_ARRAY_TYPE(aDataType)						POINTER_ARRAY_##aDataType
#define POINTER_ARRAY_PTYPE(aDataType)						PPOINTER_ARRAY_##aDataType

#define POINTER_ARRAY_TYPEDEF(aDataType)					\
	typedef struct _POINTER_ARRAY_##aDataType {				\
		size_t ValidLength;								\
		size_t AllocLength;								\
		aDataType **Data;								\
		size_t Ratio;									\
		aDataType *Storage[POINTER_ARRAY_STATIC_ALLOC];							\
	} POINTER_ARRAY_TYPE(aDataType), *PPOINTER_ARRAY_##aDataType	\


#define pointer_array_size(aArray)							((aArray)->ValidLength)

#define POINTER_ARRAY_IMPLEMENTATION(aDataType)									\
	POINTER_ARRAY_ALLOC_FUNCTION(aDataType)										\
	POINTER_ARRAY_INIT_FUNCTION(aDataType)										\
	POINTER_ARRAY_FINIT_FUNCTION(aDataType)										\
	POINTER_ARRAY_RESERVE_FUNCTION(aDataType)									\
	POINTER_ARRAY_PUSH_BACK_NO_ALLOC_FUNCTION(aDataType)						\
	POINTER_ARRAY_PUSH_BACK_FUNCTION(aDataType)									\
	POINTER_ARRAY_ITEM_FUNCTION(aDataType)										\
	POINTER_ARRAY_CONST_ITEM_FUNCTION(aDataType)								\
	POINTER_ARRAY_POP_BACK_FUNCTION(aDataType)									\
	POINTER_ARRAY_CLEAR_FUNCTION(aDataType)										\
	POINTER_ARRAY_REMOVE_FAST_FUNCTION(aDataType)								\
	POINTER_ARRAY_EXCHANGE_FUNCTION(aDataType)									\
	POINTER_ARRAY_PUSH_BACK_ARRAY_FUNCTION(aDataType)							\
	POINTER_ARRAY_CONTAINS_FUNCTION(aDataType)									\
	POINTER_ARRAY_CLEAN_COPY_FUNCTION(aDataType)								\
	POINTER_ARRAY_REMOVE_BY_ITEM_FAST(aDataType)								\
	POINTER_ARRAY_REMOVE_FUNCTION(aDataType)									\


#define POINTER_ARRAY_ALLOC_FUNCTION(aDataType)						\
	INLINE_FUNCTION ERR_VALUE pointer_array_alloc_data_##aDataType(const size_t Count, aDataType ***Data)	\
	{																			\
		return utils_calloc(Count, sizeof(aDataType *), (void **)Data);			\
	}																			\

#define POINTER_ARRAY_INIT_FUNCTION(aDataType)												\
	INLINE_FUNCTION void pointer_array_init_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const size_t Ratio)	\
	{																					\
		memset(Array, 0, sizeof(POINTER_ARRAY_TYPE(aDataType)));								\
		Array->Ratio = Ratio;															\
		Array->Data = Array->Storage;													\
		Array->AllocLength = sizeof(Array->Storage) / sizeof(Array->Storage[0]);				\
																						\
		return;																			\
	}																					\

#define POINTER_ARRAY_FINIT_FUNCTION(aDataType)												\
	INLINE_FUNCTION void pointer_array_finit_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array)			\
	{																					\
		if (Array->AllocLength > 0 && Array->Data != Array->Storage)					\
			utils_free(Array->Data);													\
																						\
		Array->Data = NULL;																\
		return;																			\
	}																					\

#define POINTER_ARRAY_RESERVE_FUNCTION(aDataType)											\
	INLINE_FUNCTION ERR_VALUE pointer_array_reserve_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const size_t Count)	\
	{																							\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;														\
		aDataType **newData = NULL;																\
																								\
		if (Count <= Array->AllocLength)														\
			return ERR_SUCCESS;																	\
																								\
		ret = pointer_array_alloc_data_##aDataType(Count, &newData);						\
		if (ret == ERR_SUCCESS) {																\
			memcpy(newData, Array->Data, Array->ValidLength*sizeof(aDataType*));					\
			if (Array->AllocLength > 0 && Array->Data != Array->Storage)						\
				utils_free(Array->Data);														\
																								\
			Array->AllocLength = Count;															\
			Array->Data = newData;																\
		}																						\
																								\
		return ret;																				\
	}																							\

#define POINTER_ARRAY_PUSH_BACK_NO_ALLOC_FUNCTION(aDataType)										\
	INLINE_FUNCTION void pointer_array_push_back_no_alloc_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const aDataType *Value)	\
	{																						\
		if (Array->ValidLength >= Array->AllocLength) exit(256);							\
		memcpy(Array->Data + Array->ValidLength, &Value, sizeof(Value));					\
		++Array->ValidLength;																\
																							\
		return;																				\
	}																						\

#define POINTER_ARRAY_PUSH_BACK_FUNCTION(aDataType)												\
	INLINE_FUNCTION ERR_VALUE pointer_array_push_back_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const aDataType *Value)	\
	{																						\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;													\
																							\
		if (Array->ValidLength < Array->AllocLength) {										\
			pointer_array_push_back_no_alloc_##aDataType(Array, Value);							\
			ret = ERR_SUCCESS;																\
		} else {																			\
			ret = pointer_array_reserve_##aDataType(Array, (uint64_t)Array->AllocLength*Array->Ratio/100);	\
			if (ret == ERR_SUCCESS)															\
				pointer_array_push_back_no_alloc_##aDataType(Array, Value);						\
		}																					\
																							\
		return ret;																			\
	}

#define POINTER_ARRAY_ITEM_FUNCTION(aDataType)													\
	INLINE_FUNCTION aDataType **pointer_array_item_##aDataType(const POINTER_ARRAY_TYPE(aDataType) *Array, const size_t Index)				\
	{																						\
		assert(Array->ValidLength > Index);													\
																							\
		return (Array->Data + Index);														\
	}																						\

#define POINTER_ARRAY_CONST_ITEM_FUNCTION(aDataType)											\
	INLINE_FUNCTION const aDataType **pointer_array_const_item_##aDataType(const POINTER_ARRAY_TYPE(aDataType) *Array, const size_t Index)				\
	{																						\
		assert(Array->ValidLength > Index);													\
																							\
		return (const aDataType **)(Array->Data + Index);														\
	}																						\

#define POINTER_ARRAY_POP_BACK_FUNCTION(aDataType)												\
	INLINE_FUNCTION aDataType **pointer_array_pop_back_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array)	\
	{																						\
		aDataType **ret = pointer_array_item_##aDataType(Array, Array->ValidLength - 1);			\
																							\
		assert(Array->ValidLength > 0);														\
		--Array->ValidLength;																\
																							\
		return ret;																			\
	}																						\

#define POINTER_ARRAY_CLEAR_FUNCTION(aDataType)													\
	INLINE_FUNCTION void pointer_array_clear_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array)				\
	{																						\
		Array->ValidLength = 0;																\
																							\
		return;																				\
	}																						\

#define POINTER_ARRAY_REMOVE_FAST_FUNCTION(aDataType)											\
	INLINE_FUNCTION void pointer_array_remove_fast_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const size_t Index) \
	{																							\
		assert(Array->ValidLength > Index);														\
		memmove(Array->Data + Index, Array->Data + Array->ValidLength - 1, sizeof(aDataType*));	\
		--Array->ValidLength;																	\
																								\
		return;																					\
	}																							\

#define POINTER_ARRAY_EXCHANGE_FUNCTION(aDataType)					\
	INLINE_FUNCTION void pointer_array_exchange_##aDataType(POINTER_ARRAY_PTYPE(aDataType) A1, POINTER_ARRAY_PTYPE(aDataType) A2)	\
	{											\
		POINTER_ARRAY_TYPE(aDataType) tmp;			\
												\
		memcpy(&tmp, A1, sizeof(tmp));			\
		memcpy(A1, A2, sizeof(tmp));			\
		memcpy(A2, &tmp, sizeof(tmp));			\
		if (A1->Data == A2->Storage)			\
			A1->Data = A1->Storage;				\
												\
		if (A2->Data == A1->Storage)			\
			A2->Data = A2->Storage;				\
												\
		return;									\
	}											\

#define POINTER_ARRAY_PUSH_BACK_ARRAY_FUNCTION(aDataType)	\
	INLINE_FUNCTION ERR_VALUE pointer_array_push_back_array_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Target, const POINTER_ARRAY_TYPE(aDataType) *Source)	\
	{																								\
		ERR_VALUE ret = ERR_SUCCESS;																\
																									\
		if (Target->AllocLength < Target->ValidLength + Source->ValidLength)						\
			ret = pointer_array_reserve_##aDataType(Target, 1 + Target->ValidLength + Source->ValidLength);	\
																									\
		if (ret == ERR_SUCCESS) {																	\
			memcpy(Target->Data + Target->ValidLength, Source->Data, Source->ValidLength*sizeof(aDataType*));	\
			Target->ValidLength += Source->ValidLength;	\
		}												\
														\
		return ret;										\
	}													\

#define POINTER_ARRAY_CONTAINS_FUNCTION(aDataType)	\
	INLINE_FUNCTION boolean pointer_array_contains_##aDataType(const POINTER_ARRAY_TYPE(aDataType) *Array, const aDataType *Item)	\
	{																	\
		boolean ret = FALSE;											\
																		\
		for (size_t i = 0; i < Array->ValidLength; ++i) {				\
			ret = (memcmp(Array->Data + i, &Item, sizeof(Item)) == 0);	\
			if (ret)													\
				break;												\
		}																\
																		\
		return ret;														\
	}																	\

#define POINTER_ARRAY_CLEAN_COPY_FUNCTION(aDataType)	\
	INLINE_FUNCTION ERR_VALUE pointer_array_clean_copy_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Dest, const POINTER_ARRAY_TYPE(aDataType) *Source) \
	{																	\
		pointer_array_clear_##aDataType(Dest);							\
		return pointer_array_push_back_array_##aDataType(Dest, Source);	\
	}																	\

#define POINTER_ARRAY_REMOVE_BY_ITEM_FAST(aDataType)	\
	INLINE_FUNCTION boolean pointer_array_remove_by_item_fast_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const aDataType *Item)	\
	{																	\
		for (size_t i = 0; i < Array->ValidLength; ++i) {				\
			if (memcmp(Array->Data + i, &Item, sizeof(Item)) == 0) {	\
				pointer_array_remove_fast_##aDataType(Array, i);		\
				return TRUE;													\
			}															\
		}																\
																		\
		return FALSE;															\
	}																	\

#define POINTER_ARRAY_REMOVE_FUNCTION(aDataType)	\
	INLINE_FUNCTION void pointer_array_remove_##aDataType(POINTER_ARRAY_PTYPE(aDataType) Array, const aDataType *Item)	\
	{																	\
		for (size_t i = 0; i < Array->ValidLength; ++i) {				\
			if (memcmp(Array->Data + i, &Item, sizeof(Item)) == 0) {	\
				memmove(Array->Data + i, Array->Data + i + 1, sizeof(Item)*(Array->ValidLength - (i + 1)));	\
				--Array->ValidLength;									\
			}															\
		}																\
																		\
		return;														\
	}																	\


POINTER_ARRAY_TYPEDEF(int8_t);
POINTER_ARRAY_IMPLEMENTATION(int8_t)
POINTER_ARRAY_TYPEDEF(uint8_t);
POINTER_ARRAY_IMPLEMENTATION(uint8_t)
POINTER_ARRAY_TYPEDEF(int16_t);
POINTER_ARRAY_IMPLEMENTATION(int16_t)
POINTER_ARRAY_TYPEDEF(uint16_t);
POINTER_ARRAY_IMPLEMENTATION(uint16_t)
POINTER_ARRAY_TYPEDEF(int32_t);
POINTER_ARRAY_IMPLEMENTATION(int32_t)
POINTER_ARRAY_TYPEDEF(uint32_t);
POINTER_ARRAY_IMPLEMENTATION(uint32_t)
POINTER_ARRAY_TYPEDEF(int64_t);
POINTER_ARRAY_IMPLEMENTATION(int64_t)
POINTER_ARRAY_TYPEDEF(uint64_t);
POINTER_ARRAY_IMPLEMENTATION(uint64_t)
POINTER_ARRAY_TYPEDEF(size_t);
POINTER_ARRAY_IMPLEMENTATION(size_t)
POINTER_ARRAY_TYPEDEF(char);
POINTER_ARRAY_IMPLEMENTATION(char)
POINTER_ARRAY_TYPEDEF(float);
POINTER_ARRAY_IMPLEMENTATION(float)
POINTER_ARRAY_TYPEDEF(double);
POINTER_ARRAY_IMPLEMENTATION(double)
POINTER_ARRAY_TYPEDEF(boolean);
POINTER_ARRAY_IMPLEMENTATION(boolean)



#endif
