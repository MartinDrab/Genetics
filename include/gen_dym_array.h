
#ifndef __GEN_DYM_ARRAY_H__
#define __GEN_DYM_ARRAY_H__


#include <assert.h>
#include "err.h"
#include "utils.h"


#define GEN_ARRAY_STATIC_ALLOC							10

#define GEN_ARRAY_TYPE(aDataType)						GEN_ARRAY_##aDataType
#define GEN_ARRAY_PTYPE(aDataType)						PGEN_ARRAY_##aDataType

#define GEN_ARRAY_TYPEDEF(aDataType)					\
	typedef struct _GEN_ARRAY_##aDataType {				\
		size_t ValidLength;								\
		size_t AllocLength;								\
		aDataType *Data;								\
		size_t Ratio;									\
		aDataType Storage[GEN_ARRAY_STATIC_ALLOC];							\
	} GEN_ARRAY_TYPE(aDataType), *PGEN_ARRAY_##aDataType	\


#define gen_array_size(aArray)							((aArray)->ValidLength)

#define GEN_ARRAY_IMPLEMENTATION(aDataType)									\
	GEN_ARRAY_ALLOC_FUNCTION(aDataType)										\
	GEN_ARRAY_INIT_FUNCTION(aDataType)										\
	GEN_ARRAY_FINIT_FUNCTION(aDataType)										\
	GEN_ARRAY_RESERVE_FUNCTION(aDataType)									\
	GEN_ARRAY_PUSH_BACK_NO_ALLOC_FUNCTION(aDataType)						\
	GEN_ARRAY_PUSH_BACK_FUNCTION(aDataType)									\
	GEN_ARRAY_ITEM_FUNCTION(aDataType)										\
	GEN_ARRAY_CONST_ITEM_FUNCTION(aDataType)								\
	GEN_ARRAY_POP_BACK_FUNCTION(aDataType)									\
	GEN_ARRAY_CLEAR_FUNCTION(aDataType)										\
	GEN_ARRAY_REMOVE_FAST_FUNCTION(aDataType)								\
	GEN_ARRAY_EXCHANGE_FUNCTION(aDataType)									\
	GEN_ARRAY_PUSH_BACK_ARRAY_FUNCTION(aDataType)							\
	GEN_ARRAY_CONTAINS_FUNCTION(aDataType)									\
	GEN_ARRAY_REMOVE_FUNCTION(aDataType)									\
	GEN_ARRAY_PUSH_BACK_ARRAY_NO_ALLOC_FUNCTION(aDataType)					\


#define GEN_ARRAY_ALLOC_FUNCTION(aDataType)						\
	INLINE_FUNCTION ERR_VALUE dym_array_alloc_data_##aDataType(const size_t Count, aDataType **Data)	\
	{																			\
		return utils_calloc(Count, sizeof(aDataType), (void **)Data);			\
	}																			\

#define GEN_ARRAY_INIT_FUNCTION(aDataType)												\
	INLINE_FUNCTION void dym_array_init_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, const size_t Ratio)	\
	{																					\
		memset(Array, 0, sizeof(GEN_ARRAY_TYPE(aDataType)));								\
		Array->Ratio = Ratio;															\
		Array->Data = Array->Storage;													\
		Array->AllocLength = sizeof(Array->Storage) / sizeof(aDataType);				\
																						\
		return;																			\
	}																					\

#define GEN_ARRAY_FINIT_FUNCTION(aDataType)												\
	INLINE_FUNCTION void dym_array_finit_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array)	\
	{																					\
		if (Array->AllocLength > 0 && Array->Data != Array->Storage)					\
			utils_free(Array->Data);													\
																						\
		Array->Data = NULL;																\
		return;																			\
	}																					\

#define GEN_ARRAY_RESERVE_FUNCTION(aDataType)											\
	INLINE_FUNCTION ERR_VALUE dym_array_reserve_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, const size_t Count)	\
	{																							\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;														\
		aDataType *newData = NULL;																\
																								\
		if (Count <= Array->AllocLength)														\
			return ERR_SUCCESS;																	\
																								\
		ret = dym_array_alloc_data_##aDataType(Count, &newData);						\
		if (ret == ERR_SUCCESS) {																\
			memcpy(newData, Array->Data, Array->ValidLength*sizeof(aDataType));					\
			if (Array->AllocLength > 0 && Array->Data != Array->Storage)						\
				utils_free(Array->Data);														\
																								\
			Array->AllocLength = Count;															\
			Array->Data = newData;																\
		}																						\
																								\
		return ret;																				\
	}																							\

#define GEN_ARRAY_PUSH_BACK_NO_ALLOC_FUNCTION(aDataType)										\
	INLINE_FUNCTION void dym_array_push_back_no_alloc_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, const aDataType Value)	\
	{																						\
		assert(Array->ValidLength < Array->AllocLength);						\
		memcpy(Array->Data + Array->ValidLength, &Value, sizeof(aDataType));				\
		++Array->ValidLength;																\
																							\
		return;																				\
	}																						\

#define GEN_ARRAY_PUSH_BACK_FUNCTION(aDataType)												\
	INLINE_FUNCTION ERR_VALUE dym_array_push_back_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, const aDataType Value)	\
	{																						\
		ERR_VALUE ret = ERR_INTERNAL_ERROR;													\
																							\
		if (Array->ValidLength < Array->AllocLength) {										\
			dym_array_push_back_no_alloc_##aDataType(Array, Value);							\
			ret = ERR_SUCCESS;																\
		} else {																			\
			ret = dym_array_reserve_##aDataType(Array, Array->AllocLength*Array->Ratio/100);	\
			if (ret == ERR_SUCCESS)															\
				dym_array_push_back_no_alloc_##aDataType(Array, Value);						\
		}																					\
																							\
		return ret;																			\
	}

#define GEN_ARRAY_ITEM_FUNCTION(aDataType)													\
	INLINE_FUNCTION aDataType *dym_array_item_##aDataType(const GEN_ARRAY_TYPE(aDataType) *Array, const size_t Index)				\
	{																						\
		assert(Array->ValidLength > Index);													\
																							\
		return (Array->Data + Index);														\
	}																						\

#define GEN_ARRAY_CONST_ITEM_FUNCTION(aDataType)											\
	INLINE_FUNCTION const aDataType *dym_array_const_item_##aDataType(const GEN_ARRAY_TYPE(aDataType) *Array, const size_t Index)				\
	{																						\
		assert(Array->ValidLength > Index);													\
																							\
		return (const aDataType *)(Array->Data + Index);														\
	}																						\

#define GEN_ARRAY_POP_BACK_FUNCTION(aDataType)												\
	INLINE_FUNCTION aDataType *dym_array_pop_back_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array)	\
	{																						\
		aDataType *ret = dym_array_item_##aDataType(Array, Array->ValidLength - 1);			\
																							\
		assert(Array->ValidLength > 0);														\
		--Array->ValidLength;																\
																							\
		return ret;																			\
	}																						\

#define GEN_ARRAY_CLEAR_FUNCTION(aDataType)													\
	INLINE_FUNCTION void dym_array_clear_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array)		\
	{																						\
		Array->ValidLength = 0;																\
																							\
		return;																				\
	}																						\

#define GEN_ARRAY_REMOVE_FAST_FUNCTION(aDataType)											\
	INLINE_FUNCTION void dym_array_remove_fast##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, const size_t Index) \
	{																							\
		assert(Array->ValidLength > Index);														\
		memmove(Array->Data + Index, Array->Data + Array->ValidLength - 1, sizeof(aDataType));	\
		--Array->ValidLength;																	\
																								\
		return;																					\
	}																							\

#define GEN_ARRAY_EXCHANGE_FUNCTION(aDataType)					\
	INLINE_FUNCTION void dym_array_exchange_##aDataType(GEN_ARRAY_PTYPE(aDataType) A1, GEN_ARRAY_PTYPE(aDataType) A2)	\
	{											\
		GEN_ARRAY_TYPE(aDataType) tmp;			\
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

#define GEN_ARRAY_PUSH_BACK_ARRAY_FUNCTION(aDataType)	\
	INLINE_FUNCTION ERR_VALUE dym_array_push_back_array_##aDataType(GEN_ARRAY_PTYPE(aDataType) Target, const GEN_ARRAY_TYPE(aDataType) *Source)	\
	{																								\
		ERR_VALUE ret = ERR_SUCCESS;																\
																									\
		if (Target->AllocLength < Target->ValidLength + Source->ValidLength)						\
			ret = dym_array_reserve_##aDataType(Target, 1 + Target->ValidLength + Source->ValidLength);	\
																									\
		if (ret == ERR_SUCCESS) {																	\
			memcpy(Target->Data + Target->ValidLength, Source->Data, Source->ValidLength*sizeof(aDataType));	\
			Target->ValidLength += Source->ValidLength;	\
		}												\
														\
		return ret;										\
	}													\

#define GEN_ARRAY_PUSH_BACK_ARRAY_NO_ALLOC_FUNCTION(aDataType)	\
	INLINE_FUNCTION void dym_array_push_back_array_no_alloc_##aDataType(GEN_ARRAY_PTYPE(aDataType) Target, const GEN_ARRAY_TYPE(aDataType) *Source)	\
	{																								\
		assert(Target->AllocLength >= Target->ValidLength + Source->ValidLength);						\
		memcpy(Target->Data + Target->ValidLength, Source->Data, Source->ValidLength*sizeof(aDataType));	\
		Target->ValidLength += Source->ValidLength;	\
														\
		return;										\
	}													\

#define GEN_ARRAY_CONTAINS_FUNCTION(aDataType)	\
	INLINE_FUNCTION boolean dym_array_contains_##aDataType(const GEN_ARRAY_TYPE(aDataType) *Array, aDataType Item)	\
	{																	\
		boolean ret = FALSE;											\
																		\
		for (size_t i = 0; i < Array->ValidLength; ++i) {				\
			ret = (memcmp(Array->Data + i, &Item, sizeof(Item)) == 0);	\
			if (ret)													\
				break;													\
		}																\
																		\
		return ret;														\
	}																	\

#define GEN_ARRAY_REMOVE_FUNCTION(aDataType)	\
	INLINE_FUNCTION void dym_array_remove_##aDataType(GEN_ARRAY_PTYPE(aDataType) Array, aDataType Item)	\
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

GEN_ARRAY_TYPEDEF(int8_t);
GEN_ARRAY_IMPLEMENTATION(int8_t)
GEN_ARRAY_TYPEDEF(uint8_t);
GEN_ARRAY_IMPLEMENTATION(uint8_t)
GEN_ARRAY_TYPEDEF(int16_t);
GEN_ARRAY_IMPLEMENTATION(int16_t)
GEN_ARRAY_TYPEDEF(uint16_t);
GEN_ARRAY_IMPLEMENTATION(uint16_t)
GEN_ARRAY_TYPEDEF(int32_t);
GEN_ARRAY_IMPLEMENTATION(int32_t)
GEN_ARRAY_TYPEDEF(uint32_t);
GEN_ARRAY_IMPLEMENTATION(uint32_t)
GEN_ARRAY_TYPEDEF(int64_t);
GEN_ARRAY_IMPLEMENTATION(int64_t)
GEN_ARRAY_TYPEDEF(uint64_t);
GEN_ARRAY_IMPLEMENTATION(uint64_t)
GEN_ARRAY_TYPEDEF(size_t);
GEN_ARRAY_IMPLEMENTATION(size_t)
GEN_ARRAY_TYPEDEF(char);
GEN_ARRAY_IMPLEMENTATION(char)
GEN_ARRAY_TYPEDEF(float);
GEN_ARRAY_IMPLEMENTATION(float)
GEN_ARRAY_TYPEDEF(double);
GEN_ARRAY_IMPLEMENTATION(double)
GEN_ARRAY_TYPEDEF(boolean);
GEN_ARRAY_IMPLEMENTATION(boolean)



#endif
