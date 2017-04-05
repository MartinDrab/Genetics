
#ifndef __UTILS_H__
#define __UTILS_H__



#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "err.h"



#ifdef _MSC_VER

#include <windows.h>

#define strcasecmp				_stricmp
#define off_t					long long
#define INLINE_FUNCTION			__inline
#define PATH_SEPARATOR			"\\"

#else 

#include <strings.h>
#undef min
#define min(a, b)				((a) < (b) ? (a) : (b))
#undef max
#define max(a, b)				((a) > (b) ? (a) : (b))
#define INLINE_FUNCTION			inline
#define PATH_SEPARATOR			"/"

#endif


#define flag_on(aValue, aFlag)					(aValue & aFlag)
#define flag_set(aValue, aFlag)					(aValue |= aFlag)
#define flag_clear(aValue, aFlag)				(aValue &= ~(aFlag))

#define in_range(aStart, aSize, aVal)			(((aStart) <= (aVal)) && ((aVal) < ((aStart) + (aSize))))

typedef uint8_t boolean;

#ifndef FALSE
#define FALSE									0
#endif
#ifndef TRUE
#define TRUE									(!FALSE)
#endif

#ifndef min
#define min(a,b)					((a) < (b) ? (a) : (b))
#endif

#define strings_equal(S1, S2)					(strcmp(S1, S2) == 0)

ERR_VALUE utils_copy_string(const char *String, char **Result);
ERR_VALUE utils_preallocate_string(const size_t Length, char **Result);
void utils_free_string(char *String);

size_t utils_ranged_rand(const size_t Begin, const size_t End);
boolean utils_prob_happened(const double Probability);

boolean utils_is_prime(const size_t Number);
size_t utils_next_prime(const size_t Number);
ERR_VALUE utils_mul_inverse(const size_t Number, const size_t Modulus, size_t *Result);
size_t utils_pow_mod(const size_t Base, const size_t Power, const size_t Modulus);

ERR_VALUE _utils_malloc(const size_t Size, void **Address);
ERR_VALUE _utils_calloc(const size_t Count, const size_t Size, void **Address);
void _utils_free(void *Address);

#define ALLOCATOR_HEADER_SIGNATURE					0xbadf00d
#define ALLOCATOR_FOOTER_SIGNATURE					0xf00df00d

typedef struct _ALLOCATOR_FOOTER {
	uint32_t Signature;
} ALLOCATOR_FOOTER, *PALLOCATOR_FOOTER;

typedef struct _ALLOCATOR_HEADER {
	struct _ALLOCATOR_HEADER *Prev;
	struct _ALLOCATOR_HEADER *Next;
	uint32_t Signature;
	const char *Function;
	uint32_t Line;
	size_t BodySize;
	const ALLOCATOR_FOOTER *Footer;
	int ThreadId;
	uint32_t Signature2;
} ALLOCATOR_HEADER, *PALLOCATOR_HEADER;


ERR_VALUE utils_allocator_malloc(const size_t Size, void **Address, const char *Function, const uint32_t Line);
ERR_VALUE utils_allocator_calloc(const size_t Count, const size_t Size, void **Address, const char *Function, const uint32_t Line);
void utils_allocator_free(void *Address);
void utils_allocator_check(void);

#define USE_DEBUG_ALLOCATOR


#define utils_malloc(aSize, aAddress)					_utils_malloc((aSize), (aAddress))
#define utils_calloc(aCount, aSize, aAddress)			_utils_calloc((aCount), (aSize), (aAddress))
#define utils_free(aAddress)							_utils_free((aAddress));
void *_utils_alloc_mark(void);
boolean _utils_alloc_diff(void *Mark);
ERR_VALUE utils_allocator_init(const size_t NumberOfThreads);

#define UTILS_TYPED_MALLOC_FUNCTION(aType)	\
	static ERR_VALUE utils_malloc_##aType(aType ** aResult)						\
	{																		\
		return utils_malloc(sizeof(aType), (void **)aResult);				\
	}																		\

#define UTILS_TYPED_CALLOC_FUNCTION(aType)									\
	INLINE_FUNCTION ERR_VALUE utils_calloc_##aType(const size_t Count, aType ** aResult)	\
	{																		\
		return utils_calloc(Count, sizeof(aType), (void **)aResult);		\
	} \

#define UTILS_NAMED_CALLOC_FUNCTION(aName, aType)									\
	INLINE_FUNCTION ERR_VALUE utils_calloc_##aName(const size_t Count, aType ** aResult)	\
	{																		\
		return utils_calloc(Count, sizeof(aType), (void **)aResult);		\
	} \


UTILS_TYPED_CALLOC_FUNCTION(uint8_t)
UTILS_TYPED_CALLOC_FUNCTION(char)
UTILS_TYPED_CALLOC_FUNCTION(size_t)
UTILS_NAMED_CALLOC_FUNCTION(puint8_t, uint8_t *)


INLINE_FUNCTION long utils_atomic_increment(long volatile *Data)
{
#ifdef _MSC_VER
	return InterlockedIncrement(Data);
#else
	return __sync_add_and_fetch(Data, 1);
#endif
}


INLINE_FUNCTION long utils_atomic_decrement(long volatile *Data)
{
#ifdef _MSC_VER
	return InterlockedDecrement(Data);
#else
	return __sync_sub_and_fetch(Data, 1);
#endif
}

INLINE_FUNCTION long utils_atomic_exchange(long volatile *Data, long Value)
{
#ifdef _MSC_VER
	return InterlockedExchange(Data, Value);
#else
	return __sync_lock_test_and_set(Data, Value);
#endif
}

INLINE_FUNCTION void utils_atomic_release(long volatile *Lock)
{
#ifdef _MSC_VER
	InterlockedExchange(Lock, 0);
#else
	__sync_lock_release(Lock);
#endif
}

#ifdef WIN32




double erand48(unsigned short xseed[3]);
double drand48();
void srand48(long seed);


#endif


#endif
