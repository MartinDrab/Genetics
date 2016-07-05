
#ifndef __UTILS_H__
#define __UTILS_H__



#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "err.h"



#ifdef _MSC_VER

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
#define INLINE_FUNCTION			static
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

#define FOPEN_MODE_READ				1
#define FOPEN_MODE_WRITE			2
#define FOPEN_MODE_APPEND			4

ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength);

ERR_VALUE utils_fopen(const char *FileName, const uint32_t Mode, FILE **Stream);
ERR_VALUE utils_fread(void *Buffer, const size_t Size, const size_t Count, FILE *Stream);
ERR_VALUE utils_fwrite(const void *Buffer, const size_t Size, const size_t Count, FILE *Stream);
ERR_VALUE utils_fclose(FILE *Stream);


ERR_VALUE _utils_malloc(const size_t Size, void **Address);
ERR_VALUE _utils_calloc(const size_t Count, const size_t Size, void **Address);
void _utils_free(void *Address);

#define ALLOCATOR_HEADER_SIGNATURE					0xbadf00d
#define ALLOCATOR_FOOTER_SIGNATURE					0xf00df00d

typedef struct _ALLOCATOR_HEADER {
	struct _ALLOCATOR_HEADER *Prev;
	struct _ALLOCATOR_HEADER *Next;
	uint32_t Signature;
	const char *Function;
	uint32_t Line;
	size_t BodySize;
	uint32_t Signature2;
} ALLOCATOR_HEADER, *PALLOCATOR_HEADER;

typedef struct _ALLOCATOR_FOOTER {
	uint32_t Signature;
} ALLOCATOR_FOOTER, *PALLOCATOR_FOOTER;

ERR_VALUE _utils_malloc_debug(const size_t Size, void **Address, const char *Function, const uint32_t Line);
ERR_VALUE _utils_calloc_debug(const size_t Count, const size_t Size, void **Address, const char *Function, const uint32_t Line);
void _utils_free_debug(void *Address);
void utils_allocator_check(void);

// #define USE_DEBUG_ALLOCATOR
#ifdef USE_DEBUG_ALLOCATOR

#define utils_malloc(aSize, aAddress)					_utils_malloc_debug((aSize), (aAddress), __FUNCTION__, __LINE__)
#define utils_calloc(aCount, aSize, aAddress)			_utils_calloc_debug((aCount), (aSize), (aAddress), __FUNCTION__, __LINE__)
#define utils_free(aAddress)							_utils_free_debug((aAddress));

#else

#define utils_malloc(aSize, aAddress)					_utils_malloc((aSize), (aAddress))
#define utils_calloc(aCount, aSize, aAddress)			_utils_calloc((aCount), (aSize), (aAddress))
#define utils_free(aAddress)							_utils_free((aAddress));

#endif

#define UTILS_TYPED_MALLOC_FUNCTION(aType)	\
	ERR_VALUE utils_malloc_##aType(aType ** aResult)						\
	{																		\
		return utils_malloc(sizeof(aType), (void **)aResult);				\
	}																		\

#define UTILS_TYPED_CALLOC_FUNCTION(aType)									\
	ERR_VALUE utils_calloc_##aType(const size_t Count, aType ** aResult)	\
	{																		\
		return utils_calloc(Count, sizeof(aType), (void **)aResult);		\
	} \




#endif
