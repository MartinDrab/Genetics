
#ifndef __UTILS_H__
#define __UTILS_H__



#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "err.h"



#ifdef _MSC_VER

#define strcasecmp				_stricmp

#else 

#include <strings.h>
#undef min
#define min(a, b)				((a) < (b) ? (a) : (b))
#undef max
#define max(a, b)				((a) > (b) ? (a) : (b))
#endif


#define flag_on(aValue, aFlag)					(aValue & aFlag)
#define flag_set(aValue, aFlag)					(aValue |= aFlag)
#define flag_clear(aValue, aFlag)				(aValue &= ~(aFlag))

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

ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength);



ERR_VALUE utils_malloc(const size_t Size, void **Address);
ERR_VALUE utils_calloc(const size_t Count, const size_t Size, void **Address);
void utils_free(void *Address);


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
