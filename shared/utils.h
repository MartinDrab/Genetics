
#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdint.h>
#include <string.h>
#include "err.h"


#define flag_on(aValue, aFlag)					(aValue & aFlag)
#define flag_set(aValue, aFlag)					(aValue |= aFlag)
#define flag_clear(aValue, aFlag)				(aValue &= ~(aFlag))

typedef uint8_t boolean;

#define FALSE									0
#define TRUE									(!FALSE)


#define strings_equal(S1, S2)					(strcmp(S1, S2) == 0)

ERR_VALUE utils_copy_string(char *String, char **Result);
void utils_free_string(char *String);

boolean utils_is_prime(const size_t Number);



#endif
