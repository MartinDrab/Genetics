
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
ERR_VALUE utils_preallocate_string(const size_t Length, char **Result);
void utils_free_string(char *String);

size_t utils_ranged_rand(const size_t Begin, const size_t End);
boolean utils_prob_happened(const double Probability);

boolean utils_is_prime(const size_t Number);
size_t utils_next_prime(const size_t Number);



#endif
