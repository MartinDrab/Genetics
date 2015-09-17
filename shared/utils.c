
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "utils.h"


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE utils_copy_string(char *String, char **Result)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *tmp = NULL;
	size_t len = ((String != NULL) ? strlen(String) : 0)*sizeof(char);
	
	tmp = malloc(len + sizeof(char));
	if (tmp != NULL) {
		memcpy(tmp, String, len + sizeof(char));
		*Result = tmp;
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

ERR_VALUE utils_preallocate_string(const size_t Length, char **Result)
{
	char *tmpResult = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	tmpResult = (char *)malloc((Length + 1)*sizeof(char));
	if (tmpResult != NULL) {
		tmpResult[Length] = '\0';
		*Result = tmpResult;
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

void utils_free_string(char *String)
{
	free(String);

	return;
}

size_t utils_ranged_rand(const size_t Begin, const size_t End)
{
	assert(End > Begin);
	return ((rand() % (End - Begin)) + Begin);
}


boolean utils_is_prime(const size_t Number)
{
	boolean ret = FALSE;
	size_t max = (size_t)sqrt((double)Number);

	ret = (Number & 1);
	if (ret) {
		size_t test = 3;

		while (test <= max) {
			ret = (Number % test != 0);
			if (!ret)
				break;

			test += 2;
		}
	}

	return ret;
}

size_t utils_next_prime(const size_t Number)
{
	size_t next = Number;

	next += (Number % 2 == 0) ? 1 : 2;

	while (!utils_is_prime(next))
		next += 2;

	return next;
}

boolean utils_prob_happened(const double Probability)
{
	return (((double)rand() / (double)RAND_MAX) < Probability);
}