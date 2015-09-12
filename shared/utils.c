
#include <stdint.h>
#include <string.h>
#include <malloc.h>
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

void utils_free_string(char *String)
{
	free(String);

	return;
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
