
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "err.h"
#include "utils.h"


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE utils_copy_string(const char *String, char **Result)
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

ERR_VALUE utils_mul_inverse(const size_t Number, const size_t Modulus, size_t *Result)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Number < Modulus) {
		ret = ERR_SUCCESS;
		size_t v = 1;
		size_t d = Number;
		size_t u = (Number == 1);
		size_t t = 1 - u;
		if (t == 1) {
			size_t c = Modulus % Number;

			if (c != 0) {
				u = (Modulus / Number); // floor()
				while (c != 1 && t == 1) {
					if (c == 0) {
						ret = ERR_NO_INVERSE;
						break;
					}

					size_t q = (d / c); // floor
					d = d % c;
					v = v + q*u;
					t = (d != 1);
					if (t == 1) {
						q = (c / d); // floor
						c = c % d;
						u = u + q*v;
					}
				}

				u = v*(1 - t) + t*(Modulus - u);
			}
			else ret = ERR_NO_INVERSE;
		}

		if (ret == ERR_SUCCESS)
			*Result = u;
	}
	else ret = ERR_MODULUS_TOO_SMALL;

	return ret;
}

size_t utils_pow_mod(const size_t Base, const size_t Power, const size_t Modulus)
{
	size_t ret = Base;

	assert(Power > 0);
	for (size_t i = 1; i < Power; ++i) {
		ret *= Base;
		ret %= Modulus;
	}

	return ret;
}

ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength)
{
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	f = fopen(FileName, "rb");
	if (f != NULL) {
		ret = fseek(f, 0, SEEK_END);
		if (ret == ERR_SUCCESS) {
			long fileSize = ftell(f);

			if (fileSize != -1L) {
				ret = fseek(f, 0, SEEK_SET);
				if (ret == ERR_SUCCESS) {
					char *tmpData = NULL;
					size_t tmpSize = (size_t)fileSize;

					tmpData = (char *)malloc(tmpSize);
					if (tmpData != NULL) {
						if (fread(tmpData, 1, tmpSize, f) == tmpSize) {
							*Data = tmpData;
							*DataLength = tmpSize;
							ret = ERR_SUCCESS;
						} else ret = ERR_FERROR;

						if (ret != ERR_SUCCESS)
							free(tmpData);
					} else ret = ERR_OUT_OF_MEMORY;
				} else ret = ERR_INTERNAL_ERROR;
			} else ret = ERR_ERRNO_VALUE;
		} else ret = ERR_INTERNAL_ERROR;

		fclose(f);
	}

	return ret;
}
