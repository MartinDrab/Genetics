
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "err.h"
#include "librandom.h"
#include "utils.h"


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE _utils_malloc(const size_t Size, void **Address)
{
	void *addr = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	addr = malloc(Size);
	if (addr != NULL) {
		*Address = addr;
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

ERR_VALUE _utils_calloc(const size_t Count, const size_t Size, void **Address)
{
	void *addr = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	addr = calloc(Count, Size);
	if (addr != NULL) {
		*Address = addr;
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

void _utils_free(void *Address)
{
	free(Address);

	return;
}

static ALLOCATOR_HEADER _allocHead = {&_allocHead, &_allocHead};


ERR_VALUE _utils_malloc_debug(const size_t Size, void **Address, const char *Function, const uint32_t Line)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	void *addr = NULL;

//	utils_allocator_check();
	ret = _utils_malloc(Size + sizeof(ALLOCATOR_FOOTER) + sizeof(ALLOCATOR_HEADER), &addr);
	if (ret == ERR_SUCCESS) {
		PALLOCATOR_HEADER h = NULL;
		PALLOCATOR_FOOTER f = NULL;

		h = (PALLOCATOR_HEADER)addr;
		addr = (unsigned char *)addr + sizeof(ALLOCATOR_HEADER);
		f = (PALLOCATOR_FOOTER)((unsigned char *)addr + Size);
		h->BodySize = Size;
		h->Function = Function;
		h->Line = Line;
		h->Signature = ALLOCATOR_HEADER_SIGNATURE;
		h->Signature2 = ALLOCATOR_HEADER_SIGNATURE;
		f->Signature = ALLOCATOR_FOOTER_SIGNATURE;
		h->Next = &_allocHead;
		h->Prev = _allocHead.Prev;
		_allocHead.Prev->Next = h;
		_allocHead.Prev = h;
		*Address = addr;
//		utils_allocator_check();
	}

	return ret;
}

ERR_VALUE _utils_calloc_debug(const size_t Count, const size_t Size, void **Address, const char *Function, const uint32_t Line)
{
	return _utils_malloc_debug(Size*Count, Address, Function, Line);
}


void _utils_free_debug(void *Address)
{
	PALLOCATOR_HEADER h = (PALLOCATOR_HEADER)Address - 1;
	PALLOCATOR_FOOTER f = (PALLOCATOR_FOOTER)((unsigned char *)Address + h->BodySize);

//	utils_allocator_check();
	if (h->Signature != ALLOCATOR_HEADER_SIGNATURE)
		__debugbreak();

	if (h->Signature2 != ALLOCATOR_HEADER_SIGNATURE)
		__debugbreak();

	if (f->Signature != ALLOCATOR_FOOTER_SIGNATURE)
		__debugbreak();

	h->Signature = 0;
	h->Signature2 = 0;
	f->Signature = 0;
	h->Prev->Next = h->Next;
	h->Next->Prev = h->Prev;
	_utils_free(h);
//	utils_allocator_check();

	return;
}


void utils_allocator_check(void)
{
	PALLOCATOR_HEADER h = NULL;

	h = _allocHead.Next;
	while (h != &_allocHead) {
		PALLOCATOR_FOOTER f = (PALLOCATOR_FOOTER)((unsigned char *)h + sizeof(ALLOCATOR_HEADER) + h->BodySize);

		if (h->Signature != ALLOCATOR_HEADER_SIGNATURE)
			__debugbreak();

		if (h->Signature2 != ALLOCATOR_HEADER_SIGNATURE)
			__debugbreak();

		if (f->Signature != ALLOCATOR_FOOTER_SIGNATURE)
			__debugbreak();

		h = h->Next;
	}

	return;
}




ERR_VALUE utils_copy_string(const char *String, char **Result)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *tmp = NULL;
	size_t len = ((String != NULL) ? strlen(String) : 0)*sizeof(char);
	
	ret = utils_calloc(len + 1, sizeof(char), &tmp);
	if (ret == ERR_SUCCESS) {
		memcpy(tmp, String, (len + 1)*sizeof(char));
		*Result = tmp;
	}

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
	utils_free(String);

	return;
}

size_t utils_ranged_rand(const size_t Begin, const size_t End)
{
	assert(End > Begin);
	return ((rand_C() % (End - Begin)) + Begin);
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


ERR_VALUE utils_fopen(const char *FileName, const uint32_t Mode, FILE **Stream)
{
	FILE *tmpStream = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *strModes[] = {
		"",
		"rb",
		"wb",
		"wb+"
		"ab",
		"rb",
		"wb",
		"wb+"
	};

#pragma warning(disable : 4996)
	tmpStream = fopen(FileName, strModes[Mode]);
	if (tmpStream != NULL) {
		*Stream = tmpStream;
		ret = ERR_SUCCESS;
	} else ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fread(void *Buffer, const size_t Size, const size_t Count, FILE *Stream)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (fread(Buffer, Size, Count, Stream) != Count)
		ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fwrite(const void *Buffer, const size_t Size, const size_t Count, FILE *Stream)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (fwrite(Buffer, Size, Count, Stream) != Count)
		ret = ERR_IO_ERROR;

	return ret;
}


ERR_VALUE utils_fclose(FILE *Stream)
{
	return (fclose(Stream) == 0) ? ERR_SUCCESS : ERR_IO_ERROR;
}


ERR_VALUE utils_file_read(const char *FileName, char **Data, size_t *DataLength)
{
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_fopen(FileName, FOPEN_MODE_READ, &f);
	if (f != NULL) {
#ifdef _MSC_VER
		ret = _fseeki64(f, 0, SEEK_END);
#else
		ret = fseeko(f, 0, SEEK_END);
#endif
		if (ret == ERR_SUCCESS) {
			off_t fileSize = 0;
#ifdef _MSC_VER
			fileSize = _ftelli64(f);
#else
			fileSize = ftello(f);
#endif
			if (fileSize != -1L) {
#ifdef _MSC_VER
				ret = _fseeki64(f, 0, SEEK_SET);
#else
				ret = fseeko(f, 0, SEEK_SET);
#endif
				if (ret == ERR_SUCCESS) {
					char *tmpData = NULL;
					size_t tmpSize = (size_t)fileSize;

					ret = utils_malloc(tmpSize + sizeof(char), &tmpData);
					if (ret == ERR_SUCCESS) {
						ret = utils_fread(tmpData, 1, tmpSize, f);
						if (ret == ERR_SUCCESS) {
							tmpData[tmpSize] = 0;
							*Data = tmpData;
							*DataLength = tmpSize;
						}

						if (ret != ERR_SUCCESS)
							utils_free(tmpData);
					}
				} else ret = ERR_IO_ERROR;
			} else ret = ERR_IO_ERROR;
		} else ret = ERR_IO_ERROR;

		utils_fclose(f);
	}

	return ret;
}
