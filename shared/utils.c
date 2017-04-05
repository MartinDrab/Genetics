
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#ifdef WIN32
#include <windows.h>
#endif
#include "err.h"
#include "utils.h"


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE _utils_malloc(const size_t Size, void **Address)
{
	void *addr = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	addr = malloc(Size);
	*Address = addr;
	if (addr != NULL) {
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

ERR_VALUE _utils_calloc(const size_t Count, const size_t Size, void **Address)
{
	void *addr = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	addr = calloc(Count, Size);
	*Address = addr;
	if (addr != NULL) {
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

void _utils_free(void *Address)
{
	free(Address);

	return;
}

static ALLOCATOR_HEADER *_allocators = NULL;
#ifdef WIN32
static HANDLE *_allocatorHeaps = NULL;
#endif

ERR_VALUE utils_allocator_init(const size_t NumberOfThreads)
{
	ERR_VALUE ret = ERR_OUT_OF_MEMORY;

	_allocators = (ALLOCATOR_HEADER *)calloc(NumberOfThreads, sizeof(ALLOCATOR_HEADER));
	if (_allocators != NULL) {
		for (size_t i = 0; i < NumberOfThreads; ++i) {
			_allocators[i].Next = _allocators + i;
			_allocators[i].Prev = _allocators + i;
		}

		ret = ERR_SUCCESS;
#ifdef WIN32
		_allocatorHeaps = calloc(NumberOfThreads, sizeof(HANDLE));
		if (_allocatorHeaps != NULL) {
			for (size_t i = 0; i < NumberOfThreads; ++i) {
				_allocatorHeaps[i] = HeapCreate(HEAP_NO_SERIALIZE, 0, 0);
				if (_allocatorHeaps[i] == NULL) {
					ret = ERR_OUT_OF_MEMORY;
					for (size_t j = 0; j < i; ++j)
						HeapDestroy(_allocatorHeaps[j]);

					break;
				}
			}

			if (ret != ERR_SUCCESS)
				free(_allocatorHeaps);
		}
#endif

		if (ret != ERR_SUCCESS)
			free(_allocators);
	}

	return ret;
}


void *_utils_alloc_mark(void)
{
	return _allocators[omp_get_thread_num()].Prev;
}

boolean _utils_alloc_diff(void *Mark)
{
	boolean ret = FALSE;
	PALLOCATOR_HEADER item = (PALLOCATOR_HEADER)Mark;
	PALLOCATOR_HEADER head = _allocators + omp_get_thread_num();

	item = item->Next;
	ret = (item != head);
	while (item != head) {
		fprintf(stderr, "[LEAK]: %zu bytes, function %s, line %u\n", item->BodySize, item->Function, item->Line);
		item = item->Next;
	}

	return ret;
}

ERR_VALUE utils_allocator_malloc(const size_t Size, void **Address, const char *Function, const uint32_t Line)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	void *addr = NULL;
	PALLOCATOR_HEADER head = _allocators + omp_get_thread_num();
#ifdef USE_DEBUG_ALLOCATOR
	const size_t realSize = Size + sizeof(ALLOCATOR_HEADER) + sizeof(ALLOCATOR_FOOTER);
#else
	const size_t realSize = Size;
#endif

	addr = malloc(realSize);
	ret = (addr != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
	if (ret == ERR_SUCCESS) {
#ifdef USE_DEBUG_ALLOCATOR
		PALLOCATOR_HEADER h = NULL;
		PALLOCATOR_FOOTER f = NULL;

		h = (PALLOCATOR_HEADER)addr;
		addr = (unsigned char *)addr + sizeof(ALLOCATOR_HEADER);
		f = (PALLOCATOR_FOOTER)((unsigned char *)addr + Size);
		h->BodySize = Size;
		h->Function = Function;
		h->Line = Line;
		h->Footer = f;
		h->ThreadId = omp_get_thread_num();
		h->Signature = ALLOCATOR_HEADER_SIGNATURE;
		h->Signature2 = ALLOCATOR_HEADER_SIGNATURE;
		f->Signature = ALLOCATOR_FOOTER_SIGNATURE;
		h->Next = head;
		h->Prev = head->Prev;
		head->Prev->Next = h;
		head->Prev = h;
#endif
		*Address = addr;
	}

	return ret;
}


ERR_VALUE utils_allocator_calloc(const size_t Count, const size_t Size, void **Address, const char *Function, const uint32_t Line)
{
	return utils_allocator_malloc(Size*Count, Address, Function, Line);
}


void utils_allocator_free(void *Address)
{
#ifdef USE_DEBUG_ALLOCATOR
	PALLOCATOR_HEADER h = (PALLOCATOR_HEADER)Address - 1;
	const ALLOCATOR_FOOTER *f = h->Footer;

//	assert(h->ThreadId == omp_get_thread_num());
	assert(h->Signature == ALLOCATOR_HEADER_SIGNATURE);
	assert(h->Signature2 == ALLOCATOR_HEADER_SIGNATURE);
	assert(f->Signature == ALLOCATOR_FOOTER_SIGNATURE);

	h->Prev->Next = h->Next;
	h->Next->Prev = h->Prev;
	memset(h, 0xbadf00d, h->BodySize + sizeof(ALLOCATOR_HEADER) + sizeof(ALLOCATOR_FOOTER));
	Address = h;
#endif
	free(Address);

	return;
}


void utils_allocator_check(void)
{
	const ALLOCATOR_HEADER *h = NULL;
	PALLOCATOR_HEADER head = _allocators + omp_get_thread_num();

	h = head->Next;
	while (h != head) {
		const ALLOCATOR_FOOTER *f = h->Footer;

		assert(h->ThreadId == omp_get_thread_num());
		assert(h->BodySize == (unsigned char *)f - (unsigned char *)h - sizeof(ALLOCATOR_HEADER));
		assert(h->Signature == ALLOCATOR_HEADER_SIGNATURE);
		assert(h->Signature2 == ALLOCATOR_HEADER_SIGNATURE);
		assert(f->Signature == ALLOCATOR_FOOTER_SIGNATURE);
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
		memcpy(tmp, String, len*sizeof(char));
		tmp[len] = '\0';
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

