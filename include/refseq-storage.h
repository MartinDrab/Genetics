
#ifndef __REFSEQ_STORAGE_H__
#define __REFSEQ_STORAGE_H__



#include "err.h"
#include "utils.h"



typedef struct _REFSEQ_STORAGE {
	size_t ValidLength;
	size_t AllocatedLength;
	char *Sequence;
	char Storage[100];
} REFSEQ_STORAGE, *PREFSEQ_STORAGE;



INLINE_FUNCTION void rs_storage_init(PREFSEQ_STORAGE Storage)
{
	Storage->ValidLength = 0;
	Storage->AllocatedLength = sizeof(Storage->Storage) / sizeof(Storage->Storage[0]);
	Storage->Sequence = Storage->Storage;

	return;
}

INLINE_FUNCTION void rs_storage_finit(PREFSEQ_STORAGE Storage)
{
	if (Storage->Sequence != Storage->Storage)
		utils_free(Storage->Sequence);

	return;
}


INLINE_FUNCTION ERR_VALUE rs_storage_add(PREFSEQ_STORAGE Storage, char Base)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Storage->ValidLength < Storage->AllocatedLength) {
		Storage->Sequence[Storage->ValidLength] = Base;
		++Storage->ValidLength;
		ret = ERR_SUCCESS;
	} else {
		char *tmp = NULL;

		ret = utils_calloc(Storage->AllocatedLength * 2, sizeof(char), (void **)&tmp);
		if (ret == ERR_SUCCESS) {
			memcpy(tmp, Storage->Sequence, Storage->ValidLength*sizeof(char));
			tmp[Storage->ValidLength] = Base;
			if (Storage->Sequence != Storage->Storage)
				utils_free(Storage->Sequence);

			Storage->Sequence = tmp;
			Storage->AllocatedLength *= 2;
			++Storage->ValidLength;
		}

	}

	return ret;
}


INLINE_FUNCTION ERR_VALUE rs_storage_add_seq(PREFSEQ_STORAGE Storage, const char *Bases, size_t Length)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Storage->ValidLength + Length <= Storage->AllocatedLength) {
		memcpy(Storage->Sequence + Storage->ValidLength, Bases, Length*sizeof(char));
		Storage->ValidLength += Length;
		ret = ERR_SUCCESS;
	} else {
		char *tmp = NULL;
		size_t newAllocLength = max(Storage->AllocatedLength * 2, Storage->ValidLength + Length + 100);

		ret = utils_calloc(newAllocLength, sizeof(char), (void **)&tmp);
		if (ret == ERR_SUCCESS) {
			memcpy(tmp, Storage->Sequence, Storage->ValidLength*sizeof(char));
			memcpy(tmp + Storage->ValidLength, Bases, Length*sizeof(char));
			if (Storage->Sequence != Storage->Storage)
				utils_free(Storage->Sequence);

			Storage->Sequence = tmp;
			Storage->AllocatedLength = newAllocLength;
			Storage->ValidLength += Length;
		}

	}

	return ret;
}


INLINE_FUNCTION void rs_storage_remove(PREFSEQ_STORAGE Storage, size_t Length)
{
	Storage->ValidLength -= Length;

	return;
}


INLINE_FUNCTION rs_storage_reset(PREFSEQ_STORAGE Storage)
{
	Storage->ValidLength = 0;

	return;
}


INLINE_FUNCTION ERR_VALUE rs_storage_create_string(const REFSEQ_STORAGE *Storage, char **String)
{
	char *tmpString = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(Storage->ValidLength + 1, sizeof(char), (void **)&tmpString);
	if (ret == ERR_SUCCESS) {
		memcpy(tmpString, Storage->Sequence, Storage->ValidLength*sizeof(char));
		tmpString[Storage->ValidLength] = '\0';
		*String = tmpString;
	}

	return ret;
}


#endif
