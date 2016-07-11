
#include "err.h"
#include "utils.h"
#include "khash.h"
#include "reads.h"
#include "paired-reads.h"


KHASH_MAP_INIT_STR(RP, PPOINTER_ARRAY_ONE_READ)


static khash_t(RP) *_table = NULL;




ERR_VALUE paired_reads_init(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	_table = kh_init(RP);
	ret = (_table != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;

	return ret;
}


void paired_reads_finit(void)
{
	PPOINTER_ARRAY_ONE_READ pairedReads = NULL;

	for (khiter_t it = kh_begin(_table); it != kh_end(_table); ++it) {
		if (kh_exist(_table, it)) {
			pairedReads = kh_val(_table, it);
			pointer_array_finit_ONE_READ(pairedReads);
			utils_free(pairedReads);
			kh_val(_table, it) = NULL;
		}
	}


	kh_destroy(RP, _table);
	_table = NULL;

	return;
}


ERR_VALUE paired_reads_insert(const ONE_READ *Read)
{
	khiter_t it = kh_end(_table);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PPOINTER_ARRAY_ONE_READ pairedReads = NULL;

	ret = ERR_SUCCESS;
	if (Read->TemplateNameLen > 0) {
		it = kh_get(RP, _table, Read->TemplateName);
		if (it == kh_end(_table)) {
			ret = utils_malloc(sizeof(POINTER_ARRAY_ONE_READ), &pairedReads);
			if (ret == ERR_SUCCESS) {
				int tmp = 0;

				pointer_array_init_ONE_READ(pairedReads, 140);
				it = kh_put(RP, _table, Read->TemplateName, &tmp);
				switch (ret) {
					case -1: ret = ERR_OUT_OF_MEMORY;  break;
					case 0: ret = ERR_ALREADY_EXISTS; break;
					case 2: ret = ERR_NOT_FOUND; break;
				}

				if (ret == ERR_SUCCESS)
					kh_value(_table, it) = pairedReads;
			}
		} else pairedReads = kh_value(_table, it);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_ONE_READ(pairedReads, Read);
	}

	return ret;
}


ERR_VALUE paired_reads_first(khiter_t *Iterator, PPOINTER_ARRAY_ONE_READ *Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	khiter_t it = kh_begin(_table);

	for (; it != kh_end(_table); ++it) {
		if (kh_exist(_table, it))
			break;
	}

	if (it != kh_end(_table)) {
		*Iterator = it;
		*Reads = kh_value(_table, it);
		ret = ERR_SUCCESS;
	} else ret = ERR_NO_MORE_ENTRIES;

	return ret;
}


ERR_VALUE paired_reads_next(khiter_t Iterator, khiter_t *NewIt, PPOINTER_ARRAY_ONE_READ *Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	for (; Iterator != kh_end(_table); ++Iterator) {
		if (kh_exist(_table, Iterator))
			break;
	}

	if (Iterator != kh_end(_table)) {
		*NewIt = Iterator;
		*Reads = kh_value(_table, Iterator);
		ret = ERR_SUCCESS;
	} else ret = ERR_NO_MORE_ENTRIES;

	return ret;
}