
#include "utils.h"
#include "file-utils.h"
#include "reads.h"
#include "sam-index.h"



static int _sam_index_comp(const SAM_INDEX_ENTRY *E1, const SAM_INDEX_ENTRY *E2)
{
	if (E1->Pos < E2->Pos)
		return -1;

	if (E1->Pos > E2->Pos)
		return 1;

	return 0;
}


ERR_VALUE sam_index_load(const char *FileName, PGEN_ARRAY_SAM_INDEX_ENTRY Index)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *data = NULL;
	size_t len = 0;

	sam_index_init(Index);
	ret = utils_file_read(FileName, &data, &len);
	if (ret == ERR_SUCCESS) {
		size_t entryCount = len / sizeof(SAM_INDEX_ENTRY);

		if (len % sizeof(SAM_INDEX_ENTRY) == 0) {
			ret = dym_array_reserve_SAM_INDEX_ENTRY(Index, entryCount);
			if (ret == ERR_SUCCESS) {
				memcpy(Index->Data, data, len);
				Index->ValidLength = entryCount;
			}			;
		} else ret = ERR_INTERNAL_ERROR;

		utils_free(data);
	}

	if (ret != ERR_SUCCESS)
		sam_index_free(Index);

	return ret;
}


ERR_VALUE sam_index_save(PGEN_ARRAY_SAM_INDEX_ENTRY Index, const char *FileName)
{
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_fopen(FileName, FOPEN_MODE_WRITE, &f);
	if (ret == ERR_SUCCESS) {
		ret = utils_fwrite(Index->Data, sizeof(SAM_INDEX_ENTRY), gen_array_size(Index), f);
		utils_fclose(f);
	}

	return ret;
}


ERR_VALUE sam_index_from_sam(const char *SAMFile, PGEN_ARRAY_SAM_INDEX_ENTRY Index)
{
	const char *line = NULL;
	const char *lineEnd = NULL;
	size_t len = 0;
	char *data = NULL;
	ONE_READ r;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_file_read(SAMFile, &data, &len);
	if (ret == ERR_SUCCESS) {
		sam_index_init(Index);
		line = data;
		while (ret == ERR_SUCCESS) {
			ret = read_create_from_sam_line(line, &r);
			if (ret == ERR_SUCCESS) {
				SAM_INDEX_ENTRY e;

				e.FileOffset = line - data;
				e.Pos = r.Pos;
				e.ReadIndex = 0;
				ret = sam_index_insert(Index, &e);
				_read_destroy_structure(&r);
			}
		}

		utils_free(data);
	}

	if (ret == ERR_SUCCESS)
		qsort(Index->Data, gen_array_size(Index), sizeof(SAM_INDEX_ENTRY), _sam_index_comp);

	return ret;
}


void sam_index_init(PGEN_ARRAY_SAM_INDEX_ENTRY Index)
{
	dym_array_init_SAM_INDEX_ENTRY(Index, 140);

	return;
}

ERR_VALUE sam_index_insert(PGEN_ARRAY_SAM_INDEX_ENTRY Index, const SAM_INDEX_ENTRY *Entry)
{
	return dym_array_push_back_SAM_INDEX_ENTRY(Index, *Entry);
}


void sam_index_free(PGEN_ARRAY_SAM_INDEX_ENTRY Index)
{
	dym_array_finit_SAM_INDEX_ENTRY(Index);

	return;
}
