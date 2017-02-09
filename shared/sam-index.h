
#ifndef __READ_SAM_INDEX_H__
#define __READ_SAM_INDEX_H__



#include "utils.h"
#include "gen_dym_array.h"


typedef struct _SAM_INDEX_ENTRY {
	uint64_t FileOffset;
	uint64_t Pos;
	uint64_t ReadIndex;
} SAM_INDEX_ENTRY, *PSAM_INDEX_ENTRY;

GEN_ARRAY_TYPEDEF(SAM_INDEX_ENTRY);
GEN_ARRAY_IMPLEMENTATION(SAM_INDEX_ENTRY);


ERR_VALUE sam_index_load(const char *FileName, PGEN_ARRAY_SAM_INDEX_ENTRY Index);
ERR_VALUE sam_index_save(PGEN_ARRAY_SAM_INDEX_ENTRY Index, const char *FileName);
ERR_VALUE sam_index_from_sam(const char *SAMFile, PGEN_ARRAY_SAM_INDEX_ENTRY Index);
void sam_index_init(PGEN_ARRAY_SAM_INDEX_ENTRY Index);
ERR_VALUE sam_index_insert(PGEN_ARRAY_SAM_INDEX_ENTRY Index, const SAM_INDEX_ENTRY *Entry);
void sam_index_free(PGEN_ARRAY_SAM_INDEX_ENTRY Index);




#endif
