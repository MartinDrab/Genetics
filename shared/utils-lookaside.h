
#ifndef __UTILS_LOOKASIDE_H__
#define __UTILS_LOOKASIDE_H__


#include "err.h"
#include "utils.h"


typedef struct _UTILS_LOOKASIDE_BLOCK {
	struct _UTILS_LOOKASIDE_BLOCK *Next;
} UTILS_LOOKASIDE_BLOCK, *PUTILS_LOOKASIDE_BLOCK;

typedef struct _UTILS_LOOKASIDE {
	size_t BlockSize;
	size_t PreallocatedCount;
	boolean Initialized;
	PUTILS_LOOKASIDE_BLOCK PreallocatedArea;
	PUTILS_LOOKASIDE_BLOCK Block;
} UTILS_LOOKASIDE, *PUTILS_LOOKASIDE;

UTILS_TYPED_CALLOC_FUNCTION(PUTILS_LOOKASIDE)


ERR_VALUE utils_lookaside_init(PUTILS_LOOKASIDE Lookaside, const size_t BlockSize, const size_t Count);
void utils_lookaside_finit(PUTILS_LOOKASIDE Lookaside);
ERR_VALUE utils_lookaside_alloc(PUTILS_LOOKASIDE Lookaside, void **Block);
void utils_lookaside_free(PUTILS_LOOKASIDE Lookaside, void *Block);



#endif
