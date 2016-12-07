
#include "err.h"
#include "utils.h"
#include "utils-lookaside.h"





ERR_VALUE utils_lookaside_init(PUTILS_LOOKASIDE Lookaside, const size_t BlockSize, const size_t Count)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	Lookaside->BlockSize = BlockSize;
	Lookaside->PreallocatedCount = Count;
	ret = utils_calloc(Count, BlockSize, &Lookaside->PreallocatedArea);
	if (ret == ERR_SUCCESS) {
		PUTILS_LOOKASIDE_BLOCK block = Lookaside->PreallocatedArea;
		
		for (size_t i = 0; i < Count - 1; ++i) {
			block->Next = (PUTILS_LOOKASIDE_BLOCK)((unsigned char *)block + BlockSize);
			block = block->Next;
		}

		block->Next = NULL;
		Lookaside->Block = Lookaside->PreallocatedArea;
	}

	return ret;
}


void utils_lookaside_finit(PUTILS_LOOKASIDE Lookaside)
{
	PUTILS_LOOKASIDE_BLOCK block = NULL;
	PUTILS_LOOKASIDE_BLOCK old = NULL;

	block = Lookaside->Block;
	while (block != NULL) {
		old = block;
		block = block->Next;
		if (!in_range(Lookaside->PreallocatedArea, Lookaside->BlockSize*Lookaside->PreallocatedCount, old))
			utils_free(old);
	}

	utils_free(Lookaside->PreallocatedArea);

	return;
}

ERR_VALUE utils_lookaside_alloc(PUTILS_LOOKASIDE Lookaside, void **Block)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Lookaside->Block != NULL) {
		*Block = Lookaside->Block;
		Lookaside->Block = Lookaside->Block->Next;
	} else ret = utils_malloc(Lookaside->BlockSize, Block);

	return ret;
}


void utils_lookaside_free(PUTILS_LOOKASIDE Lookaside, void *Block)
{
	PUTILS_LOOKASIDE_BLOCK block = (PUTILS_LOOKASIDE_BLOCK)Block;

	block->Next = Lookaside->Block;
	Lookaside->Block = block;

	return;
}
