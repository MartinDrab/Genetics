
#ifndef __GASSM_INPUT_FILE_H__
#define __GASSM_INPUT_FILE_H__


#include "reads.h"


typedef enum _EActiveRegionType {
	artUnknown,
	artValid,
} EActiveRegionType, *PEActiveRegionType;

typedef struct _ACTIVE_REGION {
	EActiveRegionType Type;
	uint64_t Offset;
	uint64_t Length;
	const char *Sequence;
} ACTIVE_REGION, *PACTIVE_REGION;


ERR_VALUE input_get_refseq(const char *FileName, const char *InputType, char **RefSeq, size_t *RefSeqLen);
void input_free_refseq(char *RefSeq, const size_t RefSeqLen);
ERR_VALUE input_get_reads(const char *Filename, const char *InputType, const uint64_t RegionStart, const uint64_t RegionSize, PONE_READ **Reads, size_t *ReadCount);
void input_free_reads(PONE_READ *Reads, const size_t Count);
ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count);
ERR_VALUE input_get_region_by_offset(const PACTIVE_REGION Regions, const uint32_t Count, const uint64_t Offset, uint32_t *Index, uint64_t *RegionOffset);
void input_free_regions(PACTIVE_REGION Regions, const size_t Count);






#endif 
