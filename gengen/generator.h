
#ifndef __GENGEN_GENERATOR_H__
#define __GENGEN_GENERATOR_H__

#include <stdint.h>
#include "err.h"
#include "utils.h"


typedef struct _GENERATOR_READ_OPTIONS {
	uint32_t MinReads;
	uint32_t MaxReads;
	uint32_t ReadLength;
} GENERATOR_READ_OPTIONS, *PGENERATOR_READ_OPTIONS;

typedef struct _GENERATOR_INDEL_OPTIONS {
	float IndelProbability;
	uint32_t MinIndels;
	uint32_t MaxIndels;
	boolean DisableIns;
	boolean DisableDels;
} GENERATOR_INDEL_OPTIONS, *PGENERATOR_INDEL_OPTIONS;

typedef struct _GENERATOR_REPLACE_OPTIONS {
	float ReplaceProbability;
	uint32_t MinReplace;
	uint32_t MaxReplace;
} GENERATOR_REPLACE_OPTIONS, *PGENERATOR_REPLACE_OPTIONS;


ERR_VALUE generate_active_region(uint32_t MinLength, uint32_t MaxLength, char *BaseTypes, char **Result, size_t *ResultLength);
void free_active_region(char *Region, const size_t Length);

ERR_VALUE generate_reads(const char *Region, const size_t RegionLength, const char *Nucleotides, const PGENERATOR_READ_OPTIONS ReadOptions, const PGENERATOR_INDEL_OPTIONS IndelOptions, const PGENERATOR_REPLACE_OPTIONS ReplaceOptions, char ***Reads, size_t *ReadCount);
void free_reads(char **Reads, const size_t Count);

#endif 
