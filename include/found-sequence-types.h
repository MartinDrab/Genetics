
#ifndef __FOUND_SEQUENCE_TYPES_H__
#define __FOUND_SEQUENCE_TYPES_H__



#include "utils.h"
#include "gen_dym_array.h"
#include "kmer-graph-base-types.h"




typedef enum _EVariantCallPhaseType {
	vcptNone,
	vcptZeroOne,
	vcptZeroTwo,
	vcptOneTwo,
	vcptTwoOne,
	// Transforms to zero-one and zero-two
	vcptBothAlt,
	vcptMax
} EVariantCallPhaseType, *PEVariantCallPhaseType;

typedef struct _VARIANT_CALL {
	char *Chrom;
	uint32_t KMerSize;
	uint32_t BinProb;
	uint32_t ReadCoverage;
	uint64_t Pos;
	char *ID;
	char *Ref;
	char *Alt;
	uint8_t Qual;
	size_t RefWeight;
	size_t AltWeight;
	GEN_ARRAY_size_t RefReads;
	GEN_ARRAY_size_t AltReads;
	boolean Valid;
	uint64_t PhasedPos;
	EVariantCallPhaseType PhaseType;
	void *Context;
} VARIANT_CALL, *PVARIANT_CALL;

GEN_ARRAY_TYPEDEF(VARIANT_CALL);
GEN_ARRAY_IMPLEMENTATION(VARIANT_CALL)
POINTER_ARRAY_TYPEDEF(VARIANT_CALL);
POINTER_ARRAY_IMPLEMENTATION(VARIANT_CALL)








#endif
