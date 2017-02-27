
#ifndef __FOUND_SEQUENCE_TYPES_H__
#define __FOUND_SEQUENCE_TYPES_H__



#include "utils.h"
#include "gen_dym_array.h"
#include "kmer-graph-base-types.h"


typedef struct _FOUND_SEQUENCE_VARIANT {
	uint32_t RefSeqStart;
	uint32_t RefSeqEnd;
	EKMerEdgeType Seq1Type;
	char *Seq1;
	size_t Seq1Len;
	size_t Seq1Weight;
	EKMerEdgeType Seq2Type;
	char *Seq2;
	size_t Seq2Len;
	size_t Seq2Weight;
	GEN_ARRAY_size_t RefReadIndices;
	GEN_ARRAY_size_t ReadIndices;
} FOUND_SEQUENCE_VARIANT, *PFOUND_SEQUENCE_VARIANT;

GEN_ARRAY_TYPEDEF(FOUND_SEQUENCE_VARIANT);
GEN_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE_VARIANT)
POINTER_ARRAY_TYPEDEF(FOUND_SEQUENCE_VARIANT);
POINTER_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE_VARIANT)

typedef struct _FOUND_SEQUENCE {
	char *Sequence;
	size_t Len;
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT Variants;
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT ReadVariants;
} FOUND_SEQUENCE, *PFOUND_SEQUENCE;

GEN_ARRAY_TYPEDEF(FOUND_SEQUENCE);
GEN_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE)
POINTER_ARRAY_TYPEDEF(FOUND_SEQUENCE);
POINTER_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE)

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
} VARIANT_CALL, *PVARIANT_CALL;

GEN_ARRAY_TYPEDEF(VARIANT_CALL);
GEN_ARRAY_IMPLEMENTATION(VARIANT_CALL)
POINTER_ARRAY_TYPEDEF(VARIANT_CALL);
POINTER_ARRAY_IMPLEMENTATION(VARIANT_CALL)








#endif
