
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

/** Describes one variant found in the de Bruijn graph. */
typedef struct _VARIANT_CALL {
	/** Chromosome. */
	char *Chrom;
	/** K-mer size used by de Bruijn graph that discovered the variant. */
	uint32_t KMerSize;
	/** The binomial test coefficient, in percents. */
	uint32_t BinProb;
	/** Variant starting position. */
	uint64_t Pos;
	/** Always - */
	char *ID;
	/** The reference part. */
	char *Ref;
	/** The alternate part. */
	char *Alt;
	/** Always 60. */
	uint8_t Qual;
	/** Not used. */
	size_t RefWeight;
	/** Not used. */
	size_t AltWeight;
	/** Indices of reads covering the reference part. */
	GEN_ARRAY_size_t RefReads;
	/** Indices of reads covering the alternate one. */
	GEN_ARRAY_size_t AltReads;
	/** Indicates whether the variant will be included in the output VCF file. */
	boolean Valid;
	/** Identifier of the variant graph component (used also for phasing). */
	uint64_t PhasedPos;
	/** Phasing type. */
	EVariantCallPhaseType PhaseType;
	void *Context;
} VARIANT_CALL, *PVARIANT_CALL;

GEN_ARRAY_TYPEDEF(VARIANT_CALL);
GEN_ARRAY_IMPLEMENTATION(VARIANT_CALL)
POINTER_ARRAY_TYPEDEF(VARIANT_CALL);
POINTER_ARRAY_IMPLEMENTATION(VARIANT_CALL)








#endif
