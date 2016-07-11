
#ifndef __FOUND_SEQUENCE_H__
#define __FOUND_SEQUENCE_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "kmer-graph-types.h"


typedef struct _FOUND_SEQUENCE_VARIANT {
	uint32_t RefSeqStart;
	uint32_t RefSeqEnd;
	EKMerEdgeType Seq1Type;
	char *Seq1;
	size_t Seq1Len;
	long Seq1Weight;
	EKMerEdgeType Seq2Type;
	char *Seq2;
	size_t Seq2Len;
	long Seq2Weight;
	size_t Reserved;
	const char *LastFPos;
	const char *LastSPos;
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

typedef struct _VARIANT_CALL {
	char *Chrom;
	uint64_t Pos;
	char *ID;
	char *Ref;
	char *Alt;
	uint8_t Qual;
	long RefWeight;
	long AltWeight;
} VARIANT_CALL, *PVARIANT_CALL;

GEN_ARRAY_TYPEDEF(VARIANT_CALL);
GEN_ARRAY_IMPLEMENTATION(VARIANT_CALL)
POINTER_ARRAY_TYPEDEF(VARIANT_CALL);
POINTER_ARRAY_IMPLEMENTATION(VARIANT_CALL)


ERR_VALUE found_sequence_alloc(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE *FS);
void found_sequence_free(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_init(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE FS);
void found_sequence_finit(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_set_variant(PFOUND_SEQUENCE FS, const size_t Index, PFOUND_SEQUENCE_VARIANT Variant);
boolean found_sequence_match(const FOUND_SEQUENCE *FS, const char *Seq, const size_t Length);
ERR_VALUE found_sequence_build_read_variants(PFOUND_SEQUENCE FS, const POINTER_ARRAY_KMER_EDGE *PathEdges);

ERR_VALUE variant_call_init(const char *Chrom, const uint64_t Pos, const char *ID, const char *Ref, const char *Alt, const uint8_t Qual, PVARIANT_CALL VC);
void variant_call_finit(PVARIANT_CALL VC);
boolean variant_call_equal(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2);
ERR_VALUE vc_array_add(PGEN_ARRAY_VARIANT_CALL Array, const VARIANT_CALL *VC);
void vc_array_finit(PGEN_ARRAY_VARIANT_CALL Array);
void vc_array_print(FILE *Stream, const GEN_ARRAY_VARIANT_CALL *Array);
void vc_array_sort(PGEN_ARRAY_VARIANT_CALL Array);
ERR_VALUE vc_array_merge(PGEN_ARRAY_VARIANT_CALL Dest, PGEN_ARRAY_VARIANT_CALL Sources, const size_t SourceCount);


#endif 
