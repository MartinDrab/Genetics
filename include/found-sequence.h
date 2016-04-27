
#ifndef __FOUND_SEQUENCE_H__
#define __FOUND_SEQUENCE_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"



typedef struct _FOUND_SEQUENCE_VARIANT {
	char *Seq1;
	size_t Seq1Len;
	char *Seq2;
	size_t Seq2Len;
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
} FOUND_SEQUENCE, *PFOUND_SEQUENCE;

GEN_ARRAY_TYPEDEF(FOUND_SEQUENCE);
GEN_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE)
POINTER_ARRAY_TYPEDEF(FOUND_SEQUENCE);
POINTER_ARRAY_IMPLEMENTATION(FOUND_SEQUENCE)


ERR_VALUE found_sequence_alloc(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE *FS);
void found_sequence_free(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_init(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE FS);
void found_sequence_finit(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_set_variant(PFOUND_SEQUENCE FS, const size_t Index, PFOUND_SEQUENCE_VARIANT Variant);
boolean found_sequence_match(const FOUND_SEQUENCE *FS, const char *Seq, const size_t Length);


#endif 
