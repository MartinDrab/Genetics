
#ifndef __FOUND_SEQUENCE_H__
#define __FOUND_SEQUENCE_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "found-sequence-types.h"
#include "kmer-graph-types.h"



void found_sequence_variant_free(const FOUND_SEQUENCE_VARIANT *FSV);
void found_sequence_variant_array_free(const FOUND_SEQUENCE_VARIANT *Array, const size_t Count);
ERR_VALUE found_sequence_variant_copy(PFOUND_SEQUENCE_VARIANT Target, const FOUND_SEQUENCE_VARIANT *Source);

ERR_VALUE found_sequence_alloc(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE *FS);
void found_sequence_free(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_init(const char *Sequence, const size_t Length, const size_t VariantCount, PFOUND_SEQUENCE FS);
void found_sequence_finit(PFOUND_SEQUENCE FS);
ERR_VALUE found_sequence_set_variant(PFOUND_SEQUENCE FS, const size_t Index, PFOUND_SEQUENCE_VARIANT Variant);
boolean found_sequence_match(const FOUND_SEQUENCE *FS, const char *Seq, const size_t Length);
ERR_VALUE found_sequence_build_read_variants(PFOUND_SEQUENCE FS, const POINTER_ARRAY_KMER_EDGE *PathEdges);

ERR_VALUE variant_call_init(const char *Chrom, const uint64_t Pos, const char *ID, const char *Ref, size_t RefLen, const char *Alt, size_t AltLen, const uint8_t Qual, PVARIANT_CALL VC);
void variant_call_finit(PVARIANT_CALL VC);
boolean variant_call_equal(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2);
ERR_VALUE vc_array_add(PGEN_ARRAY_VARIANT_CALL Array, const VARIANT_CALL *VC);
void vc_array_finit(PGEN_ARRAY_VARIANT_CALL Array);
void vc_array_print(FILE *Stream, const GEN_ARRAY_VARIANT_CALL *Array);
void vc_array_sort(PGEN_ARRAY_VARIANT_CALL Array);
ERR_VALUE vc_array_merge(PGEN_ARRAY_VARIANT_CALL Dest, PGEN_ARRAY_VARIANT_CALL Sources, const size_t SourceCount);


#endif 
