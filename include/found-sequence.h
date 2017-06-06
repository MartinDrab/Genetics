
#ifndef __FOUND_SEQUENCE_H__
#define __FOUND_SEQUENCE_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "found-sequence-types.h"
#include "kmer-graph-types.h"



ERR_VALUE variant_call_init(const char *Chrom, uint64_t Pos, const char *ID, const char *Ref, size_t RefLen, const char *Alt, size_t AltLen, const uint8_t Qual, const GEN_ARRAY_size_t *RefReads, const GEN_ARRAY_size_t *AltReads, PVARIANT_CALL VC);
void variant_call_finit(PVARIANT_CALL VC);
boolean variant_call_equal(const VARIANT_CALL *VC1, const VARIANT_CALL *VC2);
ERR_VALUE vc_array_add(PGEN_ARRAY_VARIANT_CALL Array, const VARIANT_CALL *VC, PVARIANT_CALL *Existing);
void vc_array_clear(PGEN_ARRAY_VARIANT_CALL Array);
void vc_array_finit(PGEN_ARRAY_VARIANT_CALL Array);
void vc_array_print(FILE *Stream, const char *ReferenceFile, const GEN_ARRAY_VARIANT_CALL *Array);
void vc_array_sort(PGEN_ARRAY_VARIANT_CALL Array);
ERR_VALUE vc_array_merge(PGEN_ARRAY_VARIANT_CALL Dest, PGEN_ARRAY_VARIANT_CALL Sources, const size_t SourceCount);
void vc_array_map_to_edges(PGEN_ARRAY_VARIANT_CALL VCArray);
ERR_VALUE vc_array_intersection(GEN_ARRAY_VARIANT_CALL *A1, const GEN_ARRAY_VARIANT_CALL *A2, PGEN_ARRAY_VARIANT_CALL Intersection);


#endif 
