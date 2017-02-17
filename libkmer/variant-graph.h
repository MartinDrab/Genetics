
#ifndef __VARIANT_GRAPH_H__
#define __VARIANT_GRAPH_H__



#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "khash.h"
#include "found-sequence-types.h"


typedef enum _EVariantGraphVertexColor {
	vgvcNone,
	vgvcSequence1,
	vgvcSequence2,
	vgvcBoth,
	vgvcMax,
} EVariantGraphVertexColo, *PEVariantGraphVertexColo;

GEN_ARRAY_TYPEDEF(EVariantGraphVertexColo);
GEN_ARRAY_IMPLEMENTATION(EVariantGraphVertexColo)

typedef enum _EVariangGraphVertexType {
	vgvtReference,
	vgvtAlternative,
	vgvtMax,
} EVariangGraphVertexType, *PEVariangGraphVertexType;

struct _VARIANT_GRAPH_VERTEX;

typedef struct _VARIANT_GRAPH_PAIRED_EDGE {
	struct _VARIANT_GRAPH_VERTEX *Target;
	size_t Count;
} VARIANT_GRAPH_PAIRED_EDGE, *PVARIANT_GRAPH_PAIRED_EDGE;

typedef struct _VARIANT_GRAPH_VERTEX {
	EVariangGraphVertexType Type;
	PVARIANT_CALL Variant;
	const size_t *ReadIndices;
	size_t ReadCount;
	PVARIANT_GRAPH_PAIRED_EDGE Paired;
	size_t PairedCount;
	size_t ComponentIndex;
	size_t Index;
	EVariantGraphVertexColo Color;
	boolean Uncolorable;
} VARIANT_GRAPH_VERTEX, *PVARIANT_GRAPH_VERTEX;

POINTER_ARRAY_TYPEDEF(VARIANT_GRAPH_VERTEX);
POINTER_ARRAY_IMPLEMENTATION(VARIANT_GRAPH_VERTEX)
GEN_ARRAY_TYPEDEF(VARIANT_GRAPH_VERTEX);
GEN_ARRAY_IMPLEMENTATION(VARIANT_GRAPH_VERTEX)

typedef struct _VARIANG_GRAPH_READ_EDGES {
	size_t *RefToRef;
	size_t *RefToAlt;
	size_t *AltToRef;
	size_t *AltToAlt;
} VARIANG_GRAPH_READ_EDGES, *PVARIANG_GRAPH_READ_EDGES;

KHASH_MAP_INIT_INT64(ReadToVertex, PPOINTER_ARRAY_VARIANT_GRAPH_VERTEX)

typedef struct _VARIANT_GRAPH_THRESHOLDS {
	size_t Read;
	size_t Paired;
} VARIANT_GRAPH_THRESHOLDS, *PVARIANT_GRAPH_THRESHOLDS;

typedef struct _VARIANT_GRAPH {
	union {
		struct {
			PVARIANT_GRAPH_VERTEX Reference;
			PVARIANT_GRAPH_VERTEX Alternative;
		} ByType;
		PVARIANT_GRAPH_VERTEX All[vgvtMax];
		PVARIANT_GRAPH_VERTEX OneArray;
	} Vertices;
	size_t VerticesArraySize;
	union {
		VARIANG_GRAPH_READ_EDGES ByTypes;
		size_t *All[4];
	} ReadEdges;
	khash_t(ReadToVertex) *ReadMap;
	GEN_ARRAY_size_t ComponentIndices;
	POINTER_ARRAY_VARIANT_GRAPH_VERTEX Components;
	VARIANT_GRAPH_THRESHOLDS Thresholds;
} VARIANT_GRAPH, *PVARIANT_GRAPH;



ERR_VALUE vg_graph_init(PVARIANT_CALL Variants, const size_t VariantCount, size_t Threshold, PVARIANT_GRAPH Graph);
void vg_graph_finit(PVARIANT_GRAPH Graph);
ERR_VALUE vg_graph_add_paired(PVARIANT_GRAPH Graph);
ERR_VALUE vg_graph_color(PVARIANT_GRAPH Graph);
void vg_graph_print(FILE *Stream, const VARIANT_GRAPH *Graph);
void vg_graph_finalize(PVARIANT_GRAPH Graph);


#endif
