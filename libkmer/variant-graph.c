
#include <omp.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "khash.h"
#include "found-sequence.h"
#include "paired-reads.h"
#include "variant-graph.h"



/************************************************************************/
/*                       HELPER FUNCTIONS                               */
/************************************************************************/


UTILS_TYPED_CALLOC_FUNCTION(PPOINTER_ARRAY_VARIANT_GRAPH_VERTEX)
UTILS_TYPED_CALLOC_FUNCTION(VARIANT_GRAPH_VERTEX)
UTILS_TYPED_CALLOC_FUNCTION(VARIANT_GRAPH_PAIRED_EDGE)

#define _vg_vertex_exists(aGraph, aVertex)	\
	((aVertex)->ReadCount > ((aGraph)->Thresholds.Read))

#define _vg_vertex_index(aGraph, aVertex)	\
	(size_t)(((aVertex)->Type == vgvtReference) ? ((aVertex) - (aGraph)->Vertices.ByType.Reference) : ((aVertex) - (aGraph)->Vertices.ByType.Alternative))

#define _vg_read_edge_exists(aGraph, aSourceType, aDestType, aSourceIndex)	\
	((aGraph)->ReadEdges.All[((aSourceType) << 1) + (aDestType)][(aSourceIndex)] > (aGraph)->Thresholds.Read)

#define _vg_vertex_get(aGraph, aType, aIndex)	\
	((aGraph)->Vertices.All[(aType)] + (aIndex))

#define _vg_vertex_get_variant(aGraph, aVertex)	\
	(_vg_vertex_get(aGraph, (aVertex)->Type ^ 1, _vg_vertex_index(aGraph, aVertex)))

#define _vg_opposite_color(aColor)	\
	((aColor) ^ 3)

static size_t _intersection_size(const size_t *A, const size_t ACount, const size_t *B, const size_t BCount)
{
	size_t ret = 0;
	const size_t *currentA = A;
	const size_t *currentB = B;
	size_t AIndex = 0;
	size_t BIndex = 0;

	while (AIndex < ACount && BIndex < BCount) {
		if (*currentA == *currentB) {
			++AIndex;
			++BIndex;
			++currentA;
			++currentB;
			++ret;
		} else if (*currentA < *currentB) {
			++currentA;
			++AIndex;
		} else {
			++currentB;
			++BIndex;
		}
	}

	return ret;
}


static ERR_VALUE _update_read_map(khash_t(ReadToVertex) *Map, const VARIANT_GRAPH_VERTEX *Vertex)
{
	khiter_t it;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t *ri = Vertex->ReadIndices;
	PPOINTER_ARRAY_VARIANT_GRAPH_VERTEX va = NULL;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < Vertex->ReadCount; ++i) {
		it = kh_get(ReadToVertex, Map, *ri);
		if (it == kh_end(Map)) {
			ret = utils_malloc(sizeof(POINTER_ARRAY_VARIANT_GRAPH_VERTEX), &va);
			if (ret == ERR_SUCCESS) {
				int tmp;

				pointer_array_init_VARIANT_GRAPH_VERTEX(va, 140);
				it = kh_put(ReadToVertex, Map, *ri, &tmp);
				if (tmp != -1)
					kh_val(Map, it) = va;
				else ret = ERR_OUT_OF_MEMORY;

				if (ret != ERR_SUCCESS)
					pointer_array_finit_VARIANT_GRAPH_VERTEX(va);
			}
		} else va = kh_val(Map, it);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(va, Vertex);

		if (ret != ERR_SUCCESS)
			break;

		++ri;
	}

	return ret;
}


static void _vg_vertex_init(PVARIANT_CALL Variant, const size_t Index, const EVariangGraphVertexType Type, PVARIANT_GRAPH_VERTEX Vertex)
{
	const GEN_ARRAY_size_t *indicesArray = NULL;

	Vertex->Index = Index;
	Vertex->Variant = Variant;
	Vertex->Type = Type;
	switch (Vertex->Type) {
		case vgvtReference:
			indicesArray = &Variant->RefReads;
			Vertex->Weight = Variant->RefWeight;
			break;
		case vgvtAlternative:
			indicesArray = &Variant->AltReads;
			Vertex->Weight = Variant->AltWeight;
			break;
		default:
			assert(FALSE);
			break;
	}

	Vertex->ReadCount = gen_array_size(indicesArray);
	Vertex->ReadIndices = indicesArray->Data;
	Vertex->Paired = NULL;
	Vertex->PairedCount = 0;
	Vertex->Color = vgvcNone;
	Vertex->Uncolorable = FALSE;
	Vertex->ComponentIndex = 0;

	return;
}


static void _vg_vertex_finit(PVARIANT_GRAPH_VERTEX Vertex)
{
	if (Vertex->Paired != NULL)
		utils_free(Vertex->Paired);

	return;
}


static int _pvg_vertex_comparator(const void *PV1, const void *PV2)
{
	int ret = 0;
	const VARIANT_GRAPH_VERTEX *v1 = *(const PVARIANT_GRAPH_VERTEX *)PV1;
	const VARIANT_GRAPH_VERTEX *v2 = *(const PVARIANT_GRAPH_VERTEX *)PV2;

	if (v1->Index < v2->Index)
		ret = -1;
	else if (v2->Index < v1->Index)
		ret = 1;
	else {
		if (v1->Type < v2->Type)
			ret = -1;
		else if (v2->Type < v1->Type)
			ret = 1;	
	}

	return ret;
}

static ERR_VALUE _compute_components(PVARIANT_GRAPH Graph)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PVARIANT_GRAPH_VERTEX v = NULL;
	PVARIANT_GRAPH_VERTEX var = NULL;
	size_t currentComponent = 1;
	POINTER_ARRAY_VARIANT_GRAPH_VERTEX propagationStack;

	ret = ERR_SUCCESS;
	dym_array_init_size_t(&Graph->ComponentIndices, 140);
	pointer_array_init_VARIANT_GRAPH_VERTEX(&propagationStack, 140);
	pointer_array_init_VARIANT_GRAPH_VERTEX(&Graph->Components, 140);
	ret = pointer_array_reserve_VARIANT_GRAPH_VERTEX(&Graph->Components, Graph->VerticesArraySize * 2);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < Graph->VerticesArraySize * 2; ++i) {
			v = Graph->Vertices.OneArray + i;
			if (_vg_vertex_exists(Graph, v) && v->ComponentIndex == 0) {
				ret = dym_array_push_back_size_t(&Graph->ComponentIndices, pointer_array_size(&Graph->Components));
				if (ret == ERR_SUCCESS)
					ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, v);
				
				while (pointer_array_size(&propagationStack) > 0) {
					v = *pointer_array_pop_back_VARIANT_GRAPH_VERTEX(&propagationStack);
					assert(v->ComponentIndex == 0 || v->ComponentIndex == currentComponent);
					if (_vg_vertex_exists(Graph, v) && v->ComponentIndex == 0) {
						v->ComponentIndex = currentComponent;
						pointer_array_push_back_no_alloc_VARIANT_GRAPH_VERTEX(&Graph->Components, v);

						//Variant edge
						var = _vg_vertex_get_variant(Graph, v);
						assert(var->ComponentIndex == 0 || var->ComponentIndex == currentComponent);
						if (_vg_vertex_exists(Graph, var) && var->ComponentIndex == 0)
							ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, var);

						// Read edges
						size_t vIndex = _vg_vertex_index(Graph, v);
						if (vIndex < Graph->VerticesArraySize - 1)
						{
							PVARIANT_GRAPH_VERTEX nextRef = _vg_vertex_get(Graph, vgvtReference, vIndex + 1);
							PVARIANT_GRAPH_VERTEX nextAlt = _vg_vertex_get(Graph, vgvtAlternative, vIndex + 1);

							if (_vg_read_edge_exists(Graph, v->Type, nextRef->Type, vIndex) && _vg_vertex_exists(Graph, nextRef)) {
								assert(nextRef->ComponentIndex == 0 || nextRef->ComponentIndex == currentComponent);
								if (nextRef->ComponentIndex == 0)
									ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, nextRef);
							}

							if (_vg_read_edge_exists(Graph, v->Type, nextAlt->Type, vIndex) && _vg_vertex_exists(Graph, nextAlt)) {
								assert(nextAlt->ComponentIndex == 0 || nextAlt->ComponentIndex == currentComponent);
								if (nextAlt->ComponentIndex == 0)
									ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, nextAlt);
							}
						}

						if (vIndex > 0) {
							PVARIANT_GRAPH_VERTEX prevRef = _vg_vertex_get(Graph, vgvtReference, vIndex - 1);
							PVARIANT_GRAPH_VERTEX prevAlt = _vg_vertex_get(Graph, vgvtAlternative, vIndex - 1);

							if (_vg_read_edge_exists(Graph, prevRef->Type, v->Type, vIndex - 1) && _vg_vertex_exists(Graph, prevRef)) {
								assert(prevRef->ComponentIndex == 0 || prevRef->ComponentIndex == currentComponent);
								if (prevRef->ComponentIndex == 0)
									ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, prevRef);
							}

							if (_vg_read_edge_exists(Graph, prevAlt->Type, v->Type, vIndex - 1) && _vg_vertex_exists(Graph, prevAlt)) {
								assert(prevAlt->ComponentIndex == 0 || prevAlt->ComponentIndex == currentComponent);
								if (prevAlt->ComponentIndex == 0)
									ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, prevAlt);
							}
						}

						// paired edges
						for (size_t j = 0; j < v->PairedCount; ++j) {
							PVARIANT_GRAPH_VERTEX p = v->Paired[j].Target;

							if (_vg_vertex_exists(Graph, p)) {
								assert(p->ComponentIndex == 0 || p->ComponentIndex == currentComponent);
								if (p->ComponentIndex == 0)
									ret = pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, p);
							}

							if (ret != ERR_SUCCESS)
								break;
						}
					}
				}

				++currentComponent;
			}

			if (ret != ERR_SUCCESS)
				break;
		}

		if (ret == ERR_SUCCESS)
			ret = dym_array_push_back_size_t(&Graph->ComponentIndices, pointer_array_size(&Graph->Components));
	}

	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < gen_array_size(&Graph->ComponentIndices) - 1; ++i) {
			const size_t startIndex = Graph->ComponentIndices.Data[i];
			const size_t endIndex = Graph->ComponentIndices.Data[i + 1];

			qsort(Graph->Components.Data + startIndex, endIndex - startIndex, sizeof(PVARIANT_GRAPH_VERTEX), _pvg_vertex_comparator);
		}
	}

	if (ret != ERR_SUCCESS) {
		dym_array_finit_size_t(&Graph->ComponentIndices);
		pointer_array_finit_VARIANT_GRAPH_VERTEX(&Graph->Components);
	}

	pointer_array_finit_VARIANT_GRAPH_VERTEX(&propagationStack);

	return ret;
}

typedef struct _COLOR_DECISION {
	PVARIANT_GRAPH_VERTEX Vertex;
	EVariantGraphVertexColo Color;
	boolean Propagation;
	size_t Index;
} COLOR_DECISION, *PCOLOR_DECISION;

GEN_ARRAY_TYPEDEF(COLOR_DECISION);
GEN_ARRAY_IMPLEMENTATION(COLOR_DECISION)


ERR_VALUE vg_graph_color_component(PVARIANT_GRAPH Graph, const size_t Start, const size_t End)
{
	size_t index = Start;
	GEN_ARRAY_COLOR_DECISION stack;
	const size_t totalVertices = End - Start;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PVARIANT_GRAPH_VERTEX v = NULL;
	PVARIANT_GRAPH_VERTEX var = NULL;
	boolean colors[vgvcMax];
	POINTER_ARRAY_VARIANT_GRAPH_VERTEX propagationStack;

	pointer_array_init_VARIANT_GRAPH_VERTEX(&propagationStack, 140);
	ret = pointer_array_reserve_VARIANT_GRAPH_VERTEX(&propagationStack, 4*totalVertices);
	if (ret == ERR_SUCCESS) {
		dym_array_init_COLOR_DECISION(&stack, 140);
		ret = dym_array_reserve_COLOR_DECISION(&stack, totalVertices);
		if (ret == ERR_SUCCESS) {
			EVariantGraphVertexColo color = vgvcNone;

			do {
				boolean propagation = pointer_array_size(&propagationStack) > 0;

				if (propagation)
					v = *pointer_array_pop_back_VARIANT_GRAPH_VERTEX(&propagationStack);
				else {
					if (index < End)
						v = Graph->Components.Data[index];
					else v = NULL;
				}

				if (v != NULL && _vg_vertex_exists(Graph, v) && v->Color == vgvcNone) {
					memset(colors, FALSE, sizeof(colors));
					colors[vgvcNone] = TRUE;
					colors[vgvcBoth] = TRUE;
					for (EVariantGraphVertexColo c = vgvcNone; c < color; ++c)
						colors[c] = TRUE;

					// Paired edges
					for (size_t i = 0; i < v->PairedCount; ++i) {
						PVARIANT_GRAPH_VERTEX p = v->Paired[i].Target;

						if (_vg_vertex_exists(Graph, p)) {
							assert(p->ComponentIndex == v->ComponentIndex);
							switch (p->Color) {
								case vgvcSequence1:
								case vgvcSequence2:
									colors[_vg_opposite_color(p->Color)] = TRUE;
									break;
								case vgvcNone:
									if (!p->Uncolorable)
										pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, p);
									break;
								default:
									break;
							}
						}
					}

					// Incoming read edges
					if (_vg_vertex_index(Graph, v) > 0) {
						const size_t baseIndex = _vg_vertex_index(Graph, v) - 1;
						const VARIANT_GRAPH_VERTEX *tmp = NULL;

						if (_vg_read_edge_exists(Graph, vgvtReference, v->Type, baseIndex)) {
							tmp = _vg_vertex_get(Graph, vgvtReference, baseIndex);
							if (_vg_vertex_exists(Graph, tmp)) {
								assert(tmp->ComponentIndex == v->ComponentIndex);
								switch (tmp->Color) {
									case vgvcSequence1:
									case vgvcSequence2:
										colors[_vg_opposite_color(tmp->Color)] = TRUE;
										break;
									case vgvcNone:
										if (!tmp->Uncolorable)
											pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, tmp);
										break;
									default:
										break;
								}
							}
						}

						if (_vg_read_edge_exists(Graph, vgvtAlternative, v->Type, baseIndex)) {
							tmp = _vg_vertex_get(Graph, vgvtAlternative, baseIndex);
							if (_vg_vertex_exists(Graph, tmp)) {
								assert(tmp->ComponentIndex == v->ComponentIndex);
								switch (tmp->Color) {
									case vgvcSequence1:
									case vgvcSequence2:
										colors[_vg_opposite_color(tmp->Color)] = TRUE;
										break;
									case vgvcNone:
										if (!tmp->Uncolorable)
											pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, tmp);
										break;
									default:
										break;
								}
							}
						}
					}

					// Outgoing read edges
					if (_vg_vertex_index(Graph, v) < Graph->VerticesArraySize - 1)
					{
						const size_t baseIndex = _vg_vertex_index(Graph, v);
						const VARIANT_GRAPH_VERTEX *tmp = NULL;

						if (_vg_read_edge_exists(Graph, vgvtReference, v->Type, baseIndex)) {
							tmp = _vg_vertex_get(Graph, vgvtReference, baseIndex + 1);
							if (_vg_vertex_exists(Graph, tmp)) {
								assert(tmp->ComponentIndex == v->ComponentIndex);
								switch (tmp->Color) {
									case vgvcSequence1:
									case vgvcSequence2:
										colors[_vg_opposite_color(tmp->Color)] = TRUE;
										break;
									case vgvcNone:
										if (!tmp->Uncolorable)
											pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, tmp);
										break;
									default:
										break;
								}
							}
						}

						if (_vg_read_edge_exists(Graph, vgvtAlternative, v->Type, baseIndex)) {
							tmp = _vg_vertex_get(Graph, vgvtAlternative, baseIndex + 1);
							if (_vg_vertex_exists(Graph, tmp)) {
								assert(tmp->ComponentIndex == v->ComponentIndex);
								switch (tmp->Color) {
									case vgvcSequence1:
									case vgvcSequence2:
										colors[_vg_opposite_color(tmp->Color)] = TRUE;
										break;
									case vgvcNone:
										if (!tmp->Uncolorable)
											pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, tmp);
										break;
									default:
										break;
								}
							}
						}
					}

					// Variant edge
					var = _vg_vertex_get_variant(Graph, v);
					if (_vg_vertex_exists(Graph, var)) {
						assert(var->ComponentIndex == v->ComponentIndex);
						switch (var->Color) {
							case vgvcSequence1:
							case vgvcSequence2:
								colors[var->Color] = TRUE;
								break;
							case vgvcNone:
								if (!var->Uncolorable)
									pointer_array_push_back_VARIANT_GRAPH_VERTEX(&propagationStack, var);
								break;
							default:
								break;
						}
					}

					color = vgvcNone;
					while (color < vgvcMax && colors[color])
						++color;

					if (color < vgvcMax) {
						COLOR_DECISION d;
						
						v->Color = color;
						d.Color = color;
						d.Vertex = v;
						d.Propagation = propagation;
						d.Index = index;
						dym_array_push_back_no_alloc_COLOR_DECISION(&stack, d);
						assert(gen_array_size(&stack) <= totalVertices);
						color = vgvcNone;
					} else {
						v->Uncolorable = TRUE;
						ret = ERR_CANNOT_COLOR;
						/*
						do {
							--decision;
							index = decision->Index;
							color = decision->Color;
							v = decision->Vertex;
							v->Color = vgvcNone;
						} while (decision->Propagation);

						--index;
						propagation = FALSE;
						*/
					}
				}

				if (!propagation)
					++index;
			} while (pointer_array_size(&propagationStack) > 0 || index < End);
		}

		dym_array_finit_COLOR_DECISION(&stack);
	}

	pointer_array_finit_VARIANT_GRAPH_VERTEX(&propagationStack);

	return ret;
}


static void _delete_vertex(PVARIANT_GRAPH Graph, PVARIANT_GRAPH_VERTEX Vertex)
{
	PVARIANT_GRAPH_VERTEX p = NULL;
	size_t index = _vg_vertex_index(Graph, Vertex);

	Vertex->ReadCount = 0;
	if (index > 0) {
		--index;
		switch (Vertex->Type) {
			case vgvtReference:
				Graph->ReadEdges.ByTypes.RefToRef[index] = 0;
				Graph->ReadEdges.ByTypes.AltToRef[index] = 0;
				Graph->ReadEdges.ByTypes.RefToRef[index + 1] = 0;
				Graph->ReadEdges.ByTypes.RefToAlt[index + 1] = 0;
				break;
			case vgvtAlternative:
				Graph->ReadEdges.ByTypes.RefToAlt[index] = 0;
				Graph->ReadEdges.ByTypes.AltToAlt[index] = 0;
				Graph->ReadEdges.ByTypes.AltToRef[index + 1] = 0;
				Graph->ReadEdges.ByTypes.AltToAlt[index + 1] = 0;
				break;
			default:
				break;
		}
	}

	for (size_t j = 0; j < Vertex->PairedCount; ++j) {
		p = Vertex->Paired[j].Target;
		for (size_t k = 0; k < p->PairedCount; ++k) {
			if (p->Paired[k].Target == Vertex) {
				while (p->PairedCount > k && p->Paired[k].Target == Vertex) {
					p->Paired[k] = p->Paired[p->PairedCount - 1];
					p->PairedCount--;
				}
			}
		}
	}

	Vertex->PairedCount = 0;

	return;
}


static boolean vg_graph_remove_paired_colisions(PVARIANT_GRAPH Graph, const size_t Start, const size_t End)
{
	size_t index = 0;
	PVARIANT_GRAPH_VERTEX v = NULL;
	PVARIANT_GRAPH_VERTEX p = NULL;
	PVARIANT_GRAPH_VERTEX var = NULL;
	PVARIANT_GRAPH_VERTEX *pv = Graph->Components.Data + Start;
	boolean ret = FALSE;

	for (size_t i = Start; i < End; ++i) {
		v = *pv;
		if (_vg_vertex_exists(Graph, v) && v->Color == vgvcNone && v->Type == vgvtAlternative) {
			boolean pairedSeq1 = FALSE;
			boolean pairedSeq2 = FALSE;
			
			for (size_t j = 0; j < v->PairedCount; ++j) {
				p = v->Paired[j].Target;
				pairedSeq1 |= p->Color == vgvcSequence1;
				pairedSeq2 |= p->Color == vgvcSequence2;
			}

			if (pairedSeq1 && pairedSeq2) {
				var = _vg_vertex_get_variant(Graph, v);
				if (_vg_vertex_exists(Graph, var)) {
					pairedSeq1 = FALSE;
					pairedSeq2 = FALSE;
					for (size_t j = 0; j < v->PairedCount; ++j) {
						p = v->Paired[j].Target;
						pairedSeq1 |= p->Color == vgvcSequence1;
						pairedSeq2 |= p->Color == vgvcSequence2;
					}

					if (!pairedSeq1 || !pairedSeq2) {
						_delete_vertex(Graph, v);
						v->Color = vgvcBoth;
						ret = TRUE;
					}
				}
			}
		}

		++pv;
	}

	return ret;
}


ERR_VALUE vg_graph_color(PVARIANT_GRAPH Graph)
{
	int i = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t componentCount = gen_array_size(&Graph->ComponentIndices);

	ret = ERR_SUCCESS;
#pragma omp parallel for shared(Graph)
	for (i = 0; i < (int)componentCount - 1; ++i) {
		const size_t startIndex = Graph->ComponentIndices.Data[i];
		const size_t endIndex = Graph->ComponentIndices.Data[i + 1];

		ret = vg_graph_color_component(Graph, startIndex, endIndex);
		if (ret == ERR_CANNOT_COLOR) {
			while (vg_graph_remove_paired_colisions(Graph, startIndex, endIndex))
				ret = vg_graph_color_component(Graph, startIndex, endIndex);
			
			ret = ERR_SUCCESS;
		}
	}

	return ret;
}



static void _optimize_both_paths(PVARIANT_GRAPH Graph)
{
	PVARIANT_GRAPH_VERTEX v = NULL;
	PVARIANT_GRAPH_VERTEX var = NULL;
	
	for (EVariangGraphVertexType type = vgvtReference; type < vgvtMax; ++type) {
		for (size_t i = 0; i < Graph->VerticesArraySize; ++i) {
			v = _vg_vertex_get(Graph, type, i);
			if (_vg_vertex_exists(Graph, v)) {
				boolean refExists = _vg_read_edge_exists(Graph, v->Type, vgvtReference, i);
				boolean altExists = _vg_read_edge_exists(Graph, v->Type, vgvtAlternative, i);
				boolean varExists = FALSE;

				var = _vg_vertex_get_variant(Graph, v);
				varExists = _vg_vertex_exists(Graph, var);
				if (!varExists /*|| (refExists && altExists)*/) {
					v->Color = vgvcBoth;
					if (varExists)
						_delete_vertex(Graph, var);
				}
			}
		}
	}
	
	return;
}


/************************************************************************/
/*                    PUBLIC FUNCTIONS                                  */
/************************************************************************/


ERR_VALUE vg_graph_init(PVARIANT_CALL Variants, const size_t VariantCount, size_t Threshold, PVARIANT_GRAPH Graph)
{
	int i = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	Graph->Thresholds.Paired = Threshold;
	Graph->Thresholds.Read = Threshold;
	Graph->ReadMap = kh_init(ReadToVertex);
	if (Graph->ReadMap != NULL) {
		Graph->VerticesArraySize = VariantCount;
		ret = utils_calloc_VARIANT_GRAPH_VERTEX(2 * VariantCount, &Graph->Vertices.ByType.Reference);
		if (ret == ERR_SUCCESS) {
			Graph->Vertices.ByType.Alternative = Graph->Vertices.ByType.Reference + VariantCount;
			ret = utils_calloc_size_t(4 * VariantCount, &Graph->ReadEdges.ByTypes.RefToRef);
			if (ret == ERR_SUCCESS) {
				memset(Graph->ReadEdges.ByTypes.RefToRef, 0, 4 * sizeof(size_t));
				Graph->ReadEdges.ByTypes.RefToAlt = Graph->ReadEdges.ByTypes.RefToRef + VariantCount;
				Graph->ReadEdges.ByTypes.AltToRef = Graph->ReadEdges.ByTypes.RefToAlt + VariantCount;
				Graph->ReadEdges.ByTypes.AltToAlt = Graph->ReadEdges.ByTypes.AltToRef + VariantCount;
#pragma omp parallel for shared(Graph)
				for (i = 0; i < (int)VariantCount; ++i) {
					_vg_vertex_init(Variants + i, i, vgvtReference, Graph->Vertices.ByType.Reference + i);
					_vg_vertex_init(Variants + i, i, vgvtAlternative, Graph->Vertices.ByType.Alternative + i);
				}

#pragma omp parallel for shared(Graph)
				for (i = 0; i < (int)VariantCount - 1; ++i) {
					const VARIANT_GRAPH_VERTEX *u = Graph->Vertices.ByType.Reference + i;
					const VARIANT_GRAPH_VERTEX *v = Graph->Vertices.ByType.Reference + i + 1;

					Graph->ReadEdges.ByTypes.RefToRef[i] = _intersection_size(u->ReadIndices, u->ReadCount, v->ReadIndices, v->ReadCount);
					v = Graph->Vertices.ByType.Alternative + i + 1;
					Graph->ReadEdges.ByTypes.RefToAlt[i] = _intersection_size(u->ReadIndices, u->ReadCount, v->ReadIndices, v->ReadCount);
					u = Graph->Vertices.ByType.Alternative + i;
					v = Graph->Vertices.ByType.Reference + i + 1;
					Graph->ReadEdges.ByTypes.AltToRef[i] = _intersection_size(u->ReadIndices, u->ReadCount, v->ReadIndices, v->ReadCount);
					v = Graph->Vertices.ByType.Alternative + i + 1;
					Graph->ReadEdges.ByTypes.AltToAlt[i] = _intersection_size(u->ReadIndices, u->ReadCount, v->ReadIndices, v->ReadCount);
				}

				for (size_t j = 0; j < vgvtMax; ++j) {
					for (size_t i = 0; i < VariantCount; ++i) {
						ret = _update_read_map(Graph->ReadMap, Graph->Vertices.All[j] + i);
						if (ret != ERR_SUCCESS)
							break;
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				if (ret == ERR_SUCCESS)
					_optimize_both_paths(Graph);

				if (ret != ERR_SUCCESS)
					utils_free(Graph->ReadEdges.ByTypes.RefToRef);
			}

			if (ret != ERR_SUCCESS)
				utils_free(Graph->Vertices.ByType.Reference);
		}
	
		if (ret != ERR_SUCCESS) {
			for (khiter_t it = kh_begin(Graph->ReadMap); it != kh_end(Graph->ReadMap); ++it) {
				if (kh_exist(Graph->ReadMap, it)) {
					pointer_array_finit_VARIANT_GRAPH_VERTEX(kh_val(Graph->ReadMap, it));
				}
			}

			kh_destroy(ReadToVertex, Graph->ReadMap);
		}
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}


void vg_graph_print(FILE *Stream, const VARIANT_GRAPH *Graph)
{
	const VARIANT_CALL *vc = NULL;
	const char *colors[] = {
		"yellow",
		"red",
		"blue",
		"purple",
	};
	const char *typeNames[] = {
		"Ref", 
		"Alt"
	};
	size_t vertexCounts[vgvcMax];
	size_t totalVertexCount = 0;
	size_t typeCounts[vgvtMax];

	memset(vertexCounts, 0, sizeof(vertexCounts));
	memset(typeCounts, 0, sizeof(typeCounts));
	for (size_t i = 0; i < pointer_array_size(&Graph->Components); ++i) {
		const VARIANT_GRAPH_VERTEX *v = Graph->Components.Data[i];

		if (_vg_vertex_exists(Graph, v)) {
			++vertexCounts[v->Color];
			++typeCounts[v->Type];
			++totalVertexCount;
		}
	}

	fprintf(Stream, "digraph G {\n");
	fprintf(Stream, "\t/* Total number of vertices: %zu */\n", totalVertexCount);
	fprintf(Stream, "\t/* Number of reference vertices: %zu */\n", typeCounts[vgvtReference]);
	fprintf(Stream, "\t/* Number of alternate vertices: %zu */\n", typeCounts[vgvtAlternative]);
	fprintf(Stream, "\t/* Number of non-colored vertices: %zu */\n", vertexCounts[vgvcNone]);
	fprintf(Stream, "\t/* Number of seq1-colored vertices: %zu */\n", vertexCounts[vgvcSequence1]);
	fprintf(Stream, "\t/* Number of seq2-colored vertices: %zu */\n", vertexCounts[vgvcSequence2]);
	fprintf(Stream, "\t/* Number of both-colored vertices: %zu */\n", vertexCounts[vgvcBoth]);
	for (size_t i = 0; i < Graph->VerticesArraySize; ++i) {
		const VARIANT_GRAPH_VERTEX *v;
			
		vc = Graph->Vertices.ByType.Reference[i].Variant;
		v = Graph->Vertices.ByType.Reference + i;
		if (_vg_vertex_exists(Graph, v))
			fprintf(Stream, "\tV%zu_1[label=\"Pos: %" PRId64 "\\nCount: %zu; Weight: %zu\\nType: %s\",style=filled,color=%s];\n", i, vc->Pos, v->ReadCount, v->Weight, typeNames[v->Type], colors[v->Color]);

		v = Graph->Vertices.ByType.Alternative + i;
		if (_vg_vertex_exists(Graph, v))
			fprintf(Stream, "\tV%zu_2[label=\"Pos: %" PRId64 "\\nCount: %zu;Weight: %zu\\nType: %s\",style=filled,color=%s];\n", i, vc->Pos, v->ReadCount, v->Weight, typeNames[v->Type], colors[v->Color]);

		if (_vg_vertex_exists(Graph, Graph->Vertices.ByType.Reference + i) && 
			_vg_vertex_exists(Graph, Graph->Vertices.ByType.Alternative + i))
			fprintf(Stream, "\tV%zu_1 -> V%zu_2 [arrowhead=none,color=blue];\n", i, i);
	}

	if (Graph->VerticesArraySize > 0) {
		for (size_t i = 0; i < Graph->VerticesArraySize - 1; ++i) {
			if (_vg_read_edge_exists(Graph, vgvtReference, vgvtReference, i))
				fprintf(Stream, "\tV%zu_1 -> V%zu_1[label=\"%zu\",color=green];\n", i, i + 1, Graph->ReadEdges.ByTypes.RefToRef[i]);

			if (_vg_read_edge_exists(Graph, vgvtReference, vgvtAlternative, i))
				fprintf(Stream, "\tV%zu_1 -> V%zu_2[label=\"%zu\",color=red];\n", i, i + 1, Graph->ReadEdges.ByTypes.RefToAlt[i]);

			if (_vg_read_edge_exists(Graph, vgvtAlternative, vgvtReference, i))
				fprintf(Stream, "\tV%zu_2 -> V%zu_1[label=\"%zu\",color=red];\n", i, i + 1, Graph->ReadEdges.ByTypes.AltToRef[i]);

			if (_vg_read_edge_exists(Graph, vgvtAlternative, vgvtAlternative, i))
				fprintf(Stream, "\tV%zu_2 -> V%zu_2[label=\"%zu\",color=red];\n", i, i + 1, Graph->ReadEdges.ByTypes.AltToAlt[i]);
		}
	}

	for (size_t i = 0; i < Graph->VerticesArraySize; ++i) {
		size_t pindex = 0;
		PVARIANT_GRAPH_VERTEX v = NULL;
		PVARIANT_GRAPH_VERTEX p = NULL;
		size_t c = 0;
		v = Graph->Vertices.ByType.Reference + i;
		for (size_t j = 0; j < v->PairedCount; ++j) {
			c = v->Paired[j].Count;
			p = v->Paired[j].Target;
			pindex = _vg_vertex_index(Graph, p);
			if (_vg_vertex_exists(Graph, p)) {
				switch (p->Type) {
					case vgvtReference:
				 		if (pindex > i || (pindex == i && p->Type != v->Type))
							fprintf(Stream, "\tV%zu_1 -> V%zu_1[arrowhead=none,color=yellow,label=\"%zu\"];\n", i, pindex, c);
						break;
					case vgvtAlternative:
						if (pindex > i || (pindex == i && p->Type != v->Type))
							fprintf(Stream, "\tV%zu_1 -> V%zu_2[arrowhead=none,color=yellow,label=\"%zu\"];\n", i, pindex, c);
						break;
					default:
						break;
				}
			}
		}

		v = Graph->Vertices.ByType.Alternative + i;
		for (size_t j = 0; j < v->PairedCount; ++j) {
			c = v->Paired[j].Count;
			p = v->Paired[j].Target;
			pindex = _vg_vertex_index(Graph, p);
			if (_vg_vertex_exists(Graph, p)) {
				switch (p->Type) {
				case vgvtReference:
					if (pindex > i)
						fprintf(Stream, "\tV%zu_2 -> V%zu_1[arrowhead=none,color=yellow,label=\"%zu\"];\n", i, pindex, c);
					break;
				case vgvtAlternative:
					if (pindex > i)
						fprintf(Stream, "\tV%zu_2 -> V%zu_2[arrowhead=none,color=yellow,label=\"%zu\"];\n", i, pindex, c);
					break;
				default:
					break;
				}
			}
		}
	}

	fprintf(Stream, "}\n");

	return;
}


void vg_graph_finit(PVARIANT_GRAPH Graph)
{
	pointer_array_finit_VARIANT_GRAPH_VERTEX(&Graph->Components);
	dym_array_finit_size_t(&Graph->ComponentIndices);
	utils_free(Graph->ReadEdges.ByTypes.RefToRef);
	for (size_t i = 0; i < Graph->VerticesArraySize; ++i) {
		_vg_vertex_finit(Graph->Vertices.ByType.Reference + i);
		_vg_vertex_finit(Graph->Vertices.ByType.Alternative + i);
	}

	utils_free(Graph->Vertices.ByType.Reference);
	for (khiter_t it = kh_begin(Graph->ReadMap); it != kh_end(Graph->ReadMap); ++it) {
		if (kh_exist(Graph->ReadMap, it))
			pointer_array_finit_VARIANT_GRAPH_VERTEX(kh_val(Graph->ReadMap, it));
	}

	kh_destroy(ReadToVertex, Graph->ReadMap);

	return;
}


ERR_VALUE vg_graph_add_paired(PVARIANT_GRAPH Graph)
{
	PPOINTER_ARRAY_VARIANT_GRAPH_VERTEX *readVertices = NULL;
	khiter_t pairedIter;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PPOINTER_ARRAY_ONE_READ pairedReads = NULL;

	ret = paired_reads_first(&pairedIter, &pairedReads);
	while (ret == ERR_SUCCESS) {
		const size_t prCount = pointer_array_size(pairedReads);
		
		ret = utils_calloc_PPOINTER_ARRAY_VARIANT_GRAPH_VERTEX(prCount, &readVertices);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < prCount; ++i) {
				khiter_t mapIt;
				const ONE_READ *r = pairedReads->Data[i];

				readVertices[i] = NULL;
				mapIt = kh_get(ReadToVertex, Graph->ReadMap, r->ReadIndex);
				if (mapIt != kh_end(Graph->ReadMap))
					readVertices[i] = kh_val(Graph->ReadMap, mapIt);
			}

			for (size_t i = 0; i < prCount; ++i) {
				if (readVertices[i] == NULL)
					continue;
				
				for (size_t j = 0; j < prCount; ++j) {
					if (readVertices[j] == NULL)
						continue;

					if (j == i)
						continue;

					for (size_t k = 0; k < pointer_array_size(readVertices[i]); ++k) {
						PVARIANT_GRAPH_VERTEX u = readVertices[i]->Data[k];
						
						if (!_vg_vertex_exists(Graph, u))
							continue;
						
						for (size_t l = 0; l < pointer_array_size(readVertices[j]); ++l) {
							PVARIANT_GRAPH_VERTEX v = readVertices[j]->Data[l];
							
							if (!_vg_vertex_exists(Graph, v))
								continue;

							if (u->Index == v->Index) {

							}

							PVARIANT_GRAPH_PAIRED_EDGE tmp = NULL;
							boolean exists = FALSE;

							for (size_t m = 0; m < u->PairedCount; ++m) {
								exists = (u->Paired[m].Target == v);
								if (exists) {
									u->Paired[m].Count++;
									break;
								}
							}

							if (!exists) {
								ret = utils_calloc_VARIANT_GRAPH_PAIRED_EDGE(u->PairedCount + 1, &tmp);
								if (ret == ERR_SUCCESS) {
									memcpy(tmp, u->Paired, u->PairedCount*sizeof(VARIANT_GRAPH_PAIRED_EDGE));
									tmp[u->PairedCount].Target = v;
									tmp[u->PairedCount].Count = 1;
									++u->PairedCount;
									if (u->Paired != NULL)
										utils_free(u->Paired);

									u->Paired = tmp;
								}
							}

							if (ret != ERR_SUCCESS)
								break;
						}

						if (ret != ERR_SUCCESS)
							break;
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			utils_free(readVertices);
		}

		if (ret == ERR_SUCCESS)
			ret = paired_reads_next(pairedIter, &pairedIter, &pairedReads);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = _compute_components(Graph);

	return ret;
}


void vg_graph_finalize(PVARIANT_GRAPH Graph)
{
	PVARIANT_GRAPH_VERTEX v = Graph->Vertices.ByType.Alternative;

	
	for (size_t i = 0; i < Graph->VerticesArraySize; ++i) {
		if (!_vg_vertex_exists(Graph, v))
			v->Variant->Valid = FALSE;

		if (strlen(v->Variant->Ref) >= 40)
			v->Variant->Valid = FALSE;

		++v;
	}
	
	for (size_t i = 0; i < gen_array_size(&Graph->ComponentIndices) - 1; ++i) {
		const size_t startIndex = Graph->ComponentIndices.Data[i];
		const size_t endIndex = Graph->ComponentIndices.Data[i + 1];
		const uint64_t phasePos = Graph->Components.Data[startIndex]->Variant->Pos;

		for (size_t j = startIndex; j < endIndex; ++j) {
			PVARIANT_GRAPH_VERTEX v = Graph->Components.Data[j];
			PVARIANT_GRAPH_VERTEX var = _vg_vertex_get_variant(Graph, v);

			if (_vg_vertex_exists(Graph, v) && 
				v->Color != vgvcNone &&
				v->Type == vgvtAlternative &&
				v->Variant->PhasedPos == 0) {
				v->Variant->PhasedPos = phasePos;
				switch (v->Color) {
					case vgvcSequence1:
						v->Variant->PhaseType = vcptTwoOne;
						break;
					case vgvcSequence2:
						v->Variant->PhaseType = vcptOneTwo;
						break;
					case vgvcBoth:
						v->Variant->PhaseType = vcptBothAlt;
						break;
					default:
						assert(FALSE);
						break;
				}
			}
		}
	}

	return;
}
