
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "err.h"
#include "utils.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"
#include "found-sequence.h"
#include "kmer-graph.h"


/************************************************************************/
/*                      HELPER FUNCTIONS                                */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_GRAPH)
static UTILS_TYPED_MALLOC_FUNCTION(KMER_VERTEX)
static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE)

/************************************************************************/
/*                        VERTEX BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _vertex_create(const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_VERTEX(&tmp);
	if (ret == ERR_SUCCESS) {
		ret = kmer_copy(&tmp->KMer, KMer);
		if (ret == ERR_SUCCESS) {
			tmp->DegreeIn = 0;
			tmp->degreeOut = 0;
			tmp->Type = Type;
			dym_array_create(&tmp->Successors, 140);
			dym_array_create(&tmp->Predecessors, 140);
			tmp->Finished = FALSE;
			tmp->LastRSOrder = 0;
			*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmp);
	}

	return ret;
}


static void _vertex_destroy(PKMER_VERTEX Vertex)
{
	dym_array_destroy(&Vertex->Predecessors);
	dym_array_destroy(&Vertex->Successors);
	kmer_free(Vertex->KMer);
	utils_free(Vertex);

	return;
}


static ERR_VALUE _vertex_copy(const KMER_VERTEX *Vertex, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(Vertex->KMer, Vertex->Type, &tmp);
	if (ret == ERR_SUCCESS) {
		tmp->DegreeIn = Vertex->DegreeIn;
		tmp->degreeOut = Vertex->degreeOut;
		tmp->Order = Vertex->Order;
		dym_array_create(&tmp->Successors, 140);
		ret = dym_array_copy(&tmp->Successors, &Vertex->Successors);
		if (ret == ERR_SUCCESS) {
			ret = dym_array_copy(&tmp->Predecessors, &Vertex->Predecessors);
			if (ret == ERR_SUCCESS)
				*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(tmp);
	}

	return ret;
}

/************************************************************************/
/*                        EDGE BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _edge_create(PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, const long Weight, const uint32_t Length, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Source = Source;
		tmp->Dest = Dest;
		tmp->Type = Type;
		tmp->Length = Length;
		tmp->Weight = Weight;
		tmp->PassCount = 0;
		tmp->MaxPassCount = 0;
		tmp->Order = 0;
		tmp->Seq = NULL;
		tmp->SeqLen = 0;
		tmp->Seq2 = NULL;
		tmp->Seq2Len = 0;
		read_info_init(&tmp->ReadInfo);
		*Edge = tmp;
	}

	return ret;
}


static void _edge_destroy(PKMER_EDGE Edge)
{
	read_info_finit(&Edge->ReadInfo);
	if (Edge->Seq2Len > 0)
		utils_free(Edge->Seq2);

	if (Edge->SeqLen > 0)
		utils_free(Edge->Seq);

	utils_free(Edge);

	return;
}


static ERR_VALUE _edge_copy(const KMER_EDGE *Edge, PKMER_EDGE *Result)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _edge_create(Edge->Source, Edge->Dest, Edge->Type, Edge->Weight, Edge->Length, &tmp);
	if (ret == ERR_SUCCESS) {
		if (Edge->SeqLen > 0)
			ret = utils_copy_string(Edge->Seq, &tmp->Seq);

		if (ret == ERR_SUCCESS) {
			tmp->SeqLen = Edge->SeqLen;
			if (Edge->Seq2Len > 0)
				ret = utils_copy_string(Edge->Seq2, &tmp->Seq2);
			
			if (ret == ERR_SUCCESS) {
				tmp->Seq2Len = Edge->Seq2Len;
				ret = read_info_copy(&tmp->ReadInfo, &Edge->ReadInfo);
				if (ret == ERR_SUCCESS) {
					tmp->PassCount = Edge->PassCount;
					tmp->MaxPassCount = Edge->MaxPassCount;
					tmp->Order = Edge->Order;
					*Result = tmp;
				}

				if (ret != ERR_SUCCESS) {
					if (tmp->Seq2Len > 0)
						utils_free(tmp->Seq2);
				}
			}

			if (ret != ERR_SUCCESS) {
				if (tmp->SeqLen > 0)
					utils_free(tmp->Seq);
			}
		}

		if (ret != ERR_SUCCESS)
			_edge_destroy(tmp);
	}

	return ret;
}


/************************************************************************/
/*                      VERTEX TABLE CALBLACKS                          */
/************************************************************************/


static void _vertex_table_on_insert(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order)
{
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;

	v->Order = Order;

	return;
}

static void _vertex_table_on_delete(struct _KMER_TABLE *Table, void *ItemData)
{
	_vertex_destroy((PKMER_VERTEX)ItemData);

	return;
}


static ERR_VALUE _vertex_table_on_copy(struct _KMER_TABLE *Table, void *ItemData, void **Copy)
{
	return _vertex_copy((PKMER_VERTEX)ItemData, (PKMER_VERTEX *)Copy);
}


static void _vertex_table_on_print(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream)
{
	char *colors[] = {"yellow", "green", "blue", "red"};
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;


	fprintf(Stream, "\t");
	kmer_print(Stream, v->KMer);
	fprintf(Stream, "[label=\"");
	kmer_print(Stream, v->KMer);
	fprintf(Stream, "\\nIN=%u; OUT=%u; O=%u\",style=filled,color=%s]", v->DegreeIn, v->degreeOut, v->LastRSOrder, colors[v->Type]);
	fprintf(Stream, ";\n");

	return;
}


/************************************************************************/
/*                      EDGE TABLE CALBLACKS                          */
/************************************************************************/


static void _edge_table_on_insert(struct _KMER_EDGE_TABLE *Table, void *ItemData, const uint32_t Order)
{
	PKMER_EDGE e = (PKMER_EDGE)ItemData;

	e->Order = Order;

	return;
}

static void _edge_table_on_delete(struct _KMER_EDGE_TABLE *Table, void *ItemData)
{
	_edge_destroy((PKMER_EDGE)ItemData);

	return;
}


static ERR_VALUE _edge_table_on_copy(struct _KMER_EDGE_TABLE *Table, void *ItemData, void **Copy)
{
	return _edge_copy((PKMER_EDGE)ItemData, (PKMER_EDGE *)Copy);
}


static void _edge_table_on_print(struct _KMER_EDGE_TABLE *Table, void *ItemData, FILE *Stream)
{
	PKMER_EDGE e = (PKMER_EDGE)ItemData;

	fprintf(Stream, "\t");
	kmer_print(Stream, e->Source->KMer);
	fprintf(Stream, " -> ");
	kmer_print(Stream, e->Dest->KMer);
	fprintf(Stream, " [color=");
	switch (e->Type) {
		case kmetReference:
			fprintf(Stream, "green");
			break;
		case kmetRead:
			fprintf(Stream, "red");
			break;
		case kmetVariant:
			fprintf(Stream, "blue");
			break;
	}

	fprintf(Stream, ",label=\"W: %li; M: %u; L %u\"", e->Weight, e->MaxPassCount, e->Length);
	fprintf(Stream, "];\n");

	return;
}


static boolean _kmer_vertex_no_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = TRUE;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		ret = (kmer_vertex_get_succ_edge(Vertex, i))->Type != kmetRead;
		if (!ret)
			break;
	}

	if (ret) {
		for (size_t i = 0; i < Vertex->DegreeIn; ++i) {
			ret = (kmer_vertex_get_pred_edge(Vertex, i))->Type != kmetRead;
			if (!ret)
				break;
		}
	}

	return ret;
}


/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	KMER_TABLE_CALLBACKS vCallbacks;
	KMER_EDGE_TABLE_CALLBACKS eCallbacks;

	ret = utils_malloc_KMER_GRAPH(&tmpGraph);
	if (ret == ERR_SUCCESS) {
		tmpGraph->NumberOfEdges = 0;
		tmpGraph->NumberOfVertices = 0;
		tmpGraph->KMerSize = KMerSize;
		tmpGraph->NumberOfBackwardEdges = 0;
		tmpGraph->StartingVertex = NULL;
		tmpGraph->EndingVertex = NULL;
		vCallbacks.OnCopy = _vertex_table_on_copy;
		vCallbacks.OnDelete = _vertex_table_on_delete;
		vCallbacks.OnInsert = _vertex_table_on_insert;
		vCallbacks.OnPrint = _vertex_table_on_print;
		ret = kmer_table_create(KMerSize, 2, 37, &vCallbacks, &tmpGraph->VertexTable);
		if (ret == ERR_SUCCESS) {
			eCallbacks.OnCopy = _edge_table_on_copy;
			eCallbacks.OnDelete = _edge_table_on_delete;
			eCallbacks.OnInsert = _edge_table_on_insert;
			eCallbacks.OnPrint = _edge_table_on_print;
			ret = kmer_edge_table_create(KMerSize, 2, 37, &eCallbacks, &tmpGraph->EdgeTable);
			if (ret == ERR_SUCCESS)
				*Graph = tmpGraph;

			if (ret != ERR_SUCCESS)
				kmer_table_destroy(tmpGraph->VertexTable);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpGraph);
	}

	return ret;
}


void kmer_graph_destroy(PKMER_GRAPH Graph)
{
	kmer_edge_table_destroy(Graph->EdgeTable);
	kmer_table_destroy(Graph->VertexTable);
	utils_free(Graph);

	return;
}


void kmer_graph_print(FILE *Stream, const KMER_GRAPH *Graph)
{
	fprintf(Stream, "digraph G {\n");
	fprintf(Stream, "\t/* number of vertices: %u */\n", Graph->NumberOfVertices);
	fprintf(Stream, "\t/* number of edges: %u */\n", Graph->NumberOfEdges);
	fprintf(Stream, "\t/* number of backward edges: %u */\n", Graph->NumberOfBackwardEdges);
	kmer_table_print(Stream, Graph->VertexTable);
	kmer_edge_table_print(Stream, Graph->EdgeTable);
	fprintf(Stream, "}\n");

	return;
}


ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph)
{
	PKMER_VERTEX v = NULL;
	PKMER_TABLE_ENTRY e = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	PKMER_VERTEX x = NULL;
	PKMER_VERTEX y = NULL;
	PKMER_EDGE dummy = NULL;
	PKMER_EDGE sourceEdge = NULL;
	PKMER_EDGE destEdge = NULL;
	long weight = 0;
	long sourceWeight = 0;
	long destWeight = 0;
	uint32_t length = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		boolean last = FALSE;
		
		e = iter;
		last = (kmer_table_next(Graph->VertexTable, e, &iter) == ERR_NO_MORE_ENTRIES);
		v = (PKMER_VERTEX)e->Data;
		if (v->DegreeIn == 1 && v->degreeOut == 1)
			kmer_graph_delete_vertex(Graph, v);

		if (last)
			break;
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	return ret;
}


void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold)
{
	PKMER_EDGE e = NULL;
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	boolean last = FALSE;
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX v = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (err == ERR_SUCCESS) {
		boolean last = FALSE;

		e = (PKMER_EDGE)(iter->Data);
		last = (kmer_edge_table_next(Graph->EdgeTable, iter, &iter) == ERR_NO_MORE_ENTRIES);
		if (e->Weight < Threshold)
			kmer_graph_delete_edge(Graph, e);

		if (last)
			break;
	}

	return;
}


void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer)
{
	assert(Graph->StartingVertex == NULL);
	Graph->StartingVertex = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	assert(Graph->StartingVertex != NULL);

	return;
}


void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer)
{
	assert(Graph->EndingVertex == NULL);
	Graph->EndingVertex = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	assert(Graph->EndingVertex != NULL);

	return;
}


void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph, size_t *DeletedThings)
{
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret != ERR_NO_MORE_ENTRIES) {
		boolean deleted = FALSE;
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;

		if (v->Type != kmvtRefSeqStart && v->Type != kmvtRefSeqEnd) {
			if (v->DegreeIn == 0) {
				deleted = TRUE;
				while (v->degreeOut > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_succ_edge(v, 0));
			} else if (v->degreeOut == 0) {
				deleted = TRUE;
				while (v->DegreeIn > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_pred_edge(v, 0));
			}

			if (deleted) {
				kmer_graph_delete_vertex(Graph, v);
				*DeletedThings += 1;
			}
		}

		ret = (deleted) ?
			kmer_table_first(Graph->VertexTable, &iter) :
			kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}


static ERR_VALUE _get_bubble_ends(PKMER_VERTEX FirstVertex, PDYM_ARRAY EndArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (FirstVertex->Type == kmvtRead) {
		ret = ERR_SUCCESS;
		DYM_ARRAY vertexStack;
		DYM_ARRAY edgeStack;
		PKMER_VERTEX v = FirstVertex;
		size_t edgeIndex = 0;

		ret = ERR_SUCCESS;
		dym_array_create(&vertexStack, 140);
		dym_array_create(&edgeStack, 140);
		do {
			if (edgeIndex == v->degreeOut) {
				v->Finished = FALSE;
				v = dym_array_remove_back(&vertexStack);
				edgeIndex = (size_t)dym_array_remove_back(&edgeStack);
				continue;
			}

			for (size_t i = edgeIndex; i < v->degreeOut; ++i) {
				PKMER_VERTEX tmp = kmer_vertex_get_successor(v, i);
				
				if (!tmp->Finished) {
					ret = dym_array_push_back(&vertexStack, v);
					if (ret != ERR_SUCCESS)
						break;

					ret = dym_array_push_back(&edgeStack, (void *)(i + 1));
					if (ret != ERR_SUCCESS)
						break;

					v->Finished = TRUE;
					v = tmp;
					break;
				}
			}

			if (ret == ERR_SUCCESS) {
				if (v->Type == kmvtRefSeqMiddle) {
					ret = dym_array_push_back(EndArray, v);
					if (ret != ERR_SUCCESS)
						break;

					v->Finished = FALSE;
					v = dym_array_remove_back(&vertexStack);
					edgeIndex = (size_t)dym_array_remove_back(&edgeStack);
				}
			}
		} while (ret == ERR_SUCCESS && (v != FirstVertex || edgeIndex < v->degreeOut));

		FirstVertex->Finished = FALSE;
		dym_array_destroy(&edgeStack);
		dym_array_destroy(&vertexStack);
	} else ret = dym_array_push_back(EndArray, FirstVertex);

	return ret;
}


static PKMER_EDGE _get_refseq_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_succ_edge(Vertex, i);

		if (tmp->Type == kmetReference) {
			ret = tmp;
			break;
		}
	}

	return ret;
}


static boolean _has_outgoing_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = FALSE;
	const KMER_EDGE *e = NULL;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		PKMER_EDGE e = kmer_vertex_get_succ_edge(Vertex, i);

		ret = (e->Type == kmetRead);
		if (ret)
			break;
	}

	return ret;
}


ERR_VALUE kmer_graph_resolve_bubbles(PKMER_GRAPH Graph, const uint32_t Threshold)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX bubbleStart = Graph->StartingVertex;
	size_t readEdgeIndex = (size_t)-1;

	ret = ERR_SUCCESS;
	do {		
		PKMER_EDGE rsEdge = _get_refseq_edge(bubbleStart);
		
		for (size_t i = 0; i < bubbleStart->degreeOut; ++i) {
			PKMER_EDGE tmp = kmer_vertex_get_succ_edge(bubbleStart, i);

			if (tmp->Type == kmetRead) {
				DYM_ARRAY bubbleEnds;
				
				readEdgeIndex = i;
				dym_array_create(&bubbleEnds, 140);
				ret = _get_bubble_ends(tmp->Dest, &bubbleEnds);
				if (ret == ERR_SUCCESS) {
					PKMER_VERTEX v = rsEdge->Dest;
					boolean noReadSupport = rsEdge->Weight == Threshold;

					while (!_has_outgoing_read_edges(v) && dym_array_size(&bubbleEnds) > 0 && v != Graph->EndingVertex) {
						size_t index = 0;

						if (dym_array_find(&bubbleEnds, v, &index)) {
							dym_array_remove_by_item_fast(&bubbleEnds, v);
							if (noReadSupport)
								kmer_graph_delete_edge(Graph, rsEdge);

							rsEdge = _get_refseq_edge(v);
							noReadSupport = rsEdge->Weight == Threshold;
							v = rsEdge->Dest;
						} else {
							rsEdge = _get_refseq_edge(v);
							noReadSupport &= (rsEdge->Weight == Threshold);
							v = rsEdge->Dest;
						}
					}
				}

				readEdgeIndex = (size_t)-1;
				dym_array_destroy(&bubbleEnds);
				break;
			}
		}

		if (ret == ERR_SUCCESS) {
			if (readEdgeIndex == (size_t)-1)
				bubbleStart = rsEdge->Dest;
		}
	} while (ret == ERR_SUCCESS && bubbleStart != Graph->EndingVertex);

	size_t dummy = 0;
	kmer_graph_delete_trailing_things(Graph, &dummy);

	return ret;
}


ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(KMer, Type, &v);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_table_insert(Graph->VertexTable, KMer, v);
			if (ret == ERR_TABLE_FULL) {
				ret = kmer_table_extend(Graph->VertexTable);
				if (ret == ERR_SUCCESS)
					ret = ERR_TABLE_FULL;
			}
		} while (ret == ERR_TABLE_FULL);

		if (ret == ERR_SUCCESS) {
			Graph->NumberOfVertices++;
			if (Vertex != NULL)
				*Vertex = v;
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(v);
	}

	return ret;
}


ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type)
{
	return kmer_graph_add_vertex_ex(Graph, KMer, Type, NULL);
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX v = NULL;
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	u = kmer_table_get(Graph->VertexTable, Source);
	v = kmer_table_get(Graph->VertexTable, Dest);
	if (u != NULL && v != NULL) {
		ret = _edge_create(u, v, Type, weight, Length, &edge);
		if (ret == ERR_SUCCESS) {
			do {
				ret = kmer_edge_table_insert(Graph->EdgeTable, Source, Dest, edge);
				if (ret == ERR_SUCCESS) {
					ret = dym_array_prepare_for_insert(&u->Successors, 1);
					if (ret == ERR_SUCCESS) {
						ret = dym_array_prepare_for_insert(&v->Predecessors, 1);
						if (ret == ERR_SUCCESS) {
							u->degreeOut++;
							dym_array_push_back_no_alloc(&u->Successors, edge);
							v->DegreeIn++;
							dym_array_push_back_no_alloc(&v->Predecessors, edge);
							Graph->NumberOfEdges++;
							*Edge = edge;
						}
					}
				}

				if (ret == ERR_TABLE_FULL) {
					ret = kmer_edge_table_extend(Graph->EdgeTable);
					if (ret == ERR_SUCCESS)
						ret = ERR_TABLE_FULL;
				}
			} while (ret == ERR_TABLE_FULL);

			if (ret != ERR_SUCCESS)
				_edge_destroy(edge);

			if (ret == ERR_ALREADY_EXISTS)
				*Edge = (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source, Dest);
		}
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight)
{
	PKMER_EDGE dummy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_add_edge_ex(Graph, Source, Dest, weight, 1, kmetReference, &dummy);

	return ret;
}


PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest)
{
	return (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source, Dest);
}


PKMER_VERTEX kmer_graph_get_vertex(const struct _KMER_GRAPH *Graph, const struct _KMER *KMer)
{
	return (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
}


ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PGEN_ARRAY_PKMER_VERTEX VertexArray)
{
	DYM_ARRAY tmp;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	dym_array_create(&tmp, 140);
	ret = kmer_table_get_multiple(Graph->VertexTable, KMer, &tmp);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < dym_array_size(&tmp); ++i) {
			ret = dym_array_push_back_PKMER_VERTEX(VertexArray, (PKMER_VERTEX)dym_array_get(&tmp, i));
			if (ret != ERR_SUCCESS)
				break;
		}
	}

	dym_array_destroy(&tmp);

	return ret;
}


ERR_VALUE kmer_graph_delete_vertex(PKMER_GRAPH Graph, PKMER_VERTEX Vertex)
{
	PKMER_EDGE inEdge = NULL;
	PKMER_EDGE outEdge = NULL;
	PKMER_VERTEX predVertex = NULL;
	PKMER_VERTEX succVertex = NULL;
	PKMER_EDGE newEdge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Vertex->DegreeIn <= 1 && Vertex->degreeOut <= 1) {
		if (Vertex->degreeOut == 1) {
			outEdge = kmer_vertex_get_succ_edge(Vertex, 0);
			succVertex = outEdge->Dest;
		}

		if (Vertex->DegreeIn == 1) {
			inEdge = kmer_vertex_get_pred_edge(Vertex, 0);
			predVertex = inEdge->Source;
		}

		if (predVertex == NULL || (predVertex != succVertex)) {
			if (inEdge != NULL && outEdge != NULL) {
				ret = kmer_graph_merge_edges(Graph, inEdge, outEdge);
			} else if (inEdge != NULL) {
				kmer_graph_delete_edge(Graph, inEdge);
				ret = ERR_SUCCESS;
			} else if (outEdge != NULL) {
				kmer_graph_delete_edge(Graph, outEdge);
				ret = ERR_SUCCESS;
			} else ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				kmer_table_delete(Graph->VertexTable, Vertex->KMer);
				--Graph->NumberOfVertices;
			}
		} else ret = ERR_PRED_IS_SUCC;
	} else ret = ERR_TOO_MANY_EDGES;

	return ret;
}


void kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge)
{
	ERR_VALUE err = ERR_INTERNAL_ERROR;
	PKMER_VERTEX source = Edge->Source;
	PKMER_VERTEX dest = Edge->Dest;

	err = kmer_edge_table_delete(Graph->EdgeTable, source->KMer, dest->KMer);
	if (err == ERR_SUCCESS) {
		dym_array_remove_by_item_fast(&source->Successors, Edge);
		dym_array_remove_by_item_fast(&dest->Predecessors, Edge);
		--source->degreeOut;
		--dest->DegreeIn;
		--Graph->NumberOfEdges;
	} else __debugbreak();

	return;
}


ERR_VALUE kmer_graph_merge_edges(PKMER_GRAPH Graph, PKMER_EDGE Source, PKMER_EDGE Dest)
{
	PKMER_VERTEX v = NULL;
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX w = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	u = Source->Source;
	v = Source->Dest;
	w = Dest->Dest;
	if (v == Dest->Source && u != w) {
		if (kmer_graph_get_edge(Graph, u->KMer, w->KMer) == NULL) {
			long weight = max(Source->Weight, Dest->Weight);
			uint32_t Length = Source->Length + Dest->Length;
			uint32_t mpc = Source->MaxPassCount;
			EKMerEdgeType type = (Source->Type == kmetReference && Dest->Type == kmetReference) ? kmetReference : kmetRead;
			PKMER_EDGE newEdge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, u->KMer, w->KMer, weight, Length, type, &newEdge);
			if (ret == ERR_SUCCESS) {
				newEdge->MaxPassCount = mpc;
				newEdge->SeqLen = Source->SeqLen + 1 + Dest->SeqLen;
				ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
				if (ret == ERR_SUCCESS) {
					memcpy(newEdge->Seq, Source->Seq, Source->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
					memcpy(newEdge->Seq + Source->SeqLen + 1, Dest->Seq, Dest->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen + 1 + Dest->SeqLen] = '\0';
					kmer_graph_delete_edge(Graph, Source);
					kmer_graph_delete_edge(Graph, Dest);					
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, newEdge);
			}
		} else ret = ERR_TRIANGLE;
	} else ret = (u == w) ? ERR_PRED_IS_SUCC : ERR_NOT_ADJACENT;

	return ret;
}


ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, PGEN_ARRAY_PFOUND_SEQUENCE SeqArray)
{
	size_t l = 0;
	char s[16384];
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX v = Graph->StartingVertex;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	size_t edgeIndex = 0;
	PKMER_EDGE edge = NULL;
	GEN_ARRAY_size_t edgeIndices;
	GEN_ARRAY_PKMER_EDGE edges;
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT variantStack;

	dym_array_init_FOUND_SEQUENCE_VARIANT(&variantStack, 140);
	dym_array_init_size_t(&edgeIndices, 140);
	dym_array_init_PKMER_EDGE(&edges, 140);
	memset(s, 0, sizeof(s));
	memcpy(s, v->KMer->Bases + 1, kmerSize - 1);
	l += (kmerSize - 1);
	do {
		if (edgeIndex == v->degreeOut) {
			v->Finished = FALSE;
			if (v == Graph->StartingVertex)
				break;

			edge = *dym_array_pop_back_PKMER_EDGE(&edges);
			edgeIndex = *dym_array_pop_back_size_t(&edgeIndices);
			v = edge->Source;
			if (edge->Type == kmetVariant) {
				l -= 2;
				s[l] = '\0';
				dym_array_pop_back_FOUND_SEQUENCE_VARIANT(&variantStack);
			} else {
				l -= (edge->SeqLen + 1);
				memset(s + l, 0, (edge->SeqLen + 1)*sizeof(char));
			}
		} else {
			edge = kmer_vertex_get_succ_edge(v, edgeIndex);
			if (edge->Type == kmetVariant) {
				FOUND_SEQUENCE_VARIANT variant;
				
				s[l] = '?';
				++l;
				variant.Seq1 = edge->Seq;
				variant.Seq1Len = edge->SeqLen;
				variant.Seq2 = edge->Seq2;
				variant.Seq2Len = edge->Seq2Len;
				ret = dym_array_push_back_FOUND_SEQUENCE_VARIANT(&variantStack, variant);
				if (ret != ERR_SUCCESS)
					break;
			} else {
				memcpy(s + l, edge->Seq, edge->SeqLen*sizeof(char));
				l += edge->SeqLen;
			}

			if (edge->Dest == Graph->EndingVertex) {
				PFOUND_SEQUENCE fs = NULL;
				const size_t vCount = gen_array_size(&variantStack);

				ret = found_sequence_alloc(s, l, vCount, &fs);
				if (ret == ERR_SUCCESS) {
					for (size_t i = 0; i < vCount; ++i) {
						ret = found_sequence_set_variant(fs, i, dym_array_item_FOUND_SEQUENCE_VARIANT(&variantStack, i));
						if (ret != ERR_SUCCESS)
							break;
					}

					ret = dym_array_push_back_PFOUND_SEQUENCE(SeqArray, fs);
					if (ret == ERR_SUCCESS) {
						l -= edge->SeqLen;
						memset(s + l, 0, edge->SeqLen*sizeof(char));
						++edgeIndex;
						if (gen_array_size(SeqArray) >= 300000)
							break;
					}

					if (ret != ERR_SUCCESS)
						found_sequence_free(fs);
				}
			} else {
				v = edge->Dest;
				if (!v->Finished) {
					edge->Source->Finished = TRUE;
					s[l] = kmer_get_base(v->KMer, kmerSize - 1);
					++l;
					ret = dym_array_push_back_size_t(&edgeIndices, (edgeIndex + 1));
					if (ret != ERR_SUCCESS)
						break;

					ret = dym_array_push_back_PKMER_EDGE(&edges, edge);
					if (ret != ERR_SUCCESS)
						break;

					edgeIndex = 0;
				} else {
					v = edge->Source;
					if (edge->Type == kmetVariant) {
						--l;
						s[l] = '\0';
						dym_array_pop_back_FOUND_SEQUENCE_VARIANT(&variantStack);
					} else {
						l -= edge->SeqLen;
						memset(s + l, 0, edge->SeqLen*sizeof(char));
					}

					++edgeIndex;
				}
			}
		}
	} while (ret == ERR_SUCCESS && v != Graph->StartingVertex);

	dym_array_finit_PKMER_EDGE(&edges);
	dym_array_finit_size_t(&edgeIndices);
	dym_array_finit_FOUND_SEQUENCE_VARIANT(&variantStack);
	if (ret != ERR_SUCCESS) {
		for (size_t i = 0; i < gen_array_size(SeqArray); ++i)
			found_sequence_free(*dym_array_item_PFOUND_SEQUENCE(SeqArray, i));

		dym_array_clear_PFOUND_SEQUENCE(SeqArray);
	}

	return ret;
}


ERR_VALUE kmer_vertex_get_certain_edges(const KMER_VERTEX *Vertex, const EKMerEdgeType EdgeType, const boolean Incomming, PGEN_ARRAY_PKMER_EDGE Array)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const PKMER_EDGE *edgeArray = (PKMER_EDGE *)dym_array_data((Incomming) ? (&Vertex->Predecessors) : (&Vertex->Successors));
	const size_t count = dym_array_size((Incomming) ? (&Vertex->Predecessors) : (&Vertex->Successors));

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < count; ++i) {
		const KMER_EDGE *e = *edgeArray;
	
		if (e->Type == EdgeType)
			ret = dym_array_push_back_PKMER_EDGE(Array, e);
		
		if (ret != ERR_SUCCESS)
			break;

		++edgeArray;
	}

	return ret;
}


static ERR_VALUE _capture_refseq(const KMER_EDGE *Start, const KMER_EDGE *End, char **Seq, size_t *SeqLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char seq[16384];
	size_t len = 0;
	const size_t kmerSize = kmer_get_size(Start->Dest->KMer);

	memcpy(seq + len, Start->Seq, Start->SeqLen*sizeof(char));
	len += Start->SeqLen;
	if (Start->Dest == End->Source) {
		seq[len] = kmer_get_base(Start->Dest->KMer, kmerSize - 1);
		++len;
		memcpy(seq + len, End->Seq, End->SeqLen*sizeof(char));
		len += End->SeqLen;
	} else if (Start != End) {
		do {
			seq[len] = kmer_get_base(Start->Dest->KMer, kmerSize - 1);
			++len;
			Start = _get_refseq_edge(Start->Dest);
			memcpy(seq + len, Start->Seq, Start->SeqLen*sizeof(char));
			len += Start->SeqLen;
		} while (Start != End && End->Source != Start->Dest);

		if (End->Source == Start->Dest) {
			seq[len] = kmer_get_base(Start->Dest->KMer, kmerSize - 1);
			++len;
			memcpy(seq + len, End->Seq, End->SeqLen*sizeof(char));
			len += End->SeqLen;
		}
	}

	seq[len] = '\0';
	ret = utils_copy_string(seq, Seq);
	if (ret == ERR_SUCCESS)
		*SeqLen = len;

	return ret;
}


static ERR_VALUE _remove_read_info_from_edges(PKMER_GRAPH Graph, PKMER_EDGE Edge1, PKMER_EDGE Edge2, PGEN_ARRAY_PKMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY Intersection)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_READ_INFO_ENTRY diff;

	dym_array_init_READ_INFO_ENTRY(&diff, 140);
	ret = read_info_diff(&Edge1->ReadInfo, Intersection, &diff);
	if (ret == ERR_SUCCESS) {
		Edge1->Weight -= (read_info_get_count(&Edge1->ReadInfo) - gen_array_size(&diff));
		if (Edge1->Weight < 0)
			Edge1->Weight = 0;

		ret = read_info_assign(&Edge1->ReadInfo, &diff);
	}

	if (ret == ERR_SUCCESS) {
		dym_array_clear_READ_INFO_ENTRY(&diff);
		ret = read_info_diff(&Edge2->ReadInfo, Intersection, &diff);
		if (ret == ERR_SUCCESS) {
			Edge2->Weight -= (read_info_get_count(&Edge2->ReadInfo) - gen_array_size(&diff));
			if (Edge2->Weight < 0)
				Edge2->Weight = 0;

			ret = read_info_assign(&Edge2->ReadInfo, &diff);
		}

		if (ret == ERR_SUCCESS) {
			for (size_t k = 0; k < gen_array_size(Array); ++k) {
				PKMER_EDGE tmp = *dym_array_item_PKMER_EDGE(Array, k);

				dym_array_clear_READ_INFO_ENTRY(&diff);
				ret = read_info_diff(&tmp->ReadInfo, Intersection, &diff);
				if (ret == ERR_SUCCESS) {
					tmp->Weight -= (read_info_get_count(&tmp->ReadInfo) - gen_array_size(&diff));
					if (tmp->Weight < 0)
						tmp->Weight = 0;

					ret = read_info_assign(&tmp->ReadInfo, &diff);
				}

				if (ret != ERR_SUCCESS)
					break;
			}
		}
	}

	dym_array_finit_READ_INFO_ENTRY(&diff);

	return ret;
}


ERR_VALUE kmer_graph_connect_reads_by_refseq(PKMER_GRAPH Graph, const size_t Threshold)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX passStart = Graph->StartingVertex;
	PKMER_VERTEX passEnd = NULL;
	GEN_ARRAY_PKMER_EDGE incommingReads;
	GEN_ARRAY_PKMER_EDGE outgoingReads;
	GEN_ARRAY_PKMER_EDGE rsPathEdges;
	GEN_ARRAY_READ_INFO_ENTRY intersection;
	const KMER_EDGE *e = NULL;
	const KMER_EDGE *rsLastEdge = NULL;
	const KMER_EDGE *rsNextEdge = NULL;
	const KMER_EDGE *nextRsLastEdge = NULL;

	ret = ERR_SUCCESS;
	dym_array_init_PKMER_EDGE(&incommingReads, 140);
	dym_array_init_PKMER_EDGE(&outgoingReads, 140);
	dym_array_init_PKMER_EDGE(&rsPathEdges, 140);
	dym_array_init_READ_INFO_ENTRY(&intersection, 140);
	passStart = Graph->StartingVertex;
	while (ret == ERR_SUCCESS && passEnd != Graph->EndingVertex && passStart != Graph->EndingVertex) {		
		while (ret == ERR_SUCCESS && passStart != Graph->EndingVertex) {
			ret = kmer_vertex_get_certain_edges(passStart, kmetRead, TRUE, &incommingReads);
			if (ret == ERR_SUCCESS && gen_array_size(&incommingReads) > 0)
				break;

			rsLastEdge = _get_refseq_edge(passStart);
			passStart = rsLastEdge->Dest;
		}

		if (passStart == Graph->EndingVertex)
			break;

		nextRsLastEdge = rsLastEdge;
		if (ret == ERR_SUCCESS && gen_array_size(&incommingReads) > 0) {
			GEN_ARRAY_PKMER_EDGE tmp;
			size_t readDistance = 0;

			passEnd = passStart;
			dym_array_init_PKMER_EDGE(&tmp, 140);
			while (ret == ERR_SUCCESS && passEnd != Graph->EndingVertex) {
				ret = kmer_vertex_get_certain_edges(passEnd, kmetRead, FALSE, &outgoingReads);
				if (ret == ERR_SUCCESS && gen_array_size(&outgoingReads) > 0)
					break;

				e = _get_refseq_edge(passEnd);
				nextRsLastEdge = e;
				readDistance += e->Length;
				passEnd = e->Dest;
				ret = kmer_vertex_get_certain_edges(passEnd, kmetRead, TRUE, &tmp);
				if (ret == ERR_SUCCESS && gen_array_size(&tmp) > 0)
					break;

				ret = dym_array_push_back_PKMER_EDGE(&rsPathEdges, e);
			}

			dym_array_finit_PKMER_EDGE(&tmp);
			if (ret == ERR_SUCCESS && gen_array_size(&outgoingReads) > 0) {
				const size_t inCount = gen_array_size(&incommingReads);
				const size_t outCount = gen_array_size(&outgoingReads);
				GEN_ARRAY_PKMER_EDGE edgesToDelete;

				dym_array_init_PKMER_EDGE(&edgesToDelete, 140);
				ret = dym_array_reserve_PKMER_EDGE(&edgesToDelete, inCount*outCount);
				if (ret == ERR_SUCCESS) {
					PKMER_EDGE eIn = NULL;
					PKMER_EDGE eOut = NULL;

					rsNextEdge = _get_refseq_edge(passEnd);
					for (size_t i = 0; i < inCount; ++i) {
						boolean edgeCreated = FALSE;

						eIn = *dym_array_item_PKMER_EDGE(&incommingReads, i);
						for (size_t j = 0; j < outCount; ++j) {
							edgeCreated = FALSE;
							eOut = *dym_array_item_PKMER_EDGE(&outgoingReads, j);
							if (eOut != NULL) {
								if (eIn == eOut)
									continue;

								ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, TRUE,  eIn->Length + readDistance);
								if (ret == ERR_SUCCESS && gen_array_size(&intersection) > 0) {
									char *seq = NULL;
									size_t seqLen = 0;

									ret = _capture_refseq(eIn, eOut, &seq, &seqLen);
									if (ret == ERR_SUCCESS) {
										PKMER_EDGE newEdge = NULL;

										ret = kmer_graph_add_edge_ex(Graph, eIn->Source->KMer, eOut->Dest->KMer, gen_array_size(&intersection), seqLen + 1, kmetRead, &newEdge);
										if (ret == ERR_SUCCESS) {
											putchar('E');
											edgeCreated = TRUE;
											newEdge->Seq = seq;
											newEdge->SeqLen = seqLen;
											ret = read_info_assign(&newEdge->ReadInfo, &intersection);
											if (ret == ERR_SUCCESS)
												ret = _remove_read_info_from_edges(Graph, eIn, eOut, &rsPathEdges, &intersection);
										}

										if (ret != ERR_SUCCESS)
											utils_free(seq);

										if (ret == ERR_ALREADY_EXISTS)
											ret = ERR_SUCCESS;
									}

									dym_array_clear_READ_INFO_ENTRY(&intersection);
								}

								if (edgeCreated) {
									if (read_info_get_count(&rsNextEdge->ReadInfo) != 0 ||
										read_info_get_count(&rsLastEdge->ReadInfo) == 0) {
										putchar('D');
										*dym_array_item_PKMER_EDGE(&outgoingReads, j) = NULL;
										dym_array_push_back_no_alloc_PKMER_EDGE(&edgesToDelete, eOut);
										break;
									}
								}
							}

							if (ret != ERR_SUCCESS)
								break;
						}

						if (edgeCreated) {
							if (read_info_get_count(&rsNextEdge->ReadInfo) == 0 ||
								read_info_get_count(&rsLastEdge->ReadInfo) != 0) {
								putchar('D');
								dym_array_push_back_no_alloc_PKMER_EDGE(&edgesToDelete, eIn);
							}

							continue;
						}

						if (ret != ERR_SUCCESS)
							break;
					}

					for (size_t i = 0; i < gen_array_size(&edgesToDelete); ++i)
						kmer_graph_delete_edge(Graph, *dym_array_item_PKMER_EDGE(&edgesToDelete, i));
				}

				dym_array_finit_PKMER_EDGE(&edgesToDelete);
			}
		}

		dym_array_clear_PKMER_EDGE(&incommingReads);
		dym_array_clear_PKMER_EDGE(&outgoingReads);
		dym_array_clear_PKMER_EDGE(&rsPathEdges);
		if (passEnd == passStart && passStart != Graph->EndingVertex) {
			e = _get_refseq_edge(passEnd);
			rsLastEdge = e;
			passStart = e->Dest;
		} else {
			passStart = passEnd;
			rsLastEdge = nextRsLastEdge;
		}
	}

	dym_array_finit_READ_INFO_ENTRY(&intersection);
	dym_array_finit_PKMER_EDGE(&rsPathEdges);
	dym_array_finit_PKMER_EDGE(&outgoingReads);
	dym_array_finit_PKMER_EDGE(&incommingReads);
	printf("\n");

	return ret;
}


ERR_VALUE kmer_graph_connect_reads_by_reads(PKMER_GRAPH Graph, const size_t Threshold)
{
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const kmerSize = kmer_graph_get_kmer_size(Graph);

//	kmer_graph_print(stderr, Graph);
//	fflush(stderr);
	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;

		if (v->Type == kmvtRead && v->DegreeIn > 1) {
			PKMER_VERTEX w = v;
			char midSeq[8192];
			size_t midSeqLen = 0;
			boolean cycle = FALSE;
			size_t readDistance = 0;

			midSeq[midSeqLen] = kmer_get_base(w->KMer, kmerSize - 1);
			++midSeqLen;
			while (!cycle && w->Type == kmvtRead && w->degreeOut == 1 && (v == w || w->DegreeIn == 1)) {
				PKMER_EDGE e = kmer_vertex_get_succ_edge(w, 0);
				memcpy(midSeq + midSeqLen, e->Seq, e->SeqLen*sizeof(char));
				midSeqLen += e->SeqLen;
				readDistance += e->Length;
				w = e->Dest;
				midSeq[midSeqLen] = kmer_get_base(w->KMer, kmerSize - 1);
				++midSeqLen;
				cycle = (v == w);
			}

			if (!cycle && w->Type == kmvtRead && w->degreeOut == v->DegreeIn && (v == w || w->DegreeIn == 1)) {
				size_t *edgesCreatedPerOutEdge = NULL;

				ret = utils_calloc(w->degreeOut, sizeof(size_t), &edgesCreatedPerOutEdge);
				if (ret == ERR_SUCCESS) {
					GEN_ARRAY_PKMER_EDGE edgesToDelete;
					
					memset(edgesCreatedPerOutEdge, 0, sizeof(size_t)*w->degreeOut);
					dym_array_init_PKMER_EDGE(&edgesToDelete, 140);
					ret = dym_array_reserve_PKMER_EDGE(&edgesToDelete, v->DegreeIn + w->degreeOut);
					if (ret == ERR_SUCCESS) {
						GEN_ARRAY_READ_INFO_ENTRY intersection;

						dym_array_init_READ_INFO_ENTRY(&intersection, 140);
						for (size_t i = 0; i < v->DegreeIn; ++i) {
							size_t numberOfCreatedEdges = 0;
							PKMER_EDGE eIn = kmer_vertex_get_pred_edge(v, i);

							for (size_t j = 0; j < w->degreeOut; ++j) {
								PKMER_EDGE eOut = kmer_vertex_get_succ_edge(w, j);

								ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, TRUE, eIn->Length + readDistance);
								if (ret == ERR_SUCCESS) {
									if (gen_array_size(&intersection) > 0) {
										PKMER_EDGE newEdge = NULL;

										ret = kmer_graph_add_edge_ex(Graph, eIn->Source->KMer, eOut->Dest->KMer, gen_array_size(&intersection), midSeqLen + eIn->SeqLen + eOut->SeqLen + 1, kmetRead, &newEdge);
										if (ret == ERR_SUCCESS) {
											putchar('E');
											++numberOfCreatedEdges;
											edgesCreatedPerOutEdge[j]++;
											newEdge->SeqLen = eIn->SeqLen + midSeqLen + eOut->SeqLen;
											ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
											if (ret == ERR_SUCCESS) {
												memcpy(newEdge->Seq, eIn->Seq, eIn->SeqLen*sizeof(char));
												memcpy(newEdge->Seq + eIn->SeqLen, midSeq, midSeqLen*sizeof(char));
												memcpy(newEdge->Seq + eIn->SeqLen + midSeqLen, eOut->Seq, eOut->SeqLen*sizeof(char));
												newEdge->Seq[newEdge->SeqLen] = '\0';
												ret = read_info_assign(&newEdge->ReadInfo, &intersection);
												// TODO: Remove the reads from the edges that are now avoided
											}
										}

										if (ret == ERR_ALREADY_EXISTS)
											ret = ERR_SUCCESS;
									}

									dym_array_clear_READ_INFO_ENTRY(&intersection);
								}

								if (ret != ERR_SUCCESS)
									break;
							}

							if (numberOfCreatedEdges == w->degreeOut)
								dym_array_push_back_no_alloc_PKMER_EDGE(&edgesToDelete, eIn);

							if (ret != ERR_SUCCESS)
								break;
						}

						dym_array_finit_READ_INFO_ENTRY(&intersection);
						for (size_t i = 0; i < w->degreeOut; ++i) {
							PKMER_EDGE eOut = kmer_vertex_get_succ_edge(w, i);

							if (edgesCreatedPerOutEdge[i] == v->DegreeIn)
								dym_array_push_back_no_alloc_PKMER_EDGE(&edgesToDelete, eOut);
						}

						for (size_t i = 0; i < gen_array_size(&edgesToDelete); ++i) {
							kmer_graph_delete_edge(Graph, *dym_array_item_PKMER_EDGE(&edgesToDelete, i));
							putchar('D');
						}
					}

					dym_array_finit_PKMER_EDGE(&edgesToDelete);
					utils_free(edgesCreatedPerOutEdge);
				}
			}
		}

		if (ret == ERR_SUCCESS)
			ret = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	printf("\n");

	return ret;
}

ERR_VALUE kmer_graph_detect_uncertainities(PKMER_GRAPH Graph)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE e = NULL;
	PKMER_VERTEX v = Graph->StartingVertex;

	ret = ERR_SUCCESS;
	e = _get_refseq_edge(v);
	v = e->Dest;
	while (v != Graph->EndingVertex) {
		if (v->degreeOut == 2) {
			PKMER_VERTEX start = v;
			PKMER_EDGE rsFirstEdge = NULL;
			PKMER_EDGE readFirstEdge = NULL;

			for (size_t i = 0; i < v->degreeOut; ++i) {
				e = kmer_vertex_get_succ_edge(v, i);
				if (e->Type == kmetReference)
					rsFirstEdge = e;

				if (e->Type == kmetRead)
					readFirstEdge = e;
			}

			e = rsFirstEdge;
			v = rsFirstEdge->Dest;			
			if (v->DegreeIn == 1) {
				if (v->degreeOut != 1)
					continue;

				e = _get_refseq_edge(v);
				v = e->Dest;
			}

			if (v->DegreeIn == 2) {
				PKMER_VERTEX rsEnd = v;
				PKMER_EDGE rsLastEdge = e;
				char *s1 = NULL;
				size_t s1Len = 0;

				ret = _capture_refseq(rsFirstEdge, rsLastEdge, &s1, &s1Len);
				if (ret == ERR_SUCCESS) {
					e = readFirstEdge;
					v = readFirstEdge->Dest;
					if (v->Type == kmvtRead) {
						if (v->degreeOut != 1 || v->DegreeIn != 1) {
							e = _get_refseq_edge(start);
							v = e->Dest;
							continue;
						}

						e = kmer_vertex_get_succ_edge(v, 0);
						v = e->Dest;
					}

					if (v == rsEnd) {
						PKMER_EDGE readLastEdge = e;
						char *s2 = NULL;
						size_t s2Len = 0;

						ret = _capture_refseq(readFirstEdge, readLastEdge, &s2, &s2Len);
						if (ret == ERR_SUCCESS) {
							PKMER_VERTEX tmp = NULL;
							
							if (readFirstEdge->Dest == readLastEdge->Source)
								tmp = readFirstEdge->Dest;
							else tmp = rsFirstEdge->Dest;

							kmer_graph_delete_edge(Graph, rsFirstEdge);
							if (rsLastEdge != rsFirstEdge)
								kmer_graph_delete_edge(Graph, rsLastEdge);

							kmer_graph_delete_edge(Graph, readFirstEdge);
							if (readFirstEdge != readLastEdge)
								kmer_graph_delete_edge(Graph, readLastEdge);
						
							kmer_graph_delete_vertex(Graph, tmp);
							ret = kmer_graph_add_edge_ex(Graph, start->KMer, v->KMer, 0, 1, kmetVariant, &e);
							if (ret == ERR_SUCCESS) {
								e->Seq = s1;
								e->SeqLen = s1Len;
								e->Seq2 = s2;
								e->Seq2Len = s2Len;
							}
						}
					} else {
						e = _get_refseq_edge(start);
						v = e->Dest;
					}
				}
			}
		} else {
			e = _get_refseq_edge(v);
			v = e->Dest;
		}
	}

	return ret;
}
