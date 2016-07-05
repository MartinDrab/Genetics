
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
#include "refseq-storage.h"
#include "pointer_array.h"
#include "gen_dym_array.h"
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
			tmp->Type = Type;
			pointer_array_init_KMER_EDGE(&tmp->Successors, 140);
			pointer_array_init_KMER_EDGE(&tmp->Predecessors, 140);
			tmp->Finished = FALSE;
			tmp->LastRSOrder = 0;
			tmp->RefSeqPosition = 0;
			*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmp);
	}

	return ret;
}


static void _vertex_destroy(PKMER_VERTEX Vertex)
{
	pointer_array_finit_KMER_EDGE(&Vertex->Predecessors);
	pointer_array_finit_KMER_EDGE(&Vertex->Successors);
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
		tmp->Order = Vertex->Order;
		pointer_array_init_KMER_EDGE(&tmp->Successors, 140);
		ret = pointer_array_clean_copy_KMER_EDGE(&tmp->Successors, &Vertex->Successors);
		if (ret == ERR_SUCCESS) {
			pointer_array_init_KMER_EDGE(&tmp->Predecessors, 140);
			ret = pointer_array_clean_copy_KMER_EDGE(&tmp->Predecessors, &Vertex->Predecessors);
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

static ERR_VALUE _edge_create(PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, const long Weight, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Source = Source;
		tmp->Dest = Dest;
		tmp->Type = Type;
		tmp->Weight = Weight;
		tmp->Order = 0;
		tmp->Seq = NULL;
		tmp->SeqLen = 0;
		tmp->Seq2 = NULL;
		tmp->Seq2Len = 0;
		read_info_init(&tmp->ReadInfo);
		tmp->MarkedForDelete = FALSE;
		dym_array_init_size_t(&tmp->Paths, 140);
		*Edge = tmp;
	}

	return ret;
}


static void _edge_destroy(PKMER_EDGE Edge)
{
	dym_array_finit_size_t(&Edge->Paths);
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

	ret = _edge_create(Edge->Source, Edge->Dest, Edge->Type, Edge->Weight, &tmp);
	if (ret == ERR_SUCCESS) {
		tmp->MarkedForDelete = Edge->MarkedForDelete;
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
					tmp->Order = Edge->Order;
					ret = dym_array_push_back_array_size_t(&tmp->Paths, &Edge->Paths);
					if (ret == ERR_SUCCESS)
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
	fprintf(Stream, "\\nIN=%u; OUT=%u; O=%u\",style=filled,color=%s]", kmer_vertex_in_degree(v), kmer_vertex_out_degree(v), v->Order, colors[v->Type]);
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
			fprintf(Stream, ",label=\"W: %li; L %u; P: %u\";", e->Weight, (uint32_t)e->SeqLen + 1, (uint32_t)gen_array_size(&e->Paths));
			break;
		case kmetRead:
			fprintf(Stream, "red");
			fprintf(Stream, ",label=\"W: %li; L %u; P:%u\"", e->Weight, e->SeqLen + 1, (uint32_t)gen_array_size(&e->Paths));
			break;
		case kmetVariant:
			fprintf(Stream, "blue");
			fprintf(Stream, ",label=\"W: %li; L %u; P: %u\\n1: %s\\n2: %s\"", e->Weight, e->SeqLen + 1, (uint32_t)gen_array_size(&e->Paths), e->Seq, e->Seq2);
			break;
	}

	fprintf(Stream, "];\n");

	return;
}


static boolean _kmer_vertex_no_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = TRUE;

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		ret = (kmer_vertex_get_succ_edge(Vertex, i))->Type != kmetRead;
		if (!ret)
			break;
	}

	if (ret) {
		for (size_t i = 0; i < kmer_vertex_in_degree(Vertex); ++i) {
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
		memset(tmpGraph, 0, sizeof(KMER_GRAPH));
		tmpGraph->NumberOfEdges = 0;
		tmpGraph->NumberOfVertices = 0;
		tmpGraph->KMerSize = KMerSize;
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
	fprintf(Stream, "\t/* number of reference edges: %u */\n", Graph->TypedEdgeCount[kmetReference]);
	fprintf(Stream, "\t/* number of read edges: %u */\n", Graph->TypedEdgeCount[kmetRead]);
	fprintf(Stream, "\t/* number of variant edges: %u */\n", Graph->TypedEdgeCount[kmetVariant]);
	kmer_table_print(Stream, Graph->VertexTable);
	kmer_edge_table_print(Stream, Graph->EdgeTable);
	fprintf(Stream, "}\n");
	fflush(Stream);

	return;
}


ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph)
{
	PKMER_VERTEX v = NULL;
	PKMER_TABLE_ENTRY e = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		boolean last = FALSE;
		
		e = iter;
		last = (kmer_table_next(Graph->VertexTable, e, &iter) == ERR_NO_MORE_ENTRIES);
		v = (PKMER_VERTEX)e->Data;
		if (v->Type == kmvtRefSeqMiddle && kmer_vertex_in_degree(v) == 1 && kmer_vertex_out_degree(v) == 1) {
			const KMER_EDGE *inEdge = kmer_vertex_get_pred_edge(v, 0);
			const KMER_EDGE *outEdge = kmer_vertex_get_succ_edge(v, 0);

			if (inEdge->Type != kmetVariant && outEdge->Type != kmetVariant)
				kmer_graph_delete_vertex(Graph, v);
		}

		if (last)
			break;
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	if (ret == ERR_SUCCESS) {
		ret = kmer_table_first(Graph->VertexTable, &iter);
		while (ret == ERR_SUCCESS) {
			boolean last = FALSE;

			e = iter;
			last = (kmer_table_next(Graph->VertexTable, e, &iter) == ERR_NO_MORE_ENTRIES);
			v = (PKMER_VERTEX)e->Data;
			if (v->Type == kmvtRead && kmer_vertex_in_degree(v) == 1 && kmer_vertex_out_degree(v) == 1) {
				const KMER_EDGE *inEdge = kmer_vertex_get_pred_edge(v, 0);
				const KMER_EDGE *outEdge = kmer_vertex_get_succ_edge(v, 0);

				if (inEdge->Type != kmetVariant && outEdge->Type != kmetVariant)
					kmer_graph_delete_vertex(Graph, v);
			}

			if (last)
				break;
		}
	}

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
		err = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
		if (e->Type == kmetRead && read_info_get_count(&e->ReadInfo) <= Threshold)
			kmer_graph_delete_edge(Graph, e);
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
			if (kmer_vertex_in_degree(v) == 0) {
				deleted = TRUE;
				while (kmer_vertex_out_degree(v) > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_succ_edge(v, 0));
			} else if (kmer_vertex_out_degree(v) == 0) {
				deleted = TRUE;
				while (kmer_vertex_in_degree(v) > 0)
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



static PKMER_EDGE _get_refseq_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_succ_edge(Vertex, i);

		if (tmp->Type == kmetReference) {
			ret = tmp;
			break;
		}
	}

	return ret;
}


static PKMER_EDGE _get_in_refseq_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < kmer_vertex_in_degree(Vertex); ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_pred_edge(Vertex, i);

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

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		PKMER_EDGE e = kmer_vertex_get_succ_edge(Vertex, i);

		ret = (e->Type == kmetRead);
		if (ret)
			break;
	}

	return ret;
}


ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(KMer, Type, &v);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_table_insert(Graph->VertexTable, v->KMer, v);
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


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _edge_create(Source, Dest, Type, weight, &edge);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_edge_table_insert(Graph->EdgeTable, Source->KMer, Dest->KMer, edge);
			if (ret == ERR_SUCCESS) {
				ret = pointer_array_reserve_KMER_EDGE(&Source->Successors, pointer_array_size(&Source->Successors) + 1);
				if (ret == ERR_SUCCESS) {
					ret = pointer_array_reserve_KMER_EDGE(&Dest->Predecessors, pointer_array_size(&Dest->Predecessors) + 1);
					if (ret == ERR_SUCCESS) {
						pointer_array_push_back_no_alloc_KMER_EDGE(&Source->Successors, edge);
						pointer_array_push_back_no_alloc_KMER_EDGE(&Dest->Predecessors, edge);
						Graph->NumberOfEdges++;
						Graph->TypedEdgeCount[Type]++;
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
			*Edge = (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source->KMer, Dest->KMer);
	}

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


ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PPOINTER_ARRAY_KMER_VERTEX VertexArray)
{
	DYM_ARRAY tmp;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	dym_array_create(&tmp, 140);
	ret = kmer_table_get_multiple(Graph->VertexTable, KMer, &tmp);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < dym_array_size(&tmp); ++i) {
			ret = pointer_array_push_back_KMER_VERTEX(VertexArray, (PKMER_VERTEX)dym_array_get(&tmp, i));
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

	if (kmer_vertex_in_degree(Vertex) <= 1 && kmer_vertex_out_degree(Vertex) <= 1) {
		if (kmer_vertex_out_degree(Vertex) == 1) {
			outEdge = kmer_vertex_get_succ_edge(Vertex, 0);
			succVertex = outEdge->Dest;
		}

		if (kmer_vertex_in_degree(Vertex) == 1) {
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
	EKMerEdgeType edgeType = Edge->Type;

	if (Graph->DeleteEdgeCallback != NULL)
		Graph->DeleteEdgeCallback(Graph, Edge, Graph->DeleteEdgeCallbackContext);

	err = kmer_edge_table_delete(Graph->EdgeTable, source->KMer, dest->KMer);
	if (err == ERR_SUCCESS) {
		pointer_array_remove_by_item_fast_KMER_EDGE(&source->Successors, Edge);
		pointer_array_remove_by_item_fast_KMER_EDGE(&dest->Predecessors, Edge);
		--Graph->TypedEdgeCount[edgeType];
		--Graph->NumberOfEdges;
	}

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
			boolean mfd = (Source->MarkedForDelete || Dest->MarkedForDelete);
			size_t pathCount = max(gen_array_size(&Source->Paths), gen_array_size(&Dest->Paths));
			long weight = max(Source->Weight, Dest->Weight);
			EKMerEdgeType type = (Source->Type == kmetReference && Dest->Type == kmetReference) ? kmetReference : kmetRead;
			PKMER_EDGE newEdge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, u, w, weight, 0, type, &newEdge);
			if (ret == ERR_SUCCESS) {
				newEdge->MarkedForDelete = mfd;
				newEdge->SeqLen = Source->SeqLen + 1 + Dest->SeqLen;
				ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
				if (ret == ERR_SUCCESS) {
					memcpy(newEdge->Seq, Source->Seq, Source->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
					memcpy(newEdge->Seq + Source->SeqLen + 1, Dest->Seq, Dest->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen + 1 + Dest->SeqLen] = '\0';
					ret = read_info_merge(&newEdge->ReadInfo, &Source->ReadInfo, &Dest->ReadInfo);
					if (ret == ERR_SUCCESS) {
						kmer_graph_delete_edge(Graph, Source);
						kmer_graph_delete_edge(Graph, Dest);
					}
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, newEdge);
			}
		} else ret = ERR_TRIANGLE;
	} else ret = (u == w) ? ERR_PRED_IS_SUCC : ERR_NOT_ADJACENT;

	return ret;
}


ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, PPOINTER_ARRAY_FOUND_SEQUENCE SeqArray)
{
	REFSEQ_STORAGE rsStorage;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX v = Graph->StartingVertex;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	size_t edgeIndex = 0;
	PKMER_EDGE edge = NULL;
	GEN_ARRAY_size_t edgeIndices;
	POINTER_ARRAY_KMER_EDGE edges;
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT variantStack;

	ret = ERR_SUCCESS;
	rs_storage_init(&rsStorage);
	dym_array_init_FOUND_SEQUENCE_VARIANT(&variantStack, 140);
	dym_array_init_size_t(&edgeIndices, 140);
	pointer_array_init_KMER_EDGE(&edges, 140);
		do {
			if (edgeIndex == kmer_vertex_out_degree(v)) {
				v->Finished = FALSE;
				if (v == Graph->StartingVertex)
					break;

				edge = *pointer_array_pop_back_KMER_EDGE(&edges);
				edgeIndex = *dym_array_pop_back_size_t(&edgeIndices);
				v = edge->Source;
				if (edge->Type == kmetVariant) {
					rs_storage_remove(&rsStorage, 2);
					dym_array_pop_back_FOUND_SEQUENCE_VARIANT(&variantStack);
				} else rs_storage_remove(&rsStorage, edge->SeqLen + 1);
			} else {
				edge = kmer_vertex_get_succ_edge(v, edgeIndex);
				if (edge->Type == kmetVariant) {
					FOUND_SEQUENCE_VARIANT variant;

					ret = rs_storage_add(&rsStorage, '?');
					if (ret != ERR_SUCCESS)
						break;

					variant.RefSeqStart = edge->Source->RefSeqPosition + 1;
					variant.RefSeqEnd = edge->Dest->RefSeqPosition;
					variant.Seq1Type = edge->SeqType;
					variant.Seq1 = edge->Seq;
					variant.Seq1Len = edge->SeqLen;
					variant.Seq2Type = edge->Seq2Type;
					variant.Seq2 = edge->Seq2;
					variant.Seq2Len = edge->Seq2Len;
					variant.Seq1Weight = edge->Seq1Weight;
					variant.Seq2Weight = edge->Seq2Weight;
					ret = dym_array_push_back_FOUND_SEQUENCE_VARIANT(&variantStack, variant);
					if (ret != ERR_SUCCESS)
						break;
				} else {
					ret = rs_storage_add_seq(&rsStorage, edge->Seq, edge->SeqLen);
					if (ret != ERR_SUCCESS)
						break;
				}

				if (edge->Dest == Graph->EndingVertex) {
					PFOUND_SEQUENCE fs = NULL;
					const size_t vCount = gen_array_size(&variantStack);

					ret = found_sequence_alloc(rsStorage.Sequence, rsStorage.ValidLength, vCount, &fs);
					if (ret == ERR_SUCCESS) {
						for (size_t i = 0; i < vCount; ++i) {
							ret = found_sequence_set_variant(fs, i, dym_array_item_FOUND_SEQUENCE_VARIANT(&variantStack, i));
							if (ret != ERR_SUCCESS)
								break;
						}

						for (size_t i = 0; i < pointer_array_size(&edges); ++i) {
							ret = dym_array_push_back_size_t(&(edges.Data[i]->Paths), (size_t)fs);
							if (ret != ERR_SUCCESS)
								break;
						}

						if (ret == ERR_SUCCESS) {
							ret = dym_array_push_back_size_t(&edge->Paths, (size_t)fs);
							if (ret == ERR_SUCCESS) {
								ret = pointer_array_push_back_FOUND_SEQUENCE(SeqArray, fs);
								if (ret == ERR_SUCCESS) {
									rs_storage_remove(&rsStorage, edge->SeqLen);
									++edgeIndex;
									if (pointer_array_size(SeqArray) >= 1000)
										break;
								}

								if (ret != ERR_SUCCESS)
									dym_array_pop_back_size_t(&edge->Paths);
							}

							if (ret != ERR_SUCCESS) {
								for (size_t i = 0; i < pointer_array_size(&edges); ++i)
									dym_array_pop_back_size_t(&(edges.Data[i]->Paths));
							}
						}

						if (ret != ERR_SUCCESS)
							found_sequence_free(fs);
					}
				} else {
					v = edge->Dest;
					if (v != edge->Source && !v->Finished) {
						edge->Source->Finished = TRUE;
						ret = rs_storage_add(&rsStorage, kmer_get_base(v->KMer, kmerSize - 1));
						if (ret != ERR_SUCCESS)
							break;

						ret = dym_array_push_back_size_t(&edgeIndices, (edgeIndex + 1));
						if (ret != ERR_SUCCESS)
							break;

						ret = pointer_array_push_back_KMER_EDGE(&edges, edge);
						if (ret != ERR_SUCCESS)
							break;

						edgeIndex = 0;
					} else {
						v = edge->Source;
						if (edge->Type == kmetVariant) {
							rs_storage_remove(&rsStorage, 1);
							dym_array_pop_back_FOUND_SEQUENCE_VARIANT(&variantStack);
						} else rs_storage_remove(&rsStorage, edge->SeqLen);

						++edgeIndex;
					}
				}
			}
		} while (ret == ERR_SUCCESS && v != Graph->StartingVertex);

	pointer_array_finit_KMER_EDGE(&edges);
	dym_array_finit_size_t(&edgeIndices);
	dym_array_finit_FOUND_SEQUENCE_VARIANT(&variantStack);
	if (ret != ERR_SUCCESS) {
		for (size_t i = 0; i < pointer_array_size(SeqArray); ++i)
			found_sequence_free(*pointer_array_item_FOUND_SEQUENCE(SeqArray, i));

		pointer_array_clear_FOUND_SEQUENCE(SeqArray);
	}

	rs_storage_finit(&rsStorage);

	return ret;
}


void kmer_graph_delete_seqs(PKMER_GRAPH Graph, PPOINTER_ARRAY_FOUND_SEQUENCE FoundSequences, const size_t Threshold)
{
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (ret == ERR_SUCCESS) {
		PKMER_EDGE e = (PKMER_EDGE)iter->Data;

		ret = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
		if (e->Type == kmetReference && (read_info_get_count(&e->ReadInfo) <= Threshold || e->MarkedForDelete) &&
			gen_array_size(&e->Paths) < pointer_array_size(FoundSequences)) {
			PFOUND_SEQUENCE fs = NULL;
			
			for (size_t i = 0; i < gen_array_size(&e->Paths); ++i) {
				fs = *(PFOUND_SEQUENCE *)dym_array_item_size_t(&e->Paths, i);
				if (pointer_array_remove_by_item_fast_FOUND_SEQUENCE(FoundSequences, fs))
					found_sequence_free(fs);
			}

			kmer_graph_delete_edge(Graph, e);
		}
	}

	return;
}


ERR_VALUE kmer_vertex_get_certain_edges(const KMER_VERTEX *Vertex, const EKMerEdgeType EdgeType, const boolean Incomming, PPOINTER_ARRAY_KMER_EDGE Array)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const PKMER_EDGE *edgeArray = (Incomming) ? (Vertex->Predecessors.Data) : (Vertex->Successors.Data);
	const size_t count = pointer_array_size((Incomming) ? (&Vertex->Predecessors) : (&Vertex->Successors));

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < count; ++i) {
		const KMER_EDGE *e = *edgeArray;
	
		if (e->Type == EdgeType)
			ret = pointer_array_push_back_KMER_EDGE(Array, e);
		
		if (ret != ERR_SUCCESS)
			break;

		++edgeArray;
	}

	return ret;
}


static boolean _is_triangle(const KMER_GRAPH *Graph, const KMER_VERTEX *v1, const KMER_VERTEX *v2, const KMER_VERTEX *v3)
{
	boolean ret = FALSE;

	ret = (
		kmer_graph_get_edge(Graph, v1->KMer, v2->KMer) != NULL &&
		kmer_graph_get_edge(Graph, v2->KMer, v3->KMer) != NULL &&
		kmer_graph_get_edge(Graph, v1->KMer, v3->KMer) != NULL);

	return ret;
}


void kmer_graph_resolve_db_triangles(PKMER_GRAPH Graph, const size_t Threshold)
{
	PKMER_TABLE_ENTRY entry = NULL;
	ERR_VALUE ret = ERR_SUCCESS;

	ret = kmer_table_first(Graph->VertexTable, &entry);
	while (ret == ERR_SUCCESS) {
		PKMER_VERTEX v1 = (PKMER_VERTEX)entry->Data;

		if (v1->Type == kmvtRefSeqMiddle) {
			do {
				PKMER_EDGE e = NULL;
				PKMER_VERTEX v2 = NULL;
				PKMER_VERTEX v3 = NULL;
				PKMER_VERTEX v4 = NULL;

				e = _get_refseq_edge(v1);
				if (e == NULL)
					break;

				v2 = e->Dest;
				if (v2 == Graph->EndingVertex)
					break;

				e = _get_refseq_edge(v2);
				if (e == NULL)
					break;

				v3 = e->Dest;
				if (v3 == Graph->EndingVertex)
					break;

				e = _get_refseq_edge(v3);
				if (e == NULL)
					break;

				v4 = e->Dest;
				if (v4 == Graph->EndingVertex)
					break;

				e = kmer_graph_get_edge(Graph, v2->KMer, v3->KMer);
				if (e != NULL && (e->Weight == 0 || e->Weight < (long)Threshold) && _is_triangle(Graph, v1, v2, v3) && _is_triangle(Graph, v2, v3, v4))
					kmer_graph_delete_edge(Graph, e);
			} while (FALSE);

		}

		ret = kmer_table_next(Graph->VertexTable, entry, &entry);
	}


	return;
}

static ERR_VALUE _capture_refseq(const KMER_EDGE *Start, const POINTER_ARRAY_KMER_EDGE *RSEdges, const KMER_EDGE *End, char **Seq, size_t *SeqLen)
{
	REFSEQ_STORAGE rsStorage;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t kmerSize = kmer_get_size(Start->Dest->KMer);

	rs_storage_init(&rsStorage);
	ret = rs_storage_add_edge(&rsStorage, Start, FALSE);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < pointer_array_size(RSEdges); ++i) {
			const KMER_EDGE *e = *pointer_array_const_item_KMER_EDGE(RSEdges, i);
				
			ret = rs_storage_add_edge(&rsStorage, e, FALSE);
			if (ret != ERR_SUCCESS)
				break;
		}

		if (ret == ERR_SUCCESS) {
			ret = rs_storage_add_edge(&rsStorage, End, FALSE);
			if (ret == ERR_SUCCESS) {
				rs_storage_remove(&rsStorage, 1);
				ret = rs_storage_create_string(&rsStorage, Seq);
				*SeqLen = rsStorage.ValidLength;
			}
		}
	}

	rs_storage_finit(&rsStorage);

	return ret;
}


static ERR_VALUE _remove_read_info_from_edges(PKMER_GRAPH Graph, PKMER_EDGE Edge1, PKMER_EDGE Edge2, PPOINTER_ARRAY_KMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY Intersection)
{
	size_t distance = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_READ_INFO_ENTRY diff;

	dym_array_init_READ_INFO_ENTRY(&diff, 140);
	ret = read_info_diff(&Edge1->ReadInfo, Intersection, &diff, distance);
	if (ret == ERR_SUCCESS) {
		Edge1->Weight -= (read_info_get_count(&Edge1->ReadInfo) - gen_array_size(&diff));
		if (Edge1->Weight < 0)
			Edge1->Weight = 0;

		distance += Edge1->SeqLen + 1;
		ret = read_info_assign(&Edge1->ReadInfo, &diff);
	}


	if (ret == ERR_SUCCESS) {
		for (size_t k = 0; k < gen_array_size(Array); ++k) {
			PKMER_EDGE tmp = *pointer_array_item_KMER_EDGE(Array, k);

			dym_array_clear_READ_INFO_ENTRY(&diff);
			ret = read_info_diff(&tmp->ReadInfo, Intersection, &diff, distance);
			if (ret == ERR_SUCCESS) {
				tmp->Weight -= (read_info_get_count(&tmp->ReadInfo) - gen_array_size(&diff));
				if (tmp->Weight < 0)
					tmp->Weight = 0;

				distance += tmp->SeqLen + 1;
				ret = read_info_assign(&tmp->ReadInfo, &diff);
			}

			if (ret != ERR_SUCCESS)
				break;
		}

		if (ret == ERR_SUCCESS) {
			dym_array_clear_READ_INFO_ENTRY(&diff);
			ret = read_info_diff(&Edge2->ReadInfo, Intersection, &diff, distance);
			if (ret == ERR_SUCCESS) {
				Edge2->Weight -= (read_info_get_count(&Edge2->ReadInfo) - gen_array_size(&diff));
				if (Edge2->Weight < 0)
					Edge2->Weight = 0;

				ret = read_info_assign(&Edge2->ReadInfo, &diff);
			}
		}
	}

	dym_array_finit_READ_INFO_ENTRY(&diff);

	return ret;
}


ERR_VALUE kmer_graph_find_rs_read_edges(const KMER_GRAPH *Graph, PPOINTER_ARRAY_KMER_EDGE InEdges, PPOINTER_ARRAY_KMER_EDGE OutEdges)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE_ENTRY iter = NULL;

	ret = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (ret == ERR_SUCCESS) {
		const KMER_EDGE *e = (const KMER_EDGE *)iter->Data;
		
		if (e->Type == kmetRead) {
			if (e->Source->Type == kmvtRefSeqMiddle)
				ret = pointer_array_push_back_KMER_EDGE(OutEdges, e);

			if (ret == ERR_SUCCESS) {
				if (e->Dest->Type == kmvtRefSeqMiddle)
					ret = pointer_array_push_back_KMER_EDGE(InEdges, e);
			}
		}

		if (ret == ERR_SUCCESS)
			ret = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	return ret;
}


ERR_VALUE kmer_graph_process_rs_read_edges(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_EDGE InEdges, PPOINTER_ARRAY_KMER_EDGE OutEdges, const size_t Threshold)
{
	boolean added = FALSE;
	PKMER_EDGE eIn = NULL;
	PKMER_EDGE eOut = NULL;
	size_t readDistance;
	POINTER_ARRAY_KMER_EDGE rsEdges;
	POINTER_ARRAY_KMER_EDGE edgesToDelete;
	GEN_ARRAY_READ_INFO_ENTRY intersection;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	POINTER_ARRAY_KMER_EDGE tmpInEdges;
	POINTER_ARRAY_KMER_EDGE tmpOutEdges;

	pointer_array_init_KMER_EDGE(&tmpInEdges, 140);
	pointer_array_init_KMER_EDGE(&tmpOutEdges, 140);
	pointer_array_init_KMER_EDGE(&edgesToDelete, 140);
	pointer_array_init_KMER_EDGE(&rsEdges, 140);
	dym_array_init_READ_INFO_ENTRY(&intersection, 140);
	do {
		added = FALSE;
		for (size_t i = 0; i < pointer_array_size(InEdges); ++i) {
			eIn = *pointer_array_item_KMER_EDGE(InEdges, i);
			if (eIn == NULL)
				continue;

			for (size_t j = 0; j < pointer_array_size(OutEdges); ++j) {
				eOut = *pointer_array_item_KMER_EDGE(OutEdges, j);
				if (eOut == NULL)
					continue;

				if (eOut == eIn)
					continue;

				if (eIn->Dest->RefSeqPosition > eOut->Source->RefSeqPosition)
					continue;

				readDistance = eOut->Source->RefSeqPosition - eIn->Dest->RefSeqPosition + 1;
				dym_array_clear_READ_INFO_ENTRY(&intersection);
				ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, TRUE, eIn->SeqLen + readDistance);
				if (ret == ERR_SUCCESS) {
					const size_t weight = gen_array_size(&intersection);

					if (weight > Threshold) {
						const KMER_VERTEX *v = eIn->Dest;

						pointer_array_clear_KMER_EDGE(&rsEdges);
						while (ret == ERR_SUCCESS && v != eOut->Source) {
							const KMER_EDGE *e = _get_refseq_edge(v);

							v = e->Dest;
							ret = pointer_array_push_back_KMER_EDGE(&rsEdges, e);
						}

						if (ret == ERR_SUCCESS) {
							char *seq = NULL;
							size_t seqLen = 0;

							ret = _capture_refseq(eIn, &rsEdges, eOut, &seq, &seqLen);
							if (ret == ERR_SUCCESS) {
								PKMER_EDGE newEdge = NULL;

								ret = kmer_graph_add_edge_ex(Graph, eIn->Source, eOut->Dest, weight, 0, kmetRead, &newEdge);
								if (ret == ERR_SUCCESS) {
									ret = read_info_assign(&newEdge->ReadInfo, &intersection);
									if (ret == ERR_SUCCESS) {
										newEdge->Seq = seq;
										newEdge->SeqLen = seqLen;
										newEdge->SeqType = kmetRead;
									}
								}

								if (ret == ERR_ALREADY_EXISTS && newEdge->SeqLen == seqLen && memcmp(seq, newEdge->Seq, seqLen*sizeof(char)) == 0)
									ret = read_info_union(&newEdge->ReadInfo, &intersection);

								if (ret == ERR_SUCCESS) {
									ret = _remove_read_info_from_edges(Graph, eIn, eOut, &rsEdges, &intersection);
									if (ret == ERR_SUCCESS) {
										PKMER_EDGE rsLastEdge = _get_in_refseq_edge(eIn->Dest);
										PKMER_EDGE rsNextEdge = _get_refseq_edge(eOut->Source);

										dym_array_clear_READ_INFO_ENTRY(&intersection);
										if (read_info_get_count(&rsNextEdge->ReadInfo) > Threshold ||
											read_info_get_count(&rsLastEdge->ReadInfo) <= Threshold ||
											read_info_get_count(&eOut->ReadInfo) <= Threshold) {
											ret = read_info_intersection(&rsLastEdge->ReadInfo, &eOut->ReadInfo, &intersection, FALSE, rsLastEdge->SeqLen + readDistance);
											if (ret == ERR_SUCCESS && gen_array_size(&intersection) > 0)
												ret = _remove_read_info_from_edges(Graph, rsLastEdge, eOut, &rsEdges, &intersection);

											dym_array_clear_READ_INFO_ENTRY(&intersection);
											if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eOut))
												ret = pointer_array_push_back_KMER_EDGE(&edgesToDelete, eOut);
										}

										if (read_info_get_count(&rsNextEdge->ReadInfo) <= Threshold ||
											read_info_get_count(&rsLastEdge->ReadInfo) > Threshold ||
											read_info_get_count(&eIn->ReadInfo) <= Threshold) {
											ret = read_info_intersection(&eIn->ReadInfo, &rsNextEdge->ReadInfo, &intersection, FALSE, eIn->SeqLen + readDistance);
											if (ret == ERR_SUCCESS && gen_array_size(&intersection) > 0)
												ret = _remove_read_info_from_edges(Graph, eIn, rsNextEdge, &rsEdges, &intersection);

											dym_array_clear_READ_INFO_ENTRY(&intersection);
											if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eIn))
												ret = pointer_array_push_back_KMER_EDGE(&edgesToDelete, eIn);
										}

										if (newEdge->Source->Type == kmvtRefSeqMiddle) {
											pointer_array_push_back_KMER_EDGE(&tmpOutEdges, newEdge);
											added = TRUE;
										}

										if (newEdge->Dest->Type == kmvtRefSeqMiddle) {
											pointer_array_push_back_KMER_EDGE(&tmpInEdges, newEdge);
											added = TRUE;
										}
									}
								}

								if (ret != ERR_SUCCESS)
									utils_free(seq);

								if (ret == ERR_ALREADY_EXISTS)
									ret = ERR_SUCCESS;
							}
						}
					}
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			if (ret != ERR_SUCCESS)
				break;
		}

		if (added) {
			pointer_array_push_back_array_KMER_EDGE(InEdges, &tmpInEdges);
			pointer_array_clear_KMER_EDGE(&tmpInEdges);
			pointer_array_push_back_array_KMER_EDGE(OutEdges, &tmpOutEdges);
			pointer_array_clear_KMER_EDGE(&tmpOutEdges);
		}
	} while (ret == ERR_SUCCESS && added);

	for (size_t i = 0; i < pointer_array_size(&edgesToDelete); ++i) {
		PKMER_EDGE e = *pointer_array_item_KMER_EDGE(&edgesToDelete, i);

		kmer_graph_delete_edge(Graph, e);
	}

	dym_array_finit_READ_INFO_ENTRY(&intersection);
	pointer_array_finit_KMER_EDGE(&rsEdges);
	pointer_array_finit_KMER_EDGE(&edgesToDelete);
	pointer_array_init_KMER_EDGE(&tmpOutEdges, 140);
	pointer_array_init_KMER_EDGE(&tmpInEdges, 140);

	return ret;
}


ERR_VALUE kmer_graph_connect_reads_by_pairs(PKMER_GRAPH Graph, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, size_t *ChangeCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_READ_INFO_ENTRY intersection;
	const KMER_EDGE *e = NULL;
	PKMER_EDGE rsLastEdge = NULL;
	PKMER_EDGE rsNextEdge = NULL;
	const KMER_EDGE *nextRsLastEdge = NULL;
	POINTER_ARRAY_KMER_EDGE edgesToDelete;
	POINTER_ARRAY_KMER_EDGE rsEdges;

	ret = ERR_SUCCESS;
	pointer_array_init_KMER_EDGE(&rsEdges, 140);
	dym_array_init_READ_INFO_ENTRY(&intersection, 140);
	pointer_array_init_KMER_EDGE(&edgesToDelete, 140);
	for (size_t h = 0; h < gen_array_size(PairArray); ++h) {
		boolean edgeCreated = FALSE;
		KMER_EDGE_PAIR pair = *dym_array_item_KMER_EDGE_PAIR(PairArray, h);
		PKMER_EDGE eIn = pair.U;
		PKMER_EDGE eOut = pair.V;

		rsLastEdge = _get_in_refseq_edge(eIn->Dest);
		if (rsLastEdge == NULL)
			continue;

		rsNextEdge = _get_refseq_edge(eOut->Source);
		if (rsNextEdge == NULL || eIn == eOut)
			continue;

		pointer_array_clear_KMER_EDGE(&rsEdges);
		const KMER_VERTEX *v = eIn->Dest;
		while (ret == ERR_SUCCESS && v != eOut->Source) {
			const KMER_EDGE *e = _get_refseq_edge(v);
			ret = pointer_array_push_back_KMER_EDGE(&rsEdges, e);
			v = e->Dest;
		}

		if (ret != ERR_SUCCESS)
			break;

		if (pointer_array_size(&rsEdges) != pair.ReadDistance)
			continue;

		if (pointer_array_contains_KMER_EDGE(&edgesToDelete, eIn))
			continue;

		if (pointer_array_contains_KMER_EDGE(&edgesToDelete, eOut))
			continue;

		if (ret == ERR_SUCCESS) {
			ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, TRUE, eIn->SeqLen + 1 + pair.ReadDistance);
			if (ret == ERR_SUCCESS && gen_array_size(&intersection) > Threshold) {
				char *seq = NULL;
				size_t seqLen = 0;

				ret = _capture_refseq(eIn, &rsEdges, eOut, &seq, &seqLen);
				if (ret == ERR_SUCCESS) {
					PKMER_EDGE newEdge = NULL;

					if (kmer_equal(eIn->Source->KMer, eOut->Dest->KMer))
						ret = ERR_ALREADY_EXISTS;

					if (ret == ERR_SUCCESS) {
						ret = kmer_graph_add_edge_ex(Graph, eIn->Source, eOut->Dest, gen_array_size(&intersection), seqLen + 1, kmetRead, &newEdge);
						if (ret == ERR_SUCCESS) {
							*ChangeCount++;
							edgeCreated = TRUE;
							newEdge->Seq = seq;
							newEdge->SeqLen = seqLen;
							ret = read_info_assign(&newEdge->ReadInfo, &intersection);
							if (ret == ERR_SUCCESS)
								ret = _remove_read_info_from_edges(Graph, eIn, eOut, &rsEdges, &intersection);
						}

						if (ret == ERR_ALREADY_EXISTS && newEdge->SeqLen == seqLen && memcmp(newEdge->Seq, seq, seqLen) == 0) {
							*ChangeCount++;
							edgeCreated = TRUE;
							ret = _remove_read_info_from_edges(Graph, eIn, eOut, &rsEdges, &intersection);
						}
					}

					if (ret != ERR_SUCCESS)
						utils_free(seq);

					if (ret == ERR_ALREADY_EXISTS)
						ret = ERR_SUCCESS;
				}

				dym_array_clear_READ_INFO_ENTRY(&intersection);
			}

			if (edgeCreated) {
				dym_array_clear_READ_INFO_ENTRY(&intersection);
				if (read_info_get_count(&rsNextEdge->ReadInfo) > Threshold ||
					read_info_get_count(&rsLastEdge->ReadInfo) <= Threshold ||
					read_info_get_count(&eOut->ReadInfo) <= Threshold) {
					ret = read_info_intersection(&rsLastEdge->ReadInfo, &eOut->ReadInfo, &intersection, FALSE, rsLastEdge->SeqLen + 1 + pair.ReadDistance);
					if (ret == ERR_SUCCESS && gen_array_size(&intersection) > 0)
						ret = _remove_read_info_from_edges(Graph, rsLastEdge, eOut, &rsEdges, &intersection);

					dym_array_clear_READ_INFO_ENTRY(&intersection);
					if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eOut))
						ret = pointer_array_push_back_KMER_EDGE(&edgesToDelete, eOut);				
				
//					rsLastEdge->MarkedForDelete = (read_info_get_count(&rsLastEdge->ReadInfo) <= Threshold);
				}

				if (read_info_get_count(&rsNextEdge->ReadInfo) <= Threshold ||
					read_info_get_count(&rsLastEdge->ReadInfo) > Threshold ||
					read_info_get_count(&eIn->ReadInfo) <= Threshold) {
					ret = read_info_intersection(&eIn->ReadInfo, &rsNextEdge->ReadInfo, &intersection, FALSE, eIn->SeqLen + 1 + pair.ReadDistance);
					if (ret == ERR_SUCCESS && gen_array_size(&intersection) > 0)
						ret = _remove_read_info_from_edges(Graph, eIn, rsNextEdge, &rsEdges, &intersection);

					dym_array_clear_READ_INFO_ENTRY(&intersection);
					if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eIn))
						ret = pointer_array_push_back_KMER_EDGE(&edgesToDelete, eIn);
				
//					rsNextEdge->MarkedForDelete = (read_info_get_count(&rsNextEdge->ReadInfo) <= Threshold);
				}
			}
		}

		if (ret != ERR_SUCCESS)
			break;
	}

	dym_array_finit_READ_INFO_ENTRY(&intersection);
	pointer_array_finit_KMER_EDGE(&rsEdges);

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	for (size_t i = 0; i < pointer_array_size(&edgesToDelete); ++i) {
		kmer_graph_delete_edge(Graph, *pointer_array_item_KMER_EDGE(&edgesToDelete, i));
		*ChangeCount++;
	}

	pointer_array_finit_KMER_EDGE(&edgesToDelete);

	size_t dummy = 0;
	kmer_graph_delete_trailing_things(Graph, &dummy);

	return ret;
}


ERR_VALUE kmer_graph_connect_reads_by_reads(PKMER_GRAPH Graph, const size_t Threshold)
{
	POINTER_ARRAY_KMER_EDGE edgesToDelete;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const kmerSize = kmer_graph_get_kmer_size(Graph);

	pointer_array_init_KMER_EDGE(&edgesToDelete, 140);
	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;

		if (v->Type == kmvtRead && kmer_vertex_in_degree(v) > 1) {
			PKMER_VERTEX w = v;
			boolean cycle = FALSE;
			size_t readDistance = 0;
			POINTER_ARRAY_KMER_EDGE midEdges;

			pointer_array_init_KMER_EDGE(&midEdges, 140);
			++readDistance;
			while (ret == ERR_SUCCESS && !cycle && w->Type == kmvtRead && kmer_vertex_out_degree(w) == 1 && (v == w || kmer_vertex_in_degree(w) == 1)) {
				PKMER_EDGE e = kmer_vertex_get_succ_edge(w, 0);
				
				readDistance += e->SeqLen;
				w = e->Dest;
				++readDistance;
				cycle = (v == w);
				ret = pointer_array_push_back_KMER_EDGE(&midEdges, e);
			}

			if (!cycle && w->Type == kmvtRead &&  kmer_vertex_out_degree(w) > 1 && (v == w || kmer_vertex_in_degree(w) == 1)) {
				boolean edgeCreated = FALSE;
				GEN_ARRAY_READ_INFO_ENTRY intersection;

				dym_array_init_READ_INFO_ENTRY(&intersection, 140);
				for (size_t i = 0; i < kmer_vertex_in_degree(v); ++i) {
					PKMER_EDGE eIn = kmer_vertex_get_pred_edge(v, i);

					edgeCreated = FALSE;
					for (size_t j = 0; j < kmer_vertex_out_degree(w); ++j) {
						PKMER_EDGE eOut = kmer_vertex_get_succ_edge(w, j);
							
						edgeCreated = FALSE;
						ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, TRUE, eIn->SeqLen + readDistance);
						if (ret == ERR_SUCCESS) {
							if (gen_array_size(&intersection) > Threshold) {
								PKMER_EDGE newEdge = NULL;
									
								ret = kmer_graph_add_edge_ex(Graph, eIn->Source, eOut->Dest, gen_array_size(&intersection), 0, kmetRead, &newEdge);
								if (ret == ERR_SUCCESS) {
									ret = _capture_refseq(eIn, &midEdges, eOut, &newEdge->Seq, &newEdge->SeqLen);
									edgeCreated = TRUE;
									ret = read_info_assign(&newEdge->ReadInfo, &intersection);
									if (ret == ERR_SUCCESS)
										ret = _remove_read_info_from_edges(Graph, eIn, eOut, &midEdges, &intersection);
								} else if (ret == ERR_ALREADY_EXISTS) {
									char *tmpSeq = NULL;
									size_t tmpSeqLen = 0;
									
									ret = _capture_refseq(eIn, &midEdges, eOut, &tmpSeq, &tmpSeqLen);
									if (ret == ERR_SUCCESS) {
										if (newEdge->SeqLen == tmpSeqLen && memcmp(newEdge->Seq, tmpSeq, tmpSeqLen) == 0) {
											edgeCreated = TRUE;
											ret = read_info_union(&newEdge->ReadInfo, &intersection);
											if (ret == ERR_SUCCESS)
												ret = _remove_read_info_from_edges(Graph, eIn, eOut, &midEdges, &intersection);
										}
											
										utils_free(tmpSeq);
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

					if (ret != ERR_SUCCESS)
						break;
				}

				dym_array_finit_READ_INFO_ENTRY(&intersection);
//				for (size_t i = 0; i < pointer_array_size(&midEdges); ++i) {
//					const KMER_EDGE *e = *pointer_array_const_item_KMER_EDGE(&midEdges, i);
//
//					if (read_info_get_count(&e->ReadInfo) <= Threshold && !pointer_array_contains_KMER_EDGE(&edgesToDelete, e))
//						pointer_array_push_back_KMER_EDGE(&edgesToDelete, e);
//				}
			}

			pointer_array_finit_KMER_EDGE(&midEdges);
		}

		if (ret == ERR_SUCCESS)
			ret = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

//	for (size_t i = 0; i < pointer_array_size(&edgesToDelete); ++i)
//		kmer_graph_delete_edge(Graph, *pointer_array_item_KMER_EDGE(&edgesToDelete, i));

	pointer_array_finit_KMER_EDGE(&edgesToDelete);

	return ret;
}

void kmer_graph_delete_backward_edges(PKMER_GRAPH Graph)
{
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (ret == ERR_SUCCESS) {
		PKMER_EDGE e = (PKMER_EDGE)iter->Data;
		const KMER_VERTEX *u = e->Source;
		const KMER_VERTEX *v = e->Dest;

		ret = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
		if (e->Type == kmetRead && u->Type == kmvtRefSeqMiddle &&
			v->Type == kmvtRefSeqMiddle && v->Order < u->Order)
			kmer_graph_delete_edge(Graph, e);
	}

	return;
}


ERR_VALUE kmer_graph_detect_uncertainities(PKMER_GRAPH Graph, boolean *Changed)
{
	REFSEQ_STORAGE s1;
	REFSEQ_STORAGE s2;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE e = NULL;
	PKMER_VERTEX v = Graph->StartingVertex;
	PKMER_TABLE_ENTRY entry = NULL;

	rs_storage_init(&s1);
	rs_storage_init(&s2);
	ret = kmer_table_first(Graph->VertexTable, &entry);
	while (ret == ERR_SUCCESS) {
		v = (PKMER_VERTEX)entry->Data;
		if (v->Type == kmvtRefSeqMiddle && kmer_vertex_out_degree(v) == 2) {
			PKMER_EDGE path1Start = kmer_vertex_get_succ_edge(v, 0);
			PKMER_EDGE path2Start = kmer_vertex_get_succ_edge(v, 1);
			long weight1 = 0;
			long weight2 = 0;

				if (path1Start->Type == kmetRead) {
					PKMER_EDGE tmp = path1Start;
					path1Start = path2Start;
					path2Start = tmp;
				}
				
				PKMER_VERTEX path1Vertex = path1Start->Dest;
				PKMER_VERTEX path2Vertex = path2Start->Dest;
				PKMER_EDGE e = path1Start;

				rs_storage_reset(&s1);
				rs_storage_reset(&s2);
				rs_storage_add_edge(&s1, path1Start, TRUE);
				weight1 = max(weight1, read_info_get_count(&path1Start->ReadInfo));
				rs_storage_add_edge(&s2, path2Start, TRUE);
				weight2 = max(weight2, read_info_get_count(&path2Start->ReadInfo));
				while (kmer_vertex_out_degree(path1Vertex) == 1 && kmer_vertex_in_degree(path1Vertex) == 1) {
					e = kmer_vertex_get_succ_edge(path1Vertex, 0);
					path1Vertex = e->Dest;
					rs_storage_add_edge(&s1, e, TRUE);
					weight1 = max(weight1, read_info_get_count(&e->ReadInfo));
				}

				while (kmer_vertex_out_degree(path2Vertex) == 1 && kmer_vertex_in_degree(path2Vertex) == 1) {
					e = kmer_vertex_get_succ_edge(path2Vertex, 0);
					path2Vertex = e->Dest;
					rs_storage_add_edge(&s2, e, TRUE);
					weight2 = max(weight2, read_info_get_count(&e->ReadInfo));
				}

				if (path1Vertex->Type == kmvtRefSeqMiddle && path1Vertex == path2Vertex) {
					rs_storage_remove(&s1, 1);
					rs_storage_remove(&s2, 1);
//					if (read_info_get_count(&path1Start->ReadInfo) > 0)  
					{
						EKMerEdgeType seq1Type = path1Start->Type;
						EKMerEdgeType seq2Type = path2Start->Type;
						
						kmer_graph_delete_edge(Graph, path1Start);
						kmer_graph_delete_edge(Graph, path2Start);
						ret = kmer_graph_add_edge_ex(Graph, v, path1Vertex, 1, min(s1.ValidLength, s2.ValidLength) + 1, kmetVariant, &e);
						if (ret == ERR_SUCCESS) {
							e->SeqLen = s1.ValidLength;
							ret = rs_storage_create_string(&s1, &e->Seq);
							if (ret == ERR_SUCCESS) {
								e->Seq2Len = s2.ValidLength;
								e->SeqType = seq1Type;
								e->Seq1Weight = weight1;
								e->Seq2Type = seq2Type;
								e->Seq2Weight = weight2;
								ret = rs_storage_create_string(&s2, &e->Seq2);
							}
						}
					}  
//					else kmer_graph_delete_edge(Graph, path1Start);

					*Changed = TRUE;
				}
		}

		if (ret != ERR_SUCCESS)
			break;

		ret = kmer_table_next(Graph->VertexTable, entry, &entry);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	size_t dummy = 0;
	kmer_graph_delete_trailing_things(Graph, &dummy);
	rs_storage_finit(&s2);
	rs_storage_finit(&s1);

	return ret;
}
