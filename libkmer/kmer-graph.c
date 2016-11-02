
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "err.h"
#include "utils.h"
#include "khash.h"
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

	ret = utils_malloc(sizeof(KMER_VERTEX) + kmer_get_size(KMer)*sizeof(char), &tmp);
	if (ret == ERR_SUCCESS) {
		kmer_init_from_kmer(&tmp->KMer, KMer);;
		tmp->Type = Type;
		pointer_array_init_KMER_EDGE(&tmp->Successors, 140);
		pointer_array_init_KMER_EDGE(&tmp->Predecessors, 140);
		tmp->RefSeqPosition = 0;
		tmp->Helper = FALSE;
		tmp->Lists.Next = NULL;
		tmp->Lists.Graph = NULL;
		*Result = tmp;
	}

	return ret;
}


static void _vertex_destroy(PKMER_VERTEX Vertex)
{
	PKMER_GRAPH g = Vertex->Lists.Graph;

	pointer_array_finit_KMER_EDGE(&Vertex->Predecessors);
	pointer_array_finit_KMER_EDGE(&Vertex->Successors);
	if (g != NULL) {
		Vertex->Lists.Next = g->VerticesToDeleteList;
		g->VerticesToDeleteList = Vertex;
	} else utils_free(Vertex);

	return;
}


static ERR_VALUE _vertex_copy(const KMER_VERTEX *Vertex, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(&Vertex->KMer, Vertex->Type, &tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Helper = Vertex->Helper;
		tmp->Order = Vertex->Order;
		tmp->Lists.Graph = Vertex->Lists.Graph;
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

static ERR_VALUE _edge_create(PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Source = Source;
		tmp->Dest = Dest;
		tmp->Type = Type;
		tmp->Order = 0;
		tmp->Seq = NULL;
		tmp->SeqLen = 0;
		tmp->Seq1Weight = 0;
		tmp->SeqType = kmetNone;
		tmp->Seq2 = NULL;
		tmp->Seq2Len = 0;
		tmp->Seq2Weight = 0;
		tmp->Seq2Type = kmetNone;
		read_info_init(&tmp->ReadInfo);
		tmp->MarkedForDelete = FALSE;
		dym_array_init_size_t(&tmp->Paths, 140);
		tmp->Finished = FALSE;
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

	ret = _edge_create(Edge->Source, Edge->Dest, Edge->Type, &tmp);
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
					tmp->Finished = Edge->Finished;
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
	char *colors[] = {"yellow", "lightgreen", "blue", "red", "white"};
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;


	fprintf(Stream, "\t");
	kmer_print(Stream, &v->KMer);
	fprintf(Stream, "[label=\"");
	kmer_print(Stream, &v->KMer);
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
	kmer_print(Stream, &e->Source->KMer);
	fprintf(Stream, " -> ");
	kmer_print(Stream, &e->Dest->KMer);
	fprintf(Stream, " [color=");
	switch (e->Type) {
		case kmetReference:
			fprintf(Stream, "green");
			fprintf(Stream, ",label=\"W: %lf (%Iu); L: %u; P: %u\";", read_info_weight(&e->ReadInfo, -1, -1), read_info_get_count(&e->ReadInfo), (uint32_t)e->SeqLen, (uint32_t)gen_array_size(&e->Paths));
			break;
		case kmetRead:
			fprintf(Stream, "red");
			fprintf(Stream, ",label=\"W: %lf (%Iu); L: %u; P:%u\"", read_info_weight(&e->ReadInfo, -1, -1), read_info_get_count(&e->ReadInfo),  (uint32_t)e->SeqLen, (uint32_t)gen_array_size(&e->Paths));
			break;
		case kmetVariant:
			fprintf(Stream, "blue");
			fprintf(Stream, ",label=\"W: %u; L: %u; P: %u\\n1: %s\\n2: %s\"", (uint32_t)e->Seq1Weight, (uint32_t)e->SeqLen, (uint32_t)gen_array_size(&e->Paths), e->Seq, e->Seq2);
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


static PKMER_EDGE _get_refseq_or_variant_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_succ_edge(Vertex, i);

		if (tmp->Type == kmetReference || tmp->Type == kmetVariant) {
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


/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, const size_t VerticesHint, const size_t EdgesHint, PKMER_GRAPH *Graph)
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
		tmpGraph->VerticesToDeleteList = NULL;
		vCallbacks.OnCopy = _vertex_table_on_copy;
		vCallbacks.OnDelete = _vertex_table_on_delete;
		vCallbacks.OnInsert = _vertex_table_on_insert;
		vCallbacks.OnPrint = _vertex_table_on_print;
		ret = kmer_table_create(KMerSize, 47/*utils_next_prime(VerticesHint)*/, &vCallbacks, &tmpGraph->VertexTable);
		if (ret == ERR_SUCCESS) {
			eCallbacks.OnCopy = _edge_table_on_copy;
			eCallbacks.OnDelete = _edge_table_on_delete;
			eCallbacks.OnInsert = _edge_table_on_insert;
			eCallbacks.OnPrint = _edge_table_on_print;
			ret = kmer_edge_table_create(KMerSize, 47, &eCallbacks, &tmpGraph->EdgeTable);
			if (ret == ERR_SUCCESS) {
				ret = kmer_edge_table_create(KMerSize, 47, NULL, &tmpGraph->DummyVertices);
				if (ret == ERR_SUCCESS)
					*Graph = tmpGraph;
				
				if (ret != ERR_SUCCESS)
					kmer_edge_table_destroy(tmpGraph->EdgeTable);
			}

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
	kmer_edge_table_destroy(Graph->DummyVertices);
	kmer_edge_table_destroy(Graph->EdgeTable);
	kmer_table_destroy(Graph->VertexTable);

	PKMER_VERTEX del = Graph->VerticesToDeleteList;
	PKMER_VERTEX old = Graph->VerticesToDeleteList;

	while (del != NULL) {
		old = del;
		del = del->Lists.Next;
		utils_free(old);
	}

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
	/*
	ret = ERR_SUCCESS;
	{
		REFSEQ_STORAGE rs;
		PKMER_VERTEX startVertex = NULL;
		PKMER_EDGE startingEdge = NULL;
		long w = 0;
		PKMER_EDGE edge = NULL;

		rs_storage_init(&rs);
		startVertex = _get_refseq_edge(Graph->StartingVertex)->Dest;
		startingEdge = _get_refseq_edge(startVertex);
		v = startVertex;
		while (startVertex != Graph->EndingVertex) {
			edge = _get_refseq_edge(v);
			v = edge->Dest;
			if (v == Graph->EndingVertex)
				break;

			if (kmer_vertex_in_degree(v) == 1 && kmer_vertex_out_degree(v) == 1) {
				rs_storage_add_edge(&rs, edge, FALSE);
				w = max(w, edge->Seq1Weight);
			} else {
				if (kmer_graph_get_edge(Graph, &startVertex->KMer, &v->KMer) == NULL) {
					rs_storage_add_edge(&rs, edge, FALSE);
					if (!v->Helper)
						rs_storage_remove(&rs, 1);
					
					kmer_graph_delete_edge(Graph, edge);
					if (startingEdge != edge)
						kmer_graph_delete_edge(Graph, startingEdge);
				
					ret = kmer_graph_add_edge_ex(Graph, startVertex, v, kmetReference, &edge);
					if (ret == ERR_SUCCESS) {
						edge->Seq1Weight = w;
						edge->SeqType = kmetReference;
						edge->SeqLen = rs.ValidLength;
						ret = rs_storage_create_string(&rs, &edge->Seq);
						rs_storage_reset(&rs);
					}
				}

				startVertex = v;
				startingEdge = _get_refseq_edge(startVertex);
			}
		}

		rs_storage_finit(&rs);
		
		size_t dummy = 0;

		kmer_graph_delete_trailing_things(Graph, &dummy);
	}
	*/
	
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


void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const size_t Threshold)
{
	PKMER_EDGE e = NULL;
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (err == ERR_SUCCESS) {
		e = (PKMER_EDGE)(iter->Data);
		err = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
		if (e->Type == kmetRead && read_info_get_count(&e->ReadInfo) == 0)
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


ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	v = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	if (v == NULL) {
		ret = _vertex_create(KMer, Type, &v);
		if (ret == ERR_SUCCESS) {
			do {
				ret = kmer_table_insert(Graph->VertexTable, &v->KMer, v);
				if (ret == ERR_SUCCESS) {
					Graph->NumberOfVertices++;
					v->Lists.Graph = Graph;
					*Vertex = v;
				} else if (ret == ERR_TABLE_FULL) {
					ret = kmer_table_extend(Graph->VertexTable);
					if (ret == ERR_SUCCESS)
						ret = ERR_TABLE_FULL;
				}
			} while (ret == ERR_TABLE_FULL);

			if (ret != ERR_SUCCESS)
				_vertex_destroy(v);
		}
	} else {
		*Vertex = v;
		ret = ERR_ALREADY_EXISTS;
	}

	return ret;
}


ERR_VALUE kmer_graph_add_helper_vertex(PKMER_GRAPH Graph, const KMER *KMer1, const KMER *KMer2, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *dummySeq = alloca(kmerSize*sizeof(char));
	PKMER kmer = NULL;

	ret = ERR_SUCCESS;
	memset(dummySeq, 'D', kmerSize*sizeof(char));
	KMER_STACK_ALLOC(kmer, 0, kmerSize, dummySeq);
	v = (PKMER_VERTEX)kmer_edge_table_get(Graph->DummyVertices, KMer1, KMer2);
	if (v == NULL) {
		do {
			ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
			if (ret == ERR_ALREADY_EXISTS)
				kmer_get_number(kmer) += 1;
		} while (ret == ERR_ALREADY_EXISTS);
	
		if (ret == ERR_SUCCESS) {
			v->Helper = TRUE;
			do {
				ret = kmer_edge_table_insert(Graph->DummyVertices, KMer1, KMer2, v);
				if (ret == ERR_TABLE_FULL) {
					ret = kmer_edge_table_extend(Graph->DummyVertices);
					if (ret == ERR_SUCCESS)
						ret = ERR_TABLE_FULL;
				}
			} while (ret == ERR_TABLE_FULL);

			if (ret == ERR_SUCCESS)
				*Vertex = v;

			if (ret != ERR_SUCCESS)
				kmer_graph_delete_vertex(Graph, v);
		}
	} else {
		ret = ERR_ALREADY_EXISTS;
		*Vertex = v;
	}


	return ret;
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	edge = (PKMER_EDGE)kmer_edge_table_get(Graph->EdgeTable, &Source->KMer, &Dest->KMer);
	if (edge == NULL) {
		ret = _edge_create(Source, Dest, Type, &edge);
		if (ret == ERR_SUCCESS) {
			do {
				ret = kmer_edge_table_insert(Graph->EdgeTable, &Source->KMer, &Dest->KMer, edge);
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
		}
	} else {
		*Edge = edge;
		ret = ERR_ALREADY_EXISTS;
	}

	return ret;
}


PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest)
{
	return (PKMER_EDGE)kmer_edge_table_get(Graph->EdgeTable, Source, Dest);
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
	ret = dym_array_reserve(&tmp, 20);
	if (ret == ERR_SUCCESS) {
		ret = kmer_table_get_multiple(Graph->VertexTable, KMer, &tmp);
		if (ret == ERR_SUCCESS) {
			ret = pointer_array_reserve_KMER_VERTEX(VertexArray, dym_array_size(&tmp));
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < dym_array_size(&tmp); ++i) {
					pointer_array_push_back_no_alloc_KMER_VERTEX(VertexArray, (PKMER_VERTEX)dym_array_get(&tmp, i));
				}
			}
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
				ret = kmer_table_delete(Graph->VertexTable, &Vertex->KMer);
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

	err = kmer_edge_table_delete(Graph->EdgeTable, &source->KMer, &dest->KMer);
	if (err == ERR_SUCCESS) {
		pointer_array_remove_by_item_fast_KMER_EDGE(&source->Successors, Edge);
		pointer_array_remove_by_item_fast_KMER_EDGE(&dest->Predecessors, Edge);
		--Graph->TypedEdgeCount[edgeType];
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
		if (kmer_graph_get_edge(Graph, &u->KMer, &w->KMer) == NULL) {
			boolean mfd = (Source->MarkedForDelete || Dest->MarkedForDelete);
			size_t pathCount = max(gen_array_size(&Source->Paths), gen_array_size(&Dest->Paths));
			EKMerEdgeType type = (Source->Type == kmetReference && Dest->Type == kmetReference) ? kmetReference : kmetRead;
			PKMER_EDGE newEdge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, u, w, type, &newEdge);
			if (ret == ERR_SUCCESS) {
				newEdge->MarkedForDelete = mfd;
				newEdge->SeqLen = Source->SeqLen + 1 + Dest->SeqLen;
				ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
				if (ret == ERR_SUCCESS) {
					memcpy(newEdge->Seq, Source->Seq, Source->SeqLen*sizeof(char));
					if (!v->Helper) {
						newEdge->Seq[Source->SeqLen] = kmer_get_base((&v->KMer), Graph->KMerSize - 1);
						memcpy(newEdge->Seq + Source->SeqLen + 1, Dest->Seq, Dest->SeqLen*sizeof(char));
					} else {
						memcpy(newEdge->Seq + Source->SeqLen, Dest->Seq, Dest->SeqLen*sizeof(char));
						--newEdge->SeqLen;
					}

					newEdge->Seq[newEdge->SeqLen] = '\0';
					newEdge->SeqType = Source->SeqType;
					newEdge->Seq1Weight = max(Source->Seq1Weight, Dest->Seq1Weight);
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


ERR_VALUE kmer_graph_get_splitted_edge(PKMER_GRAPH Graph, const KMER_VERTEX *Source, const KMER_VERTEX *Dest, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX helperVertex = NULL;
	PKMER_EDGE es = NULL;
	PKMER_EDGE ed = NULL;

	helperVertex = (PKMER_VERTEX)kmer_edge_table_get(Graph->DummyVertices, &Source->KMer, &Dest->KMer);
	if (helperVertex != NULL) {
		es = kmer_graph_get_edge(Graph, &Source->KMer, &helperVertex->KMer);
		ed = kmer_graph_get_edge(Graph, &helperVertex->KMer, &Dest->KMer);
		ret = ERR_SUCCESS;
		if (SplitVertex != NULL)
			*SplitVertex = helperVertex;

		if (SourceEdge != NULL)
			*SourceEdge = es;

		if (DestEdge != NULL)
			*DestEdge = ed;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_graph_split_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex)
{
	PKMER_EDGE es = NULL;
	PKMER_EDGE ed = NULL;
	PKMER_VERTEX helperVertex = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_add_helper_vertex(Graph, &Edge->Source->KMer, &Edge->Dest->KMer, &helperVertex);
	if (ret == ERR_SUCCESS) {
		ret = kmer_graph_add_edge_ex(Graph, Edge->Source, helperVertex, Edge->Type, &es);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_add_edge_ex(Graph, helperVertex, Edge->Dest, Edge->Type, &ed);
			if (ret == ERR_SUCCESS) {
				char *seq = NULL;
				size_t seqLen = Edge->SeqLen;

				if (seqLen > 0) {
					ret = utils_calloc(seqLen + 1, sizeof(char), &seq);
					if (ret == ERR_SUCCESS) {
						memcpy(seq, Edge->Seq, seqLen*sizeof(char));
						seq[seqLen] = '\0';
						kmer_edge_add_seq(es, Edge->SeqType, seq, seqLen);
						seq = NULL;
					}
				} else es->Seq1Weight = 1;

				es->Seq1Weight = Edge->Seq1Weight;
				ed->Seq1Weight = Edge->Seq1Weight;
				if (ret == ERR_SUCCESS) {
					ret = read_info_copy(&es->ReadInfo, &Edge->ReadInfo);
					if (ret == ERR_SUCCESS) {
						ret = read_info_copy(&ed->ReadInfo, &Edge->ReadInfo);
						if (ret == ERR_SUCCESS) {
							kmer_graph_delete_edge(Graph, Edge);
							if (SourceEdge != NULL)
								*SourceEdge = es;

							if (DestEdge != NULL)
								*DestEdge = ed;

							if (SplitVertex != NULL)
								*SplitVertex = helperVertex;
						}
					}
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, ed);
			}

			if (ret != ERR_SUCCESS)
				kmer_graph_delete_edge(Graph, es);
		}

		if (ret != ERR_SUCCESS)
			kmer_graph_delete_vertex(Graph, helperVertex);
	} else if (ret == ERR_ALREADY_EXISTS) {
		if (SplitVertex != NULL)
			*SplitVertex = helperVertex;

		if (SourceEdge != NULL)
			*SourceEdge = kmer_graph_get_edge(Graph, &Edge->Source->KMer, &helperVertex->KMer);

		if (DestEdge != NULL)
			*DestEdge = kmer_graph_get_edge(Graph, &helperVertex->KMer, &Edge->Dest->KMer);
	}

	assert(*SourceEdge != NULL);
	assert(*DestEdge != NULL);
	assert(*SplitVertex != NULL);

	return ret;
}


ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, const char *RefSeq, PPOINTER_ARRAY_FOUND_SEQUENCE SeqArray)
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
	rs_storage_add_seq(&rsStorage, RefSeq, kmerSize - 1);
	dym_array_init_FOUND_SEQUENCE_VARIANT(&variantStack, 140);
	dym_array_init_size_t(&edgeIndices, 140);
	pointer_array_init_KMER_EDGE(&edges, 140);
		do {
			if (edgeIndex == kmer_vertex_out_degree(v)) {
				if (v == Graph->StartingVertex)
					break;

				edge = *pointer_array_pop_back_KMER_EDGE(&edges);
				edge->Finished = FALSE;
				edgeIndex = *dym_array_pop_back_size_t(&edgeIndices);
				v = edge->Source;
				rs_storage_remove_edge(&rsStorage, edge);				
				if (edge->Type == kmetVariant)
					dym_array_pop_back_FOUND_SEQUENCE_VARIANT(&variantStack);
			} else {
				edge = kmer_vertex_get_succ_edge(v, edgeIndex);
				if (!edge->Finished) {
					if (edge->Type == kmetVariant) {
						FOUND_SEQUENCE_VARIANT variant;

						ret = rs_storage_add(&rsStorage, '?');
						if (ret != ERR_SUCCESS)
							break;

						assert(edge->Dest->Type == kmvtRefSeqMiddle);
						ret = rs_storage_add_vertex(&rsStorage, edge->Dest);
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
						ret = rs_storage_add_edge(&rsStorage, edge, FALSE);
						if (ret != ERR_SUCCESS)
							break;
					}

					edge->Finished = TRUE;
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

							ret = found_sequence_build_read_variants(fs, &edges);
							if (ret == ERR_SUCCESS) {
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
											edge->Finished = FALSE;
											rs_storage_remove_edge(&rsStorage, edge);
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
							}

							if (ret != ERR_SUCCESS)
								found_sequence_free(fs);
						}
					} else {
						v = edge->Dest;
						ret = dym_array_push_back_size_t(&edgeIndices, (edgeIndex + 1));
						if (ret != ERR_SUCCESS)
							break;

						ret = pointer_array_push_back_KMER_EDGE(&edges, edge);
						if (ret != ERR_SUCCESS)
							break;

						edgeIndex = 0;
					}
				} else ++edgeIndex;
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


void kmer_edge_add_seq(PKMER_EDGE Edge, EKMerEdgeType Type, const char *Seq, const size_t Length)
{
	if (Edge->SeqType == kmetNone) {
		Edge->SeqType = Type;
		Edge->Seq = Seq;
		Edge->SeqLen = Length;
		Edge->Seq1Weight = 1;
	} else if (Edge->Type == kmetVariant) {
		if (Edge->Seq2Type == kmetNone) {
			Edge->Seq2Type = Type;
			Edge->Seq2 = Seq;
			Edge->Seq2Len = Length;
			Edge->Seq2Weight = 1;
		} else {
			printf("Attempt to add second sequence to full variant edge\n");
			exit(0);
		}
	} else {
		printf("Attempt to add second sequence to a non-variant edge\n");
		exit(0);
	}

	return;
}


static ERR_VALUE _capture_refseq(const KMER_EDGE *Start, const POINTER_ARRAY_KMER_EDGE *RSEdges, const KMER_EDGE *End, char **Seq, size_t *SeqLen)
{
	REFSEQ_STORAGE rsStorage;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t kmerSize = kmer_get_size(&Start->Dest->KMer);

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
				if (!End->Dest->Helper &&
					End->Dest->Type != kmvtRefSeqEnd)
					rs_storage_remove(&rsStorage, 1);
				
				ret = rs_storage_create_string(&rsStorage, Seq);
				*SeqLen = rsStorage.ValidLength;
			}
		}
	}

	rs_storage_finit(&rsStorage);

	return ret;
}


static void _remove_read_info_from_edge_array(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY ToRemove, size_t Distance, size_t *NewDistance)
{
	for (size_t k = 0; k < pointer_array_size(Array); ++k) {
		PKMER_EDGE tmp = *pointer_array_item_KMER_EDGE(Array, k);

		read_info_subtract(&tmp->ReadInfo, ToRemove, Distance);
		Distance += tmp->SeqLen;
		if (!tmp->Dest->Helper)
			Distance++;
	}

	if (NewDistance != NULL)
		*NewDistance = Distance;

	return;
}


static void _remove_read_info_from_edges(PKMER_GRAPH Graph, PKMER_EDGE Edge1, PKMER_EDGE Edge2, PPOINTER_ARRAY_KMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY Intersection)
{
	size_t distance = 0;

	read_info_subtract(&Edge1->ReadInfo, Intersection, distance);
	distance += Edge1->SeqLen;
	if (!Edge1->Dest->Helper)
		++distance;

	_remove_read_info_from_edge_array(Graph, Array, Intersection, distance, &distance);
	read_info_subtract(&Edge2->ReadInfo, Intersection, distance);

	return;
}


typedef struct _EDGE_REMOVE_CONTEXT {
	POINTER_ARRAY_KMER_EDGE RSEdges;
	PKMER_EDGE TargetEdge;
	size_t ReadDistance;
} EDGE_REMOVE_CONTEXT, *PEDGE_REMOVE_CONTEXT;

POINTER_ARRAY_TYPEDEF(EDGE_REMOVE_CONTEXT);
POINTER_ARRAY_IMPLEMENTATION(EDGE_REMOVE_CONTEXT)


static ERR_VALUE _remove_context_create(PPOINTER_ARRAY_KMER_EDGE RSEdges, PKMER_EDGE TargetEdge, const size_t ReadDistance, PEDGE_REMOVE_CONTEXT *Context)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PEDGE_REMOVE_CONTEXT tmpContext = NULL;

	ret = utils_malloc(sizeof(EDGE_REMOVE_CONTEXT), &tmpContext);
	if (ret == ERR_SUCCESS) {
		tmpContext->ReadDistance = ReadDistance;
		tmpContext->TargetEdge = TargetEdge;
		pointer_array_init_KMER_EDGE(&tmpContext->RSEdges, 140);
		ret = pointer_array_push_back_array_KMER_EDGE(&tmpContext->RSEdges, RSEdges);
		if (ret == ERR_SUCCESS)
			*Context = tmpContext;
		
		if (ret != ERR_SUCCESS) {
			pointer_array_finit_KMER_EDGE(&tmpContext->RSEdges);
			utils_free(tmpContext);
		}
	}

	return ret;
}


static void _remove_context_destroy(PEDGE_REMOVE_CONTEXT Context)
{
	pointer_array_finit_KMER_EDGE(&Context->RSEdges);
	utils_free(Context);

	return;
}

static void _remove_context_apply(PKMER_GRAPH Graph, PEDGE_REMOVE_CONTEXT Context)
{
	_remove_read_info_from_edge_array(Graph, &Context->RSEdges, &Context->TargetEdge->ReadInfo.Array, Context->ReadDistance, NULL);
	_remove_context_destroy(Context);

	return;
}

ERR_VALUE kmer_graph_connect_reads_by_pairs(PKMER_GRAPH Graph, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, size_t *ChangeCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_READ_INFO_ENTRY intersection;
	PKMER_EDGE rsLastEdge = NULL;
	PKMER_EDGE rsNextEdge = NULL;
	const KMER_EDGE *nextRsLastEdge = NULL;
	POINTER_ARRAY_KMER_EDGE edgesToDelete;
	POINTER_ARRAY_KMER_EDGE rsEdges;
	size_t dummy = 0;
	boolean deleteEOut = FALSE;
	boolean deleteeIn = FALSE;
	POINTER_ARRAY_TYPE(EDGE_REMOVE_CONTEXT) removeContexts;

	pointer_array_init_EDGE_REMOVE_CONTEXT(&removeContexts, 140);
	ret = pointer_array_reserve_EDGE_REMOVE_CONTEXT(&removeContexts, 2 * gen_array_size(PairArray));
	if (ret == ERR_SUCCESS) {
		dym_array_init_READ_INFO_ENTRY(&intersection, 140);
		pointer_array_init_KMER_EDGE(&edgesToDelete, 140);
		ret = pointer_array_reserve_KMER_EDGE(&edgesToDelete, 2 * gen_array_size(PairArray));
		if (ret == ERR_SUCCESS) {
			for (size_t h = 0; h < gen_array_size(PairArray); ++h) {
				boolean edgeCreated = FALSE;
				KMER_EDGE_PAIR pair = *dym_array_item_KMER_EDGE_PAIR(PairArray, h);
				PKMER_EDGE eIn = pair.U;
				PKMER_EDGE eOut = pair.V;

				if (eIn == eOut)
					continue;

				rsLastEdge = _get_in_refseq_edge(eIn->Dest);
				rsNextEdge = _get_refseq_edge(eOut->Source);
				rsEdges.Data = pair.Edges;
				rsEdges.ValidLength = pair.EdgeCount;
				rsEdges.AllocLength = pair.EdgeCount;
				if (pointer_array_size(&rsEdges) != pair.ReadDistance)
					continue;

					ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, eIn->SeqLen + (!eIn->Dest->Helper ? 1 : 0) + pair.ReadDistance);
					if (ret == ERR_SUCCESS && gen_array_size(&intersection) > Threshold) {
						if (kmer_equal(&eIn->Source->KMer, &eOut->Dest->KMer))
							ret = ERR_ALREADY_EXISTS;

						if (ret == ERR_SUCCESS) {
							*ChangeCount++;
							edgeCreated = TRUE;
							_remove_read_info_from_edges(Graph, eIn, eOut, &rsEdges, &intersection);
						}
						
						ret = ERR_SUCCESS;
					} else if (ret == ERR_SUCCESS) {
						PKMER_EDGE e = kmer_graph_get_edge(Graph, &eIn->Source->KMer, &eOut->Dest->KMer);

						if (e != NULL)
							kmer_graph_delete_edge(Graph, e);
					}

					if (edgeCreated) {
						deleteeIn = FALSE;
						deleteEOut = FALSE;
						if (rsNextEdge != NULL) {
							ret = read_info_intersection(&eIn->ReadInfo, &rsNextEdge->ReadInfo, &intersection, eIn->SeqLen + 1 + pair.ReadDistance);
							if (ret == ERR_SUCCESS)
								deleteeIn = (gen_array_size(&intersection) <= Threshold);
						} else deleteeIn = TRUE;

						if (ret == ERR_SUCCESS) {
							if (rsLastEdge != NULL) {
								ret = read_info_intersection(&rsLastEdge->ReadInfo, &eOut->ReadInfo, &intersection, rsLastEdge->SeqLen + 1 + pair.ReadDistance);
								if (ret == ERR_SUCCESS)
									deleteEOut = (gen_array_size(&intersection) <= Threshold);
							} else deleteEOut = TRUE;
						}

						if (deleteeIn) {
							PEDGE_REMOVE_CONTEXT removeContext = NULL;

							ret = _remove_context_create(&rsEdges, eIn, eIn->SeqLen + 1, &removeContext);
							if (ret == ERR_SUCCESS) {
								pointer_array_push_back_no_alloc_EDGE_REMOVE_CONTEXT(&removeContexts, removeContext);
								if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eIn))
									pointer_array_push_back_no_alloc_KMER_EDGE(&edgesToDelete, eIn);
							}
						}

						if (deleteEOut) {
							PEDGE_REMOVE_CONTEXT removeContext = NULL;

							ret = _remove_context_create(&rsEdges, eOut, -pair.ReadDistance, &removeContext);
							if (ret == ERR_SUCCESS) {
								pointer_array_push_back_no_alloc_EDGE_REMOVE_CONTEXT(&removeContexts, removeContext);
								if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eOut))
									pointer_array_push_back_no_alloc_KMER_EDGE(&edgesToDelete, eOut);
							}
						}
					}

				if (ret != ERR_SUCCESS)
					break;
			}

			dym_array_finit_READ_INFO_ENTRY(&intersection);

			if (ret == ERR_NO_MORE_ENTRIES)
				ret = ERR_SUCCESS;

			PEDGE_REMOVE_CONTEXT *pctx = removeContexts.Data;
			for (size_t i = 0; i < pointer_array_size(&removeContexts); ++i) {
				_remove_context_apply(Graph, *pctx);
				++pctx;
			}

			for (size_t i = 0; i < pointer_array_size(&edgesToDelete); ++i) {
				kmer_graph_delete_edge(Graph, *pointer_array_item_KMER_EDGE(&edgesToDelete, i));
				*ChangeCount++;
			}
		}

		pointer_array_finit_KMER_EDGE(&edgesToDelete);
	}

	pointer_array_finit_EDGE_REMOVE_CONTEXT(&removeContexts);

	kmer_graph_delete_trailing_things(Graph, &dummy);

	return ret;
}


KHASH_SET_INIT_INT(es);

ERR_VALUE kmer_graph_detect_uncertainities(PKMER_GRAPH Graph, boolean *Changed)
{
	boolean edgeCreated = FALSE;
	REFSEQ_STORAGE s1;
	REFSEQ_STORAGE s2;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE e = NULL;
	PKMER_VERTEX v = NULL;

	ret = ERR_SUCCESS;
	rs_storage_init(&s1);
	rs_storage_init(&s2);
	v = _get_refseq_edge(Graph->StartingVertex)->Dest;
	while (v != Graph->EndingVertex) {
		edgeCreated = FALSE;
		if (v->Type == kmvtRefSeqMiddle && kmer_vertex_out_degree(v) == 2) {
			PKMER_EDGE path1Start = kmer_vertex_get_succ_edge(v, 0);
			PKMER_EDGE path2Start = kmer_vertex_get_succ_edge(v, 1);
			size_t weight1 = 0;
			size_t weight2 = 0;

			if (path1Start->Type == kmetRead) {
				PKMER_EDGE tmp = path1Start;
				path1Start = path2Start;
				path2Start = tmp;
			}
				
			EKMerEdgeType seq1Type = path1Start->Type;
			EKMerEdgeType seq2Type = path2Start->Type;
			PKMER_VERTEX path1Vertex = path1Start->Dest;
			PKMER_VERTEX path2Vertex = path2Start->Dest;

			{
				khash_t(es)	*table = kh_init(es);

				rs_storage_reset(&s2);
				rs_storage_add_edge(&s2, path2Start, TRUE);
				weight2 = max(weight2, path2Start->Seq1Weight);
				while (kmer_vertex_in_degree(path2Vertex) == 1 && kmer_vertex_out_degree(path2Vertex) == 1 && path2Vertex->Type == kmvtRead) {
					int r;
					PKMER_EDGE e = NULL;

					e = kmer_vertex_get_succ_edge(path2Vertex, 0);
					if (kh_get(es, table, e->Order) != kh_end(table))
						break;

					path2Vertex = e->Dest;
					rs_storage_add_edge(&s2, e, TRUE);
					weight2 = max(weight2, e->Seq1Weight);
					kh_put(es, table, e->Order, &r);
				}

				kh_destroy(es, table);
			}

			if (path2Vertex->Type == kmvtRefSeqMiddle) {
				rs_storage_reset(&s1);
				rs_storage_add_edge(&s1, path1Start, TRUE);
				weight1 = max(weight1, path1Start->Seq1Weight);
				while (ret == ERR_SUCCESS && kmer_vertex_in_degree(path1Vertex) == 1 && kmer_vertex_out_degree(path1Vertex) == 1 && path1Vertex->Type == kmvtRefSeqMiddle && path1Vertex != path2Vertex) {
					PKMER_EDGE e = NULL;

					e = kmer_vertex_get_succ_edge(path1Vertex, 0);					
					path1Vertex = e->Dest;
					ret = rs_storage_add_edge(&s1, e, TRUE);
					if (e->Type == kmetVariant)
						seq1Type = kmetRead;

					weight1 = max(weight1, e->Seq1Weight);
				}

				if (ret == ERR_SUCCESS && path1Vertex == path2Vertex) {
					PKMER_EDGE e = NULL;

					if (!path1Vertex->Helper)
						rs_storage_remove(&s1, 1);
					
					if (!path2Vertex->Helper)
						rs_storage_remove(&s2, 1);
					
					kmer_graph_delete_edge(Graph, path1Start);
					kmer_graph_delete_edge(Graph, path2Start);
					ret = kmer_graph_add_edge_ex(Graph, v, path1Vertex, kmetVariant, &e);
					if (ret == ERR_SUCCESS) {
						char *tmpSeq = NULL;
						
						v = _get_refseq_or_variant_edge(path1Vertex)->Dest;
						edgeCreated = TRUE;
						ret = rs_storage_create_string(&s1, &tmpSeq);
						if (ret == ERR_SUCCESS) {
							kmer_edge_add_seq(e, seq1Type, tmpSeq, s1.ValidLength);
							e->Seq1Weight = weight1;
							ret = rs_storage_create_string(&s2, &tmpSeq);
							if (ret == ERR_SUCCESS) {
								kmer_edge_add_seq(e, seq2Type, tmpSeq, s2.ValidLength);
								e->Seq2Weight = weight2;
							}
						}
					}

					*Changed = TRUE;
				}

				if (ret == ERR_TWO_READ_SEQUENCES)
					ret = ERR_SUCCESS;
			}
		}

		if (ret != ERR_SUCCESS)
			break;

		if (!edgeCreated)
			v = _get_refseq_or_variant_edge(v)->Dest;
	}

	size_t dummy = 0;
	kmer_graph_delete_trailing_things(Graph, &dummy);
	rs_storage_finit(&s2);
	rs_storage_finit(&s1);

	return ret;
}
