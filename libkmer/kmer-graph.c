
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
#include "kmer-graph.h"


/************************************************************************/
/*                      HELPER FUNCTIONS                                */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_GRAPH)
static UTILS_TYPED_MALLOC_FUNCTION(KMER_VERTEX)

/************************************************************************/
/*                        VERTEX BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _vertex_create(const PKMER KMer, EKMerVertexType Type, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_VERTEX(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->KMer = kmer_copy(KMer);
		ret = (tmp->KMer != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		if (ret == ERR_SUCCESS) {
			tmp->DegreeIn = 0;
			tmp->degreeOut = 0;
			tmp->Type = Type;
			dym_array_create(&tmp->Successors, 140);
			dym_array_create(&tmp->Predecessors, 140);
			*Result = tmp;
		}
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
	fprintf(Stream, "\\nIN=%u; OUT=%u\",style=filled,color=%s]", v->DegreeIn, v->degreeOut, colors[v->Type]);
	fprintf(Stream, ";\n");

	return;
}


/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	KMER_TABLE_CALLBACKS vCallbacks;

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
			ret = kmer_edge_table_create(KMerSize, 2, 37, &tmpGraph->EdgeTable);
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


ERR_VALUE kmer_graph_copy(const PKMER_GRAPH Source, PKMER_GRAPH *Copied)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_GRAPH tmpGraph = NULL;

	ret = utils_malloc_KMER_GRAPH(&tmpGraph);
	if (ret == ERR_SUCCESS) {
		tmpGraph->KMerSize = Source->KMerSize;
		tmpGraph->NumberOfEdges = Source->NumberOfEdges;
		tmpGraph->NumberOfVertices = Source->NumberOfVertices;
		// TODO: Copy also the starting and ending vertices
		tmpGraph->StartingVertex = NULL;
		tmpGraph->EndingVertex = NULL;
		ret = kmer_table_copy(Source->VertexTable, &tmpGraph->VertexTable);
		if (ret == ERR_SUCCESS) {
			ret = kmer_edge_table_copy(Source->EdgeTable, &tmpGraph->EdgeTable);
			if (ret == ERR_SUCCESS)
				*Copied = tmpGraph;

			if (ret != ERR_SUCCESS)
				kmer_table_destroy(tmpGraph->VertexTable);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpGraph);
	}

	return ret;
}


void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph)
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
	PKMER x = NULL;
	PKMER y = NULL;
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
		if (v->DegreeIn == 1 && v->degreeOut == 1) {
			x = (PKMER)dym_array_data(&v->Predecessors)[0];
			assert(x != NULL);
			y = (PKMER)dym_array_data(&v->Successors)[0];
			assert(y != NULL);
			if (!kmer_equal(x, y) && kmer_edge_table_get(Graph->EdgeTable, x, y) == NULL) {
				PKMER oldV = v->KMer;
				double prob = 0;
				uint32_t maxPassCount = 0;

				sourceEdge = kmer_edge_table_get(Graph->EdgeTable, x, v->KMer);
				assert(sourceEdge != NULL);
				destEdge = kmer_edge_table_get(Graph->EdgeTable, v->KMer, y);
				assert(destEdge != NULL);
				sourceWeight = sourceEdge->Weight;
				destWeight = destEdge->Weight;
				if (sourceWeight == 0)
					sourceWeight = destWeight;

				if (destWeight == 0)
					destWeight = sourceWeight;

				maxPassCount = min(sourceEdge->MaxPassCount, destEdge->MaxPassCount);
				weight = (sourceWeight + destWeight) / 2;
				prob = min(sourceEdge->Probability, destEdge->Probability);
				length = sourceEdge->Length + destEdge->Length;
				kmer_edge_table_delete(Graph->EdgeTable, sourceEdge->Source, sourceEdge->Dest);
				kmer_edge_table_delete(Graph->EdgeTable, destEdge->Source, destEdge->Dest);
				kmer_table_delete(Graph->VertexTable, v->KMer);
				ret = kmer_edge_table_insert_ex(Graph->EdgeTable, x, y, &dummy);
				if (ret == ERR_SUCCESS) {
					dummy->Length = length;
					dummy->Weight = weight;
					dummy->Probability = prob;
					dummy->MaxPassCount = maxPassCount;
					dummy->PassCount = 0;
					v = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, x);
					assert(v != NULL);
					dym_array_replace(&v->Successors, oldV, y);
					v = kmer_table_get(Graph->VertexTable, y);
					assert(v != NULL);
					dym_array_replace(&v->Predecessors, oldV, x);
				}

				--Graph->NumberOfVertices;
				--Graph->NumberOfEdges;
			}
		}

		if (last)
			break;
	}

	return ret;
}


void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold)
{
	PKMER_EDGE e = NULL;
	PKMER_EDGE iter = NULL;
	boolean last = FALSE;
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX v = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (err == ERR_SUCCESS) {
		boolean last = FALSE;

		e = iter;
		last = (kmer_edge_table_next(Graph->EdgeTable, e, &iter) == ERR_NO_MORE_ENTRIES);
		if (e->Weight < Threshold) {
			// TODO: Delete the edge
			u = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, e->Source);
			assert(u != NULL);
			v = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, e->Dest);
			assert(v != NULL);
			if (u->Type != kmvtRefSeqStart && u->Type != kmvtRefSeqEnd &&
				v->Type != kmvtRefSeqStart && v->Type != kmvtRefSeqEnd &&
				(e->Type != kmetReference)
				) {
				kmer_edge_table_delete_by_entry(Graph->EdgeTable, e);
				dym_array_remove_by_item_fast(&u->Successors, v->KMer);
				u->degreeOut--;
				dym_array_remove_by_item_fast(&v->Predecessors, u->KMer);
				v->DegreeIn--;
				Graph->NumberOfEdges--;
				if (u->DegreeIn == 0 && u->degreeOut == 0) {
					kmer_table_delete(Graph->VertexTable, u->KMer);
					Graph->NumberOfVertices--;
				}

				if (v->DegreeIn == 0 && v->degreeOut == 0) {
					kmer_table_delete(Graph->VertexTable, v->KMer);
					Graph->NumberOfVertices--;
				}
			}
		}

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


void kmer_graph_compute_edge_probablities(PKMER_GRAPH Graph)
{
	long total = 0;
	PKMER_VERTEX v = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret != ERR_NO_MORE_ENTRIES) {
		total = 0;
		v = (PKMER_VERTEX)iter->Data;
		for (size_t i = 0; i < v->degreeOut; ++i) {
			PKMER_EDGE e = kmer_edge_table_get(Graph->EdgeTable, v->KMer, (PKMER)dym_array_data(&v->Successors)[i]);
		
			total += (e->Weight + 1);
		}

		for (size_t i = 0; i < v->degreeOut; ++i) {
			PKMER_EDGE e = kmer_edge_table_get(Graph->EdgeTable, v->KMer, (PKMER)dym_array_data(&v->Successors)[i]);

			e->Probability = (total != 0) ? ((double)(e->Weight + 1) / (double)total) : (1/(double)v->degreeOut);
			if (e->Probability != 0)
				e->Probability = log(e->Probability);
		}

		ret = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}


ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer, const EKMerVertexType Type)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(KMer, Type, &v);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_table_insert(Graph->VertexTable, KMer, v);
			if (ret == ERR_TABLE_FULL)
				kmer_table_extend(Graph->VertexTable);
		} while (ret == ERR_TABLE_FULL);

		if (ret == ERR_SUCCESS)
			Graph->NumberOfVertices++;

		if (ret != ERR_SUCCESS)
			_vertex_destroy(v);
	}

	return ret;
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, PKMER_EDGE *Edge)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	do {
		ret = kmer_edge_table_insert_ex(Graph->EdgeTable, Source, Dest, &edge);
		if (ret == ERR_SUCCESS) {
			PKMER_VERTEX u = NULL;
			PKMER_VERTEX v = NULL;

			edge->Weight = weight;
			edge->Length = Length;
			edge->PassCount = 0;
			edge->MaxPassCount = 0;
			u = kmer_table_get(Graph->VertexTable, edge->Source);
			v = kmer_table_get(Graph->VertexTable, edge->Dest);
			assert(u != NULL);
			assert(v != NULL);
			ret = dym_array_prepare_for_insert(&u->Successors, 1);
			if (ret == ERR_SUCCESS) {
				ret = dym_array_prepare_for_insert(&v->Predecessors, 1);
				if (ret == ERR_SUCCESS) {
					u->degreeOut++;
					dym_array_push_back_no_alloc(&u->Successors, v->KMer);
					v->DegreeIn++;
					dym_array_push_back_no_alloc(&v->Predecessors, u->KMer);
					if (v->Order <= u->Order)
						Graph->NumberOfBackwardEdges++;

					*Edge = edge;
				}
			}
		}

		if (ret == ERR_TABLE_FULL)
			kmer_edge_table_extend(Graph->EdgeTable);
	} while (ret == ERR_TABLE_FULL);

	if (ret == ERR_ALREADY_EXISTS)
		*Edge = edge;

	if (ret == ERR_SUCCESS)
		Graph->NumberOfEdges++;

	return ret;
}


ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight)
{
	PKMER_EDGE dummy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_add_edge_ex(Graph, Source, Dest, weight, 1, &dummy);

	return ret;
}


PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest)
{
	return kmer_edge_table_get(Graph->EdgeTable, Source, Dest);
}
