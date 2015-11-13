
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "kmer-graph.h"


/************************************************************************/
/*                      HELPER FUNCTIONS                                */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_GRAPH)

/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_GRAPH(&tmpGraph);
	if (ret == ERR_SUCCESS) {
		tmpGraph->NumberOfEdges = 0;
		tmpGraph->NumberOfVertices = 0;
		tmpGraph->KMerSize = KMerSize;
		tmpGraph->NumberOfBackwardEdges = 0;
		ret = kmer_table_create(KMerSize, 2, 37, &tmpGraph->VertexTable);
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
	PKMER_TABLE_ENTRY v = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	PKMER x = NULL;
	PKMER y = NULL;
	PKMER_EDGE dummy = NULL;
	PKMER_EDGE sourceEdge = NULL;
	PKMER_EDGE destEdge = NULL;
	long weight = 0;
	uint32_t length = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		boolean last = FALSE;
		
		v = iter;
		last = (kmer_table_next(Graph->VertexTable, v, &iter) == ERR_NO_MORE_ENTRIES);
		if (v->DegreeIn == 1 && v->degreeOut == 1) {
			x = v->AdvancedInfo.PassThroughVertex.Input;
			assert(x != NULL);
			y = v->AdvancedInfo.PassThroughVertex.Output;
			assert(y != NULL);
			if (!kmer_equal(x, y) && kmer_edge_table_get(Graph->EdgeTable, x, y) == NULL) {
				sourceEdge = kmer_edge_table_get(Graph->EdgeTable, x, v->KMer);
				assert(sourceEdge != NULL);
				destEdge = kmer_edge_table_get(Graph->EdgeTable, v->KMer, y);
				assert(destEdge != NULL);
				weight = sourceEdge->Weight;
				length = sourceEdge->Length + destEdge->Length;
				kmer_edge_table_delete(Graph->EdgeTable, sourceEdge->Source, sourceEdge->Dest);
				kmer_edge_table_delete(Graph->EdgeTable, destEdge->Source, destEdge->Dest);
				kmer_table_delete(Graph->VertexTable, v->KMer);
				ret = kmer_edge_table_insert_ex(Graph->EdgeTable, x, y, &dummy);
				if (ret == ERR_SUCCESS) {
					dummy->Length = length;
					dummy->Weight = weight;
					v = kmer_table_get(Graph->VertexTable, x);
					assert(v != NULL);
					v->AdvancedInfo.PassThroughVertex.Output = y;
					v = kmer_table_get(Graph->VertexTable, y);
					assert(v != NULL);
					v->AdvancedInfo.PassThroughVertex.Input = x;
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


ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	do {
		ret = kmer_table_insert(Graph->VertexTable, KMer);
		if (ret == ERR_TABLE_FULL)
			kmer_table_extend(Graph->VertexTable);
	} while (ret == ERR_TABLE_FULL);

	if (ret == ERR_SUCCESS)
		Graph->NumberOfVertices++;

	return ret;
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, PKMER_EDGE *Edge)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	do {
		ret = kmer_edge_table_insert_ex(Graph->EdgeTable, Source, Dest, &edge);
		if (ret == ERR_SUCCESS) {
			PKMER_TABLE_ENTRY u = NULL;
			PKMER_TABLE_ENTRY v = NULL;

			edge->Weight = weight;
			edge->Length = Length;
			u = kmer_table_get(Graph->VertexTable, edge->Source);
			v = kmer_table_get(Graph->VertexTable, edge->Dest);
			assert(u != NULL);
			assert(v != NULL);
			u->degreeOut++;
			u->AdvancedInfo.PassThroughVertex.Output = v->KMer;
			v->DegreeIn++;
			v->AdvancedInfo.PassThroughVertex.Input = u->KMer;
			if (v->Order <= u->Order)
				Graph->NumberOfBackwardEdges++;

			*Edge = edge;
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
