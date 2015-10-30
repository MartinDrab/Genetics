
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


void kmer_graph_print(const PKMER_GRAPH Graph)
{
	printf("digraph G {\n");
	printf("\t/* number of vertices: %u */\n", Graph->NumberOfVertices);
	printf("\t/* number of edges: %u */\n", Graph->NumberOfEdges);
	kmer_table_print(Graph->VertexTable);
	kmer_edge_table_print(Graph->EdgeTable);
	printf("}\n");

	return;
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
			PKMER_TABLE_ENTRY v = NULL;

			edge->Weight = weight;
			edge->Length = Length;
			v = kmer_table_get(Graph->VertexTable, edge->Source);
			assert(v != NULL);
			v->degreeOut++;
			v = kmer_table_get(Graph->VertexTable, edge->Dest);
			assert(v != NULL);
			v->DegreeIn++;
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
