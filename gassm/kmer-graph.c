
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
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/

ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	tmpGraph = (PKMER_GRAPH)malloc(sizeof(KMER_GRAPH));
	if (tmpGraph != NULL) {
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
			free(tmpGraph);
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

void kmer_graph_destroy(PKMER_GRAPH Graph)
{
	kmer_edge_table_destroy(Graph->EdgeTable);
	kmer_table_destroy(Graph->VertexTable);
	free(Graph);

	return;
}

void kmer_graph_print(const PKMER_GRAPH Graph)
{
	printf("digraph G {\n");
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

	return ret;
}

ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE edge = NULL;

	do {
		ret = kmer_edge_table_insert(Graph->EdgeTable, Source, Dest);
		if (ret == ERR_SUCCESS) {
			edge = kmer_edge_table_get(Graph->EdgeTable, Source, Dest);
			assert(edge != NULL);
			edge->Weight = weight;
		}

		if (ret == ERR_TABLE_FULL)
			kmer_edge_table_extend(Graph->EdgeTable);
	} while (ret == ERR_TABLE_FULL);

	return ret;
}

PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest)
{
	return kmer_edge_table_get(Graph->EdgeTable, Source, Dest);
}
