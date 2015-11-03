
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"




typedef struct _KMER_GRAPH {
	uint32_t KMerSize;
	uint32_t NumberOfVertices;
	uint32_t NumberOfEdges;
	uint32_t NumberOfBackwardEdges;
	PKMER_TABLE VertexTable;
	PKMER_EDGE_TABLE EdgeTable;
} KMER_GRAPH, *PKMER_GRAPH;

#define kmer_graph_get_kmer_size(aGraph)					((aGraph)->KMerSize)

ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_copy(const PKMER_GRAPH Source, PKMER_GRAPH *Copied);
void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer);
ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, PKMER_EDGE *Edge);

PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest);


#endif 
