
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"




typedef struct _KMER_GRAPH {
	PKMER_TABLE VertexTable;
	PKMER_EDGE_TABLE EdgeTable;
} KMER_GRAPH, *PKMER_GRAPH;


ERR_VALUE kmer_graph_create(const size_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer);
ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight);
PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest);


#endif 
