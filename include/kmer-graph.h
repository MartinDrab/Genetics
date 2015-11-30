
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"


typedef enum _EKMerVertexType {
	kmvtRefSeqStart,
	kmvtRefSeqMiddle,
	kmvtRefSeqEnd,
	kmvtRead,
} EKMerVertexType, PEKMerVertexType;

typedef struct _KMER_VERTEX {
	PKMER KMer;
	uint32_t DegreeIn;
	uint32_t degreeOut;
	uint32_t Order;
	EKMerVertexType Type;
	DYM_ARRAY Successors;
	DYM_ARRAY Predecessors;
} KMER_VERTEX, *PKMER_VERTEX;

typedef struct _KMER_GRAPH {
	uint32_t KMerSize;
	uint32_t NumberOfVertices;
	uint32_t NumberOfEdges;
	uint32_t NumberOfBackwardEdges;
	PKMER_TABLE VertexTable;
	PKMER_EDGE_TABLE EdgeTable;
	PKMER_VERTEX StartingVertex;
	PKMER_VERTEX EndingVertex;
} KMER_GRAPH, *PKMER_GRAPH;

#define kmer_graph_get_kmer_size(aGraph)					((aGraph)->KMerSize)

ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_copy(const PKMER_GRAPH Source, PKMER_GRAPH *Copied);
void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);
void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold);
void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_compute_edge_probablities(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer, const EKMerVertexType Type);
ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, PKMER_EDGE *Edge);

PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest);


#endif 
