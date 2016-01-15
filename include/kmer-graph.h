
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

typedef enum _EKMerEdgeType {
	kmetReference,
	kmetRead,
} EKMerEdgeType, *PEKMerEdgeType;


/** Represents one edge in a kmer graph. */
typedef struct _KMER_EDGE {
	PKMER_VERTEX Source;
	PKMER_VERTEX Dest;
	/** Number of reads going through the edge. */
	long Weight;
	/** Number of kmers the edge skips in order to avoid cycles. Maximum value is
	the kmer size minus one. */
	uint32_t Length;
	/** Edge creation order. */
	unsigned int Order;
	EKMerEdgeType Type;
	double Probability;
	uint32_t PassCount;
	uint32_t MaxPassCount;
	void *Shortcut;
} KMER_EDGE, *PKMER_EDGE;

typedef struct _KMER_GRAPH_SHORTCUT {
	PKMER_VERTEX StartVertex;
	PKMER_VERTEX EndVertex;
	char *Sequence;
	uint32_t Length;
	double Probability;
	uint32_t PassCount;
	uint32_t MaxPassCount;
} KMER_GRAPH_SHORTCUT, *PKMER_GRAPH_SHORTCUT;


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

#define kmer_vertex_get_succ_edge(aVertex, aIndex)			((PKMER_EDGE)(dym_array_data(&(aVertex)->Successors)[(aIndex)]))
#define kmer_vertex_get_pred_edge(aVertex, aIndex)			((PKMER_EDGE)(dym_array_data(&(aVertex)->Predecessors)[(aIndex)]))
#define kmer_vertex_get_successor(aVertex, aIndex)			(kmer_vertex_get_succ_edge((aVertex), (aIndex))->Dest)
#define kmer_vertex_get_predecessor(aVertex, aIndex)		(kmer_vertex_get_pred_edge((aVertex), (aIndex))->Source)


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);
void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold);
void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_compute_edge_probablities(PKMER_GRAPH Graph);
void kmer_graph_compute_shurtcuts(PKMER_GRAPH Graph, const size_t MaxLength);
void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer, const EKMerVertexType Type);
ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge);

PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest);


#endif 
