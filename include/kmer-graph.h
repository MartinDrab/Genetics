
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "read-info.h"
#include "found-sequence.h"


typedef struct _KMER_VERTEX;

typedef enum _EKMerVertexType {
	kmvtRefSeqStart,
	kmvtRefSeqMiddle,
	kmvtRefSeqEnd,
	kmvtRead,
} EKMerVertexType, PEKMerVertexType;


typedef enum _EKMerEdgeType {
	kmetReference,
	kmetRead,
	kmetVariant,
	kmetMax,
} EKMerEdgeType, *PEKMerEdgeType;


/** Represents one edge in a kmer graph. */
typedef struct _KMER_EDGE {
	struct _KMER_VERTEX *Source;
	struct _KMER_VERTEX *Dest;
	/** Number of reads going through the edge. */
	long Weight;
	/** Number of kmers the edge skips in order to avoid cycles. Maximum value is
	the kmer size minus one. */
	uint32_t Length;
	/** Edge creation order. */
	unsigned int Order;
	EKMerEdgeType Type;
	uint32_t PassCount;
	uint32_t MaxPassCount;
	char *Seq;
	size_t SeqLen;
	char *Seq2;
	size_t Seq2Len;
	READ_INFO ReadInfo;
} KMER_EDGE, *PKMER_EDGE;

POINTER_ARRAY_TYPEDEF(KMER_EDGE);
POINTER_ARRAY_IMPLEMENTATION(KMER_EDGE)

typedef struct _KMER_VERTEX {
	PKMER KMer;
	uint32_t Order;
	EKMerVertexType Type;
	POINTER_ARRAY_KMER_EDGE Successors;
	POINTER_ARRAY_KMER_EDGE Predecessors;
	boolean Finished;
	uint32_t LastRSOrder;
} KMER_VERTEX, *PKMER_VERTEX;

POINTER_ARRAY_TYPEDEF(KMER_VERTEX);
POINTER_ARRAY_IMPLEMENTATION(KMER_VERTEX)

typedef struct _KMER_VERTEX_PAIR {
	PKMER_VERTEX U;
	PKMER_VERTEX V;
} KMER_VERTEX_PAIR, *PKMER_VERTEX_PAIR;

GEN_ARRAY_TYPEDEF(KMER_VERTEX_PAIR);
GEN_ARRAY_IMPLEMENTATION(KMER_VERTEX_PAIR)

typedef struct _KMER_EDGE_PAIR {
	PKMER_EDGE U;
	PKMER_EDGE V;
	uint32_t ReadDistance;
} KMER_EDGE_PAIR, *PKMER_EDGE_PAIR;

GEN_ARRAY_TYPEDEF(KMER_EDGE_PAIR);
GEN_ARRAY_IMPLEMENTATION(KMER_EDGE_PAIR)

typedef struct _KMER_GRAPH {
	uint32_t KMerSize;
	uint32_t NumberOfVertices;
	uint32_t NumberOfEdges;
	uint32_t TypedEdgeCount[kmetMax];
	PKMER_TABLE VertexTable;
	PKMER_EDGE_TABLE EdgeTable;
	PKMER_VERTEX StartingVertex;
	PKMER_VERTEX EndingVertex;
} KMER_GRAPH, *PKMER_GRAPH;


#define kmer_graph_get_kmer_size(aGraph)					((aGraph)->KMerSize)
#define kmer_graph_get_edge_count(aGraph)					((aGraph)->NumberOfEdges)
#define kmer_graph_get_vertex_count(aGraph)					((aGraph)->NumberOfVertices)
#define kmer_graph_get_cycle_count(aGraph)					((aGraph)->NumberOfBackwardEdges)

#define kmer_vertex_get_succ_edge(aVertex, aIndex)			(*pointer_array_item_KMER_EDGE(&(aVertex)->Successors, (aIndex)))
#define kmer_vertex_get_pred_edge(aVertex, aIndex)			(*pointer_array_item_KMER_EDGE(&(aVertex)->Predecessors, (aIndex)))
#define kmer_vertex_get_successor(aVertex, aIndex)			(kmer_vertex_get_succ_edge((aVertex), (aIndex))->Dest)
#define kmer_vertex_get_predecessor(aVertex, aIndex)		(kmer_vertex_get_pred_edge((aVertex), (aIndex))->Source)
#define kmer_vertex_in_degree(aVertex)						(pointer_array_size(&(aVertex)->Predecessors))
#define kmer_vertex_out_degree(aVertex)						(pointer_array_size(&(aVertex)->Successors))


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
void kmer_graph_print(FILE *Stream, const KMER_GRAPH *Graph);
void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer);

ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);
void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold);
void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph, size_t *DeletedThings);
ERR_VALUE kmer_graph_resolve_bubbles(PKMER_GRAPH Graph, const uint32_t Threshold);
ERR_VALUE kmer_graph_connect_reads_by_refseq(PKMER_GRAPH Graph, const size_t Threshold, size_t *ChangeCount);
ERR_VALUE kmer_graph_connect_reads_by_reads(PKMER_GRAPH Graph, const size_t Threshold);
ERR_VALUE kmer_graph_connect_reads_by_pairs(PKMER_GRAPH Graph, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, size_t *ChangeCount);
ERR_VALUE kmer_graph_detect_uncertainities(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type);
ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge);
ERR_VALUE kmer_graph_delete_vertex(PKMER_GRAPH Graph, PKMER_VERTEX Vertex);
void kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge);
ERR_VALUE kmer_graph_merge_edges(PKMER_GRAPH Graph, PKMER_EDGE Source, PKMER_EDGE Dest);
ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, PPOINTER_ARRAY_FOUND_SEQUENCE SeqArray);
ERR_VALUE kmer_vertex_get_certain_edges(const KMER_VERTEX *Vertex, const EKMerEdgeType EdgeType, const boolean Incomming, PPOINTER_ARRAY_KMER_EDGE Array);
void kmer_graph_resolve_db_triangles(PKMER_GRAPH Graph, const size_t Threshold);

PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest);
PKMER_VERTEX kmer_graph_get_vertex(const struct _KMER_GRAPH *Graph, const struct _KMER *KMer);
ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PPOINTER_ARRAY_KMER_VERTEX VertexArray);



#endif 
