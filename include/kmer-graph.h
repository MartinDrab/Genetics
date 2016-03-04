
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"

typedef struct _KMER_EDGE;

typedef enum _EKMerVertexType {
	kmvtRefSeqStart,
	kmvtRefSeqMiddle,
	kmvtRefSeqEnd,
	kmvtRead,
} EKMerVertexType, PEKMerVertexType;

typedef enum _EVertexPassType {
	vptRefSeq,
	vptRead,
} EVertexPassType, *PEVertexPassType;

typedef struct _KMER_VERTEX_PASS {
	const struct _KMER_EDGE *Incomming;
	const struct _KMER_EDGE *Outgoing;
	EVertexPassType PassType;
	union {
		struct {
			size_t ReadIndex;
		} ReadPass;
	} Data;
} KMER_VERTEX_PASS, *PKMER_VERTEX_PASS;

typedef struct _KMER_VERTEX {
	PKMER KMer;
	uint32_t DegreeIn;
	uint32_t degreeOut;
	uint32_t Order;
	EKMerVertexType Type;
	DYM_ARRAY Successors;
	DYM_ARRAY Predecessors;
	
	PKMER_VERTEX_PASS Passes;
	size_t PassCount;
	PKMER_VERTEX_PASS CurrentPass;
	size_t StartIndex;
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
	char *Seq;
	size_t SeqLen;
	DYM_ARRAY ReadIndices;
	DYM_ARRAY ReadPositions;
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
#define kmer_graph_get_edge_count(aGraph)					((aGraph)->NumberOfEdges)
#define kmer_graph_get_vertex_count(aGraph)					((aGraph)->NumberOfVertices)
#define kmer_graph_get_cycle_count(aGraph)					((aGraph)->NumberOfBackwardEdges)

#define kmer_vertex_get_succ_edge(aVertex, aIndex)			((PKMER_EDGE)(dym_array_data(&(aVertex)->Successors)[(aIndex)]))
#define kmer_vertex_get_pred_edge(aVertex, aIndex)			((PKMER_EDGE)(dym_array_data(&(aVertex)->Predecessors)[(aIndex)]))
#define kmer_vertex_get_successor(aVertex, aIndex)			(kmer_vertex_get_succ_edge((aVertex), (aIndex))->Dest)
#define kmer_vertex_get_predecessor(aVertex, aIndex)		(kmer_vertex_get_pred_edge((aVertex), (aIndex))->Source)
#define kmer_vertex_get_pass_count(aVertex)					((aVertex)->PassCount)
#define kmer_vertex_get_pass(aVertex, aIndex)				((aVertex)->Passes + (aIndex))

ERR_VALUE kmer_vertex_add_pass(PKMER_VERTEX Vertex, const struct _KMER_EDGE *Incomming, const struct _KMER_EDGE *Outgoing, const EVertexPassType PassType, const size_t Data);
void kmer_vertex_remove_pass(PKMER_VERTEX Vertex, const size_t Index);
void kmer_vertex_remove_passes(PKMER_VERTEX Vertex, const struct _KMER_EDGE *Source, const struct _KMER_EDGE *Dest);

ERR_VALUE kmer_edge_add_read(PKMER_EDGE Edge, const size_t ReadIndex, const size_t ReadPosition);
ERR_VALUE kmer_edge_connecting_reads(const struct _KMER_EDGE **EdgeSet, const size_t EdgeSetSize, PDYM_ARRAY ReadIndices);
#define kmer_edge_get_read_count(aEdge)				dym_array_size(&(aEdge)->ReadIndices)
void kmer_edge_get_read(const KMER_EDGE *Edge, const size_t Index, size_t *ReadIndex, size_t *ReadPosition);
ERR_VALUE kmer_edge_find_read(const KMER_EDGE *Edge, const size_t StartIndex, size_t *EndIndex, const size_t ReadIndex, size_t *ReadPosition);
void kmer_edge_remove_read(PKMER_EDGE Edge, const size_t StartIndex, const size_t ReadIndex);

ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph);
void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer);

ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);
void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold);
void kmer_graph_compute_edge_probablities(PKMER_GRAPH Graph);
void kmer_graph_compute_shurtcuts(PKMER_GRAPH Graph, const size_t MaxLength);
void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph);
void kmer_graph_separate_distinct_passes(PKMER_GRAPH Graph);
void kmer_graph_resolve_repetitions(PKMER_GRAPH Graph);
ERR_VALUE kmer_graph_connect_reads(PKMER_GRAPH Graph, const size_t Threshold, const ONE_READ *Reads);
ERR_VALUE kmer_graph_merge_unbranched_refseq(PKMER_GRAPH Graph);

ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer, const EKMerVertexType Type);
ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge);
ERR_VALUE kmer_graph_get_seq(const KMER_GRAPH *Graph, char **Seq, size_t *SeqLen);
ERR_VALUE kmer_graph_delete_vertex(PKMER_GRAPH Graph, PKMER_VERTEX Vertex);
ERR_VALUE kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge);
ERR_VALUE kmer_graph_merge_edges(PKMER_GRAPH Graph, PKMER_EDGE Source, PKMER_EDGE Dest);
ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, PDYM_ARRAY SeqArray);
ERR_VALUE kmer_graph_find_bubble(PKMER_GRAPH Graph, PKMER_EDGE RefSeqEdge, PKMER_EDGE ReadEdge, PKMER_VERTEX *EndVertex, boolean *RefSeqSupported, const uint32_t Threshold);

PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest);
PKMER_VERTEX kmer_graph_get_vertex(const struct _KMER_GRAPH *Graph, const struct _KMER *KMer);



#endif 
