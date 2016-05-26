
#ifndef __KMER_GRAPH_TYPES_H__
#define __KMER_GRAPH_TYPES_H__


#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "read-info.h"

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
	char *Seq;
	size_t SeqLen;
	EKMerEdgeType SeqType;
	long Seq1Weight;
	char *Seq2;
	size_t Seq2Len;
	EKMerEdgeType Seq2Type;
	long Seq2Weight;
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
	uint32_t RefSeqPosition;
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

typedef struct _KMER_GRAPH;
typedef void (GRAPH_ON_DELETE_EDGE_CALLBACK)(const struct _KMER_GRAPH *Graph, const KMER_EDGE *Edge, void *Context);

typedef struct _KMER_GRAPH {
	uint32_t KMerSize;
	uint32_t NumberOfVertices;
	uint32_t NumberOfEdges;
	uint32_t TypedEdgeCount[kmetMax];
	PKMER_TABLE VertexTable;
	PKMER_EDGE_TABLE EdgeTable;
	PKMER_VERTEX StartingVertex;
	PKMER_VERTEX EndingVertex;
	GRAPH_ON_DELETE_EDGE_CALLBACK *DeleteEdgeCallback;
	void *DeleteEdgeCallbackContext;
} KMER_GRAPH, *PKMER_GRAPH;






#endif
