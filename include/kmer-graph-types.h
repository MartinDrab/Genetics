
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
typedef struct _KMER_GRAPH;

typedef enum _EKMerVertexType {
	kmvtRefSeqStart,
	kmvtRefSeqMiddle,
	kmvtRefSeqEnd,
	kmvtRead,
	kmvtMax
} EKMerVertexType, PEKMerVertexType;


typedef enum _EKMerEdgeType {
	kmetReference,
	kmetRead,
	kmetVariant,
	kmetNone,
	kmetMax,
} EKMerEdgeType, *PEKMerEdgeType;


/** Represents one edge in a kmer graph. */
typedef struct _KMER_EDGE {
	struct _KMER_VERTEX *Source;
	struct _KMER_VERTEX *Dest;
	/** Edge creation order. */
	unsigned int Order;
	EKMerEdgeType Type;
	char *Seq;
	size_t SeqLen;
	EKMerEdgeType SeqType;
	size_t Seq1Weight;
	char *Seq2;
	size_t Seq2Len;
	EKMerEdgeType Seq2Type;
	size_t Seq2Weight;
	READ_INFO ReadInfo;
	boolean MarkedForDelete;
	boolean Finished;
	GEN_ARRAY_size_t Paths;
} KMER_EDGE, *PKMER_EDGE;

POINTER_ARRAY_TYPEDEF(KMER_EDGE);
POINTER_ARRAY_IMPLEMENTATION(KMER_EDGE)

typedef struct _KMER_VERTEX {
	uint32_t Order;
	boolean Helper;
	EKMerVertexType Type;
	POINTER_ARRAY_KMER_EDGE Successors;
	POINTER_ARRAY_KMER_EDGE Predecessors;
	uint32_t RefSeqPosition;
	union {
		struct _KMER_VERTEX *Next;
		struct _KMER_GRAPH *Graph;
	} Lists;
	KMER KMer;
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
	PKMER_EDGE ConnectingEdge;
	uint32_t ReadDistance;
	PKMER_EDGE *Edges;
	size_t EdgeCount;
} KMER_EDGE_PAIR, *PKMER_EDGE_PAIR;

typedef struct _KMER_EDGE_PAIR_KEY {
	PKMER_EDGE U;
	PKMER_EDGE V;
	PKMER_EDGE ConnectingEdge;
	uint32_t ReadDistance;
} KMER_EDGE_PAIR_KEY, *PKMER_EDGE_PAIR_KEY;

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
	PKMER_EDGE_TABLE DummyVertices;
	PKMER_VERTEX StartingVertex;
	PKMER_VERTEX EndingVertex;
	GRAPH_ON_DELETE_EDGE_CALLBACK *DeleteEdgeCallback;
	void *DeleteEdgeCallbackContext;
	PKMER_VERTEX VerticesToDeleteList;
} KMER_GRAPH, *PKMER_GRAPH;






#endif
