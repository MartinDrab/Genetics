
#ifndef __KMER_GRAPH_TYPES_H__
#define __KMER_GRAPH_TYPES_H__


#include "utils.h"
#include "kmer-graph-base-types.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "read-info.h"
#include "found-sequence-types.h"

typedef struct _KMER_VERTEX;
typedef struct _KMER_GRAPH;


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
	READ_INFO ReadInfo;
	boolean MarkedForDelete;
	boolean Finished;
	GEN_ARRAY_size_t Paths;
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT Variants;
} KMER_EDGE, *PKMER_EDGE;

POINTER_ARRAY_TYPEDEF(KMER_EDGE);
POINTER_ARRAY_IMPLEMENTATION(KMER_EDGE)

typedef struct _KMER_VERTEX {
	uint32_t Order;
	boolean Helper;
	boolean LongEdgeAllowed;
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

typedef struct _KMER_GRAPH;
typedef PKMER_VERTEX (KMER_GRAPH_VERTEX_ALLOCATOR)(struct _KMER_GRAPH *Graph, void *Context);
typedef PKMER_EDGE (KMER_GRAPH_EDGE_ALLOCATOR)(struct _KMER_GRAPH *Graph, void *Context);
typedef void(KMER_GRAPH_VERTEX_FREER)(struct _KMER_GRAPH *Graph, PKMER_VERTEX Vertex, void *Context);
typedef void(KMER_GRAPH_EDGE_FREER)(struct _KMER_GRAPH *Graph, PKMER_EDGE Edge, void *Context);


typedef struct _KMER_GRAPH_ALLOCATOR {
	KMER_GRAPH_VERTEX_ALLOCATOR *VertexAllocator;
	KMER_GRAPH_VERTEX_FREER *VertexFreer;
	void *VertexAllocatorContext;
	KMER_GRAPH_EDGE_ALLOCATOR *EdgeAllocator;
	KMER_GRAPH_EDGE_FREER *EdgeFreer;
	void *EdgeAllocatorContext;
} KMER_GRAPH_ALLOCATOR, *PKMER_GRAPH_ALLOCATOR;

typedef struct _KMER_LIST {
	POINTER_ARRAY_KMER_VERTEX Vertices;
	KMER Kmer;
} KMER_LIST, *PKMER_LIST;

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
	PKMER_TABLE KmerListTable;
	uint8_t QualityTable[256];
	KMER_GRAPH_ALLOCATOR Allocator;
} KMER_GRAPH, *PKMER_GRAPH;






#endif
