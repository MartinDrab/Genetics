
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
#include "variant-types.h"

struct _KMER_VERTEX;
struct _KMER_GRAPH;


/** Represents one edge in a kmer graph. */
typedef struct _KMER_EDGE {
	struct _KMER_VERTEX *Source;
	struct _KMER_VERTEX *Dest;
	/** Edge creation order. */
	unsigned int Order;
	/** Edge type (read, reference, variant). */
	EKMerEdgeType Type;
	/** Sequence remembered on the edge. */
	const char *Seq;
	/** Length of the edge sequence. */
	size_t SeqLen;
	/** Not used. */
	EKMerEdgeType SeqType;
	/** Not used. */
	size_t Seq1Weight;
	/** Indices, positions and base qualities of reads covering the edge. */
	READ_INFO ReadInfo;
	boolean MarkedForDelete;
	/** Variants owned by the edge. */
	POINTER_ARRAY_VARIANT_CALL VCs;
	GEN_ARRAY_size_t Weights;
	/** Indices of reads covering the edge (this information is used
	    later when the edge is merged with other edges). */
	POINTER_ARRAY_READ_INFO ReadIndices;
	struct {
		/** Indicates whether the edge is a long/connecting edge. */
		boolean LongEdge;
		/** Reference position where the connecting edge end. */
		uint32_t RefSeqEnd;
		/** Reference position where the connecting edge ends. */
		uint32_t RefSeqStart;
	} LongData;
} KMER_EDGE, *PKMER_EDGE;

POINTER_ARRAY_TYPEDEF(KMER_EDGE);
POINTER_ARRAY_IMPLEMENTATION(KMER_EDGE)


/** Represents one vertex of the de Bruijn graph. */
typedef struct _KMER_VERTEX {
	uint32_t Order;
	/** Indicates whether this is a helper vertex. */
	boolean Helper;
	/** Can connecting/long edges start here? */
	boolean LongEdgeAllowed;
	/** Can this be the first vertex of a read? */
	boolean ReadStartAllowed;
	/** Indicates whether the vertex k-mer string is unique within the de Bruijn graph. */
	boolean Unique;
	/** Not used. */
	boolean ShortVariant;
	/** Type of the vertex (start, end, reference, read). */
	EKMerVertexType Type;
	/** Reference edge pointing to the next reference vertex. */
	PKMER_EDGE RefEdge;
	/** Reference or variant edge leading to the next reference vertex. */
	PKMER_EDGE RefVarEdge;
	/** Outgoing edges. */
	POINTER_ARRAY_KMER_EDGE Successors;
	/** Incoming edges. */
	POINTER_ARRAY_KMER_EDGE Predecessors;
	/** Reference position of the vertex, relative to the active region. */
	uint32_t RefSeqPosition;
	/** Absolute reference position of the vertex. */
	uint64_t AbsPos;
	union {
		struct _KMER_VERTEX *Next;
		struct _KMER_GRAPH *Graph;
	} Lists;
	/** K-mer representing the vertex. Must be always last in the structure, since
	    it may occupy more bytes than this definition indicates. */
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

/** Maintains information about a connecting edge. */
typedef struct _KMER_EDGE_PAIR {
	/** The input edge. */
	PKMER_EDGE U;
	/** The output edge. */
	PKMER_EDGE V;
	/** The connecting edge. */
	PKMER_EDGE ConnectingEdge;
	/** Number of bases the connection edge skips. */
	size_t ReadDistance;
	/** Edges (not including the input and output one) the connection edge skips/replaces. */
	PKMER_EDGE *Edges;
	/** Number of skipped edges (not including the input and output one). */
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

struct _KMER_GRAPH;
typedef void (GRAPH_ON_DELETE_EDGE_CALLBACK)(const struct _KMER_GRAPH *Graph, const KMER_EDGE *Edge, void *Context);

struct _KMER_GRAPH;
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

/** Stores k-mers equal by sequence (but differrent by their context numbers). */
typedef struct _KMER_LIST {
	/** List of the vertices sharing the same k-mer string. */
	POINTER_ARRAY_KMER_VERTEX Vertices;
	/** The k-mer string. */
	KMER Kmer;
} KMER_LIST, *PKMER_LIST;


/** Represents the de Bruijn graph. */
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
	POINTER_ARRAY_KMER_VERTEX RefVertices;
	KMER_GRAPH_ALLOCATOR Allocator;
} KMER_GRAPH, *PKMER_GRAPH;






#endif
