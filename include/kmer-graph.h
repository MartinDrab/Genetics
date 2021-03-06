
#ifndef __KMER_GRAPH_H__
#define __KMER_GRAPH_H__


#include <stdint.h>
#include "err.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "read-info.h"
#include "kmer-graph-types.h"
#include "variant.h"
#include "assembly.h"


UTILS_TYPED_CALLOC_FUNCTION(PKMER_VERTEX)
UTILS_TYPED_CALLOC_FUNCTION(PKMER_EDGE)
UTILS_NAMED_CALLOC_FUNCTION(PPKMER_VERTEX, PKMER_VERTEX *)
UTILS_NAMED_CALLOC_FUNCTION(PPKMER_EDGE, PKMER_EDGE *)


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
PKMER_EDGE kmer_vertex_get_edge_by_base(PKMER_VERTEX Vertex, const char Base);


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, const size_t VerticesHint, const size_t EdgesHint, PKMER_GRAPH *Graph);
void kmer_graph_destroy(PKMER_GRAPH Graph);
void kmer_graph_print(FILE *Stream, const KMER_GRAPH *Graph);
void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer);
void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer);

void kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph);
void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const size_t Threshold);
void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph, size_t *DeletedThings);
ERR_VALUE kmer_graph_connect_reads_by_pairs(PKMER_GRAPH Graph, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, size_t *ChangeCount);
ERR_VALUE kmer_graph_detect_variant(PKMER_GRAPH Graph, PGEN_ARRAY_VARIANT_CALL VCArray, const char *CHrom, const PARSE_OPTIONS *Options, boolean *Changed);

ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex);
ERR_VALUE kmer_graph_add_helper_vertex(PKMER_GRAPH Graph, const KMER *KMer1, const KMER *KMer2, PKMER_VERTEX *Vertex);
ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, PKMER_EDGE *Edge);
ERR_VALUE kmer_graph_delete_vertex(PKMER_GRAPH Graph, PKMER_VERTEX Vertex);
void kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge);
ERR_VALUE kmer_graph_merge_edges(PKMER_GRAPH Graph, PKMER_EDGE Source, PKMER_EDGE Dest);
ERR_VALUE kmer_graph_get_splitted_edge(PKMER_GRAPH Graph, const KMER_VERTEX *Source, const KMER_VERTEX *Dest, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex);
ERR_VALUE kmer_graph_split_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex);
void kmer_edge_add_seq(PKMER_EDGE Edge, EKMerEdgeType Type, const char *Seq, const size_t Length);

PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest);
ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PPOINTER_ARRAY_KMER_VERTEX *VertexArray);

PKMER_EDGE _get_refseq_edge(const KMER_VERTEX *Vertex);
PKMER_EDGE _get_refseq_or_variant_edge(const KMER_VERTEX *Vertex);

ERR_VALUE kmer_graph_compute_weights(PKMER_GRAPH Graph);

ERR_VALUE kmer_edge_add_read(PKMER_EDGE Edge, size_t ReadIndex, size_t ReadPosition, uint8_t Quality);


#endif 
