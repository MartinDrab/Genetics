
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "err.h"
#include "utils.h"
#include "dym-array.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"
#include "kmer-graph.h"


/************************************************************************/
/*                      HELPER FUNCTIONS                                */
/************************************************************************/


static UTILS_TYPED_MALLOC_FUNCTION(KMER_GRAPH)
static UTILS_TYPED_MALLOC_FUNCTION(KMER_VERTEX)
static UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE)

/************************************************************************/
/*                        VERTEX BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _vertex_create(const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_VERTEX(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->KMer = kmer_copy(KMer);
		ret = (tmp->KMer != NULL) ? ERR_SUCCESS : ERR_OUT_OF_MEMORY;
		if (ret == ERR_SUCCESS) {
			tmp->DegreeIn = 0;
			tmp->degreeOut = 0;
			tmp->Type = Type;
			dym_array_create(&tmp->Successors, 140);
			dym_array_create(&tmp->Predecessors, 140);
			tmp->Finished = FALSE;
			*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmp);
	}

	return ret;
}


static void _vertex_destroy(PKMER_VERTEX Vertex)
{
	dym_array_destroy(&Vertex->Predecessors);
	dym_array_destroy(&Vertex->Successors);
	kmer_free(Vertex->KMer);
	utils_free(Vertex);

	return;
}


static ERR_VALUE _vertex_copy(const KMER_VERTEX *Vertex, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(Vertex->KMer, Vertex->Type, &tmp);
	if (ret == ERR_SUCCESS) {
		tmp->DegreeIn = Vertex->DegreeIn;
		tmp->degreeOut = Vertex->degreeOut;
		tmp->Order = Vertex->Order;
		dym_array_create(&tmp->Successors, 140);
		ret = dym_array_copy(&tmp->Successors, &Vertex->Successors);
		if (ret == ERR_SUCCESS) {
			ret = dym_array_copy(&tmp->Predecessors, &Vertex->Predecessors);
			if (ret == ERR_SUCCESS)
				*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(tmp);
	}

	return ret;
}

/************************************************************************/
/*                        EDGE BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _edge_create(PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, const long Weight, const uint32_t Length, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE(&tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Source = Source;
		tmp->Dest = Dest;
		tmp->Type = Type;
		tmp->Length = Length;
		tmp->Weight = Weight;
		tmp->PassCount = 0;
		tmp->MaxPassCount = 0;
		tmp->Probability = 0;
		tmp->Order = 0;
		tmp->Shortcut = NULL;
		tmp->Seq = NULL;
		tmp->SeqLen = 0;
		*Edge = tmp;
	}

	return ret;
}


static void _edge_destroy(PKMER_EDGE Edge)
{
	if (Edge->SeqLen > 0)
		utils_free(Edge->Seq);

	utils_free(Edge);

	return;
}


static ERR_VALUE _edge_copy(const KMER_EDGE *Edge, PKMER_EDGE *Result)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _edge_create(Edge->Source, Edge->Dest, Edge->Type, Edge->Weight, Edge->Length, &tmp);
	if (ret == ERR_SUCCESS) {
		if (Edge->SeqLen > 0)
			ret = utils_copy_string(Edge->Seq, &tmp->Seq);

		if (ret == ERR_SUCCESS) {
			tmp->SeqLen = Edge->SeqLen;
			tmp->PassCount = Edge->PassCount;
			tmp->MaxPassCount = Edge->MaxPassCount;
			tmp->Probability = Edge->Probability;
			tmp->Order = Edge->Order;
			*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			_edge_destroy(tmp);
	}

	return ret;
}


/************************************************************************/
/*                      VERTEX TABLE CALBLACKS                          */
/************************************************************************/


static void _vertex_table_on_insert(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order)
{
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;

	v->Order = Order;

	return;
}

static void _vertex_table_on_delete(struct _KMER_TABLE *Table, void *ItemData)
{
	_vertex_destroy((PKMER_VERTEX)ItemData);

	return;
}


static ERR_VALUE _vertex_table_on_copy(struct _KMER_TABLE *Table, void *ItemData, void **Copy)
{
	return _vertex_copy((PKMER_VERTEX)ItemData, (PKMER_VERTEX *)Copy);
}


static void _vertex_table_on_print(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream)
{
	char *colors[] = {"yellow", "green", "blue", "red"};
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;


	fprintf(Stream, "\t");
	kmer_print(Stream, v->KMer);
	fprintf(Stream, "[label=\"");
	kmer_print(Stream, v->KMer);
	fprintf(Stream, "\\nIN=%u; OUT=%u\",style=filled,color=%s]", v->DegreeIn, v->degreeOut, colors[v->Type]);
	fprintf(Stream, ";\n");

	return;
}


/************************************************************************/
/*                      EDGE TABLE CALBLACKS                          */
/************************************************************************/


static void _edge_table_on_insert(struct _KMER_EDGE_TABLE *Table, void *ItemData, const uint32_t Order)
{
	PKMER_EDGE e = (PKMER_EDGE)ItemData;

	e->Order = Order;

	return;
}

static void _edge_table_on_delete(struct _KMER_EDGE_TABLE *Table, void *ItemData)
{
	_edge_destroy((PKMER_EDGE)ItemData);

	return;
}


static ERR_VALUE _edge_table_on_copy(struct _KMER_EDGE_TABLE *Table, void *ItemData, void **Copy)
{
	return _edge_copy((PKMER_EDGE)ItemData, (PKMER_EDGE *)Copy);
}


static void _edge_table_on_print(struct _KMER_EDGE_TABLE *Table, void *ItemData, FILE *Stream)
{
	PKMER_EDGE e = (PKMER_EDGE)ItemData;

	fprintf(Stream, "\t");
	kmer_print(Stream, e->Source->KMer);
	fprintf(Stream, " -> ");
	kmer_print(Stream, e->Dest->KMer);
	fprintf(Stream, " [");
	if (e->Type == kmetReference)
		fprintf(Stream, "color=green");
	else fprintf(Stream, "color=red");

	fprintf(Stream, ",label=\"W: %li; M: %u; L %u\"", e->Weight, e->MaxPassCount, e->Length);
	fprintf(Stream, "];\n");

	return;
}


static boolean _kmer_vertex_no_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = TRUE;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		ret = (kmer_vertex_get_succ_edge(Vertex, i))->Type != kmetRead;
		if (!ret)
			break;
	}

	if (ret) {
		for (size_t i = 0; i < Vertex->DegreeIn; ++i) {
			ret = (kmer_vertex_get_pred_edge(Vertex, i))->Type != kmetRead;
			if (!ret)
				break;
		}
	}

	return ret;
}


/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	KMER_TABLE_CALLBACKS vCallbacks;
	KMER_EDGE_TABLE_CALLBACKS eCallbacks;

	ret = utils_malloc_KMER_GRAPH(&tmpGraph);
	if (ret == ERR_SUCCESS) {
		tmpGraph->NumberOfEdges = 0;
		tmpGraph->NumberOfVertices = 0;
		tmpGraph->KMerSize = KMerSize;
		tmpGraph->NumberOfBackwardEdges = 0;
		tmpGraph->StartingVertex = NULL;
		tmpGraph->EndingVertex = NULL;
		vCallbacks.OnCopy = _vertex_table_on_copy;
		vCallbacks.OnDelete = _vertex_table_on_delete;
		vCallbacks.OnInsert = _vertex_table_on_insert;
		vCallbacks.OnPrint = _vertex_table_on_print;
		ret = kmer_table_create(KMerSize, 2, 37, &vCallbacks, &tmpGraph->VertexTable);
		if (ret == ERR_SUCCESS) {
			eCallbacks.OnCopy = _edge_table_on_copy;
			eCallbacks.OnDelete = _edge_table_on_delete;
			eCallbacks.OnInsert = _edge_table_on_insert;
			eCallbacks.OnPrint = _edge_table_on_print;
			ret = kmer_edge_table_create(KMerSize, 2, 37, &eCallbacks, &tmpGraph->EdgeTable);
			if (ret == ERR_SUCCESS)
				*Graph = tmpGraph;

			if (ret != ERR_SUCCESS)
				kmer_table_destroy(tmpGraph->VertexTable);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpGraph);
	}

	return ret;
}


void kmer_graph_destroy(PKMER_GRAPH Graph)
{
	kmer_edge_table_destroy(Graph->EdgeTable);
	kmer_table_destroy(Graph->VertexTable);
	utils_free(Graph);

	return;
}


void kmer_graph_print(FILE *Stream, const PKMER_GRAPH Graph)
{
	fprintf(Stream, "digraph G {\n");
	fprintf(Stream, "\t/* number of vertices: %u */\n", Graph->NumberOfVertices);
	fprintf(Stream, "\t/* number of edges: %u */\n", Graph->NumberOfEdges);
	fprintf(Stream, "\t/* number of backward edges: %u */\n", Graph->NumberOfBackwardEdges);
	kmer_table_print(Stream, Graph->VertexTable);
	kmer_edge_table_print(Stream, Graph->EdgeTable);
	fprintf(Stream, "}\n");

	return;
}


ERR_VALUE kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph)
{
	PKMER_VERTEX v = NULL;
	PKMER_TABLE_ENTRY e = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	PKMER_VERTEX x = NULL;
	PKMER_VERTEX y = NULL;
	PKMER_EDGE dummy = NULL;
	PKMER_EDGE sourceEdge = NULL;
	PKMER_EDGE destEdge = NULL;
	long weight = 0;
	long sourceWeight = 0;
	long destWeight = 0;
	uint32_t length = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret == ERR_SUCCESS) {
		boolean last = FALSE;
		
		e = iter;
		last = (kmer_table_next(Graph->VertexTable, e, &iter) == ERR_NO_MORE_ENTRIES);
		v = (PKMER_VERTEX)e->Data;
		if (v->DegreeIn == 1 && v->degreeOut == 1)
			kmer_graph_delete_vertex(Graph, v);

		if (last)
			break;
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	return ret;
}


void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const long Threshold)
{
	PKMER_EDGE e = NULL;
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	boolean last = FALSE;
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX v = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (err == ERR_SUCCESS) {
		boolean last = FALSE;

		e = (PKMER_EDGE)(iter->Data);
		last = (kmer_edge_table_next(Graph->EdgeTable, iter, &iter) == ERR_NO_MORE_ENTRIES);
		if (e->Weight < Threshold)
			kmer_graph_delete_edge(Graph, e);

		if (last)
			break;
	}

	return;
}


void kmer_graph_set_starting_vertex(PKMER_GRAPH Graph, const KMER *KMer)
{
	assert(Graph->StartingVertex == NULL);
	Graph->StartingVertex = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	assert(Graph->StartingVertex != NULL);

	return;
}


void kmer_graph_set_ending_vertex(PKMER_GRAPH Graph, const KMER *KMer)
{
	assert(Graph->EndingVertex == NULL);
	Graph->EndingVertex = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	assert(Graph->EndingVertex != NULL);

	return;
}


static int _edge_prob_compare(const void *e1, const void *e2)
{
	const KMER_EDGE *edge1 = (PKMER_EDGE)e1;
	const KMER_EDGE *edge2 = (PKMER_EDGE)e2;

	if (edge1->Probability < edge2->Probability)
		return 1;
	else if (edge1->Probability > edge2->Probability)
		return -1;

	return 0;
}

void kmer_graph_compute_edge_probablities(PKMER_GRAPH Graph)
{
	long total = 0;
	PKMER_VERTEX v = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret != ERR_NO_MORE_ENTRIES) {
		total = 0;
		v = (PKMER_VERTEX)iter->Data;
		for (size_t i = 0; i < v->degreeOut; ++i) {
			PKMER_EDGE e = kmer_vertex_get_succ_edge(v , i);

			total += (e->Weight + 1);
		}

		for (size_t i = 0; i < v->degreeOut; ++i) {
			PKMER_EDGE e = kmer_vertex_get_succ_edge(v, i);

			e->Probability = (total != 0) ? ((double)(e->Weight + 1) / (double)total) : (1/(double)v->degreeOut);
			if (e->Probability != 0)
				e->Probability = log(e->Probability);
		}

		qsort(dym_array_data(&v->Successors), v->degreeOut, sizeof(PKMER_EDGE), _edge_prob_compare);
		ret = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}


void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph, size_t *DeletedThings)
{
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_table_first(Graph->VertexTable, &iter);
	while (ret != ERR_NO_MORE_ENTRIES) {
		boolean deleted = FALSE;
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;

		if (v->Type != kmvtRefSeqStart && v->Type != kmvtRefSeqEnd) {
			if (v->DegreeIn == 0) {
				deleted = TRUE;
				while (v->degreeOut > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_succ_edge(v, 0));
			} else if (v->degreeOut == 0) {
				deleted = TRUE;
				while (v->DegreeIn > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_pred_edge(v, 0));
			}

			if (deleted) {
				kmer_graph_delete_vertex(Graph, v);
				*DeletedThings += 1;
			}
		}

		ret = (deleted) ?
			kmer_table_first(Graph->VertexTable, &iter) :
			kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}


static ERR_VALUE _get_bubble_ends(PKMER_VERTEX FirstVertex, PDYM_ARRAY EndArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (FirstVertex->Type == kmvtRead) {
		ret = ERR_SUCCESS;
		DYM_ARRAY vertexStack;
		DYM_ARRAY edgeStack;
		PKMER_VERTEX v = FirstVertex;
		size_t edgeIndex = 0;

		ret = ERR_SUCCESS;
		dym_array_create(&vertexStack, 140);
		dym_array_create(&edgeStack, 140);
		do {
			if (edgeIndex == v->degreeOut) {
				v->Finished = FALSE;
				v = dym_array_remove_back(&vertexStack);
				edgeIndex = (size_t)dym_array_remove_back(&edgeStack);
				continue;
			}

			for (size_t i = edgeIndex; i < v->degreeOut; ++i) {
				PKMER_VERTEX tmp = kmer_vertex_get_successor(v, i);
				
				if (!tmp->Finished) {
					ret = dym_array_push_back(&vertexStack, v);
					if (ret != ERR_SUCCESS)
						break;

					ret = dym_array_push_back(&edgeStack, (void *)(i + 1));
					if (ret != ERR_SUCCESS)
						break;

					v->Finished = TRUE;
					v = tmp;
					break;
				}
			}

			if (ret == ERR_SUCCESS) {
				if (v->Type == kmvtRefSeqMiddle) {
					ret = dym_array_push_back(EndArray, v);
					if (ret != ERR_SUCCESS)
						break;

					v->Finished = FALSE;
					v = dym_array_remove_back(&vertexStack);
					edgeIndex = (size_t)dym_array_remove_back(&edgeStack);
				}
			}
		} while (ret == ERR_SUCCESS && (v != FirstVertex || edgeIndex < v->degreeOut));

		FirstVertex->Finished = FALSE;
		dym_array_destroy(&edgeStack);
		dym_array_destroy(&vertexStack);
	} else ret = dym_array_push_back(EndArray, FirstVertex);

	return ret;
}


static PKMER_EDGE _get_refseq_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_succ_edge(Vertex, i);

		if (tmp->Type == kmetReference) {
			ret = tmp;
			break;
		}
	}

	return ret;
}


static boolean _has_outgoing_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = FALSE;
	const KMER_EDGE *e = NULL;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		PKMER_EDGE e = kmer_vertex_get_succ_edge(Vertex, i);

		ret = (e->Type == kmetRead);
		if (ret)
			break;
	}

	return ret;
}


ERR_VALUE kmer_graph_resolve_bubbles(PKMER_GRAPH Graph, const uint32_t Threshold)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX bubbleStart = Graph->StartingVertex;
	size_t readEdgeIndex = (size_t)-1;

	ret = ERR_SUCCESS;
	do {		
		PKMER_EDGE rsEdge = _get_refseq_edge(bubbleStart);
		
		for (size_t i = 0; i < bubbleStart->degreeOut; ++i) {
			PKMER_EDGE tmp = kmer_vertex_get_succ_edge(bubbleStart, i);

			if (tmp->Type == kmetRead) {
				DYM_ARRAY bubbleEnds;
				
				readEdgeIndex = i;
				dym_array_create(&bubbleEnds, 140);
				ret = _get_bubble_ends(tmp->Dest, &bubbleEnds);
				if (ret == ERR_SUCCESS) {
					PKMER_VERTEX v = rsEdge->Dest;
					boolean noReadSupport = rsEdge->Weight == Threshold;

					while (!_has_outgoing_read_edges(v) && dym_array_size(&bubbleEnds) > 0 && v != Graph->EndingVertex) {
						size_t index = 0;

						if (dym_array_find(&bubbleEnds, v, &index)) {
							dym_array_remove_by_item_fast(&bubbleEnds, v);
							if (noReadSupport)
								kmer_graph_delete_edge(Graph, rsEdge);

							rsEdge = _get_refseq_edge(v);
							noReadSupport = rsEdge->Weight == Threshold;
							v = rsEdge->Dest;
						} else {
							rsEdge = _get_refseq_edge(v);
							noReadSupport &= (rsEdge->Weight == Threshold);
							v = rsEdge->Dest;
						}
					}
				}

				readEdgeIndex = (size_t)-1;
				dym_array_destroy(&bubbleEnds);
				break;
			}
		}

		if (ret == ERR_SUCCESS) {
			if (readEdgeIndex == (size_t)-1)
				bubbleStart = rsEdge->Dest;
		}
	} while (ret == ERR_SUCCESS && bubbleStart != Graph->EndingVertex);

	size_t dummy = 0;
	kmer_graph_delete_trailing_things(Graph, &dummy);

	return ret;
}


void kmer_graph_compute_shurtcuts(PKMER_GRAPH Graph, const size_t MaxLength)
{
	PKMER_GRAPH_SHORTCUT tmpShortcut = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	char *seqBuffer = NULL;

	ret = utils_calloc(MaxLength + 1, sizeof(char), &seqBuffer);
	if (ret == ERR_SUCCESS) {
		memset(seqBuffer, 0, (MaxLength + 1)*sizeof(char));
		ret = kmer_table_first(Graph->VertexTable, &iter);
		while (ret != ERR_NO_MORE_ENTRIES) {
			PKMER_VERTEX shortcutStart = NULL;
			PKMER_EDGE shortcutEdge = NULL;

			shortcutStart = (PKMER_VERTEX)iter->Data;
			if ((shortcutStart->DegreeIn == 0 && shortcutStart->degreeOut == 1) || (shortcutStart->degreeOut > 1 || shortcutStart->DegreeIn > 1)) {
				for (size_t i = 0; i < shortcutStart->degreeOut; ++i) {
					PKMER_VERTEX u = NULL;

					shortcutEdge = kmer_vertex_get_succ_edge(shortcutStart, i);
					u = shortcutEdge->Dest;
					if (u->DegreeIn == 1 && u->degreeOut == 1) {
						PKMER_EDGE e = shortcutEdge;
						uint32_t sLength = 0;
						uint32_t maxPassCount = (uint32_t)-1;
						size_t seqIndex = 0;

						do {
							sLength += e->Length;
							maxPassCount = min(maxPassCount, e->MaxPassCount);
							seqBuffer[seqIndex] = kmer_get_base(u->KMer, kmer_graph_get_kmer_size(Graph) - 1);
							++seqIndex;
							e = kmer_vertex_get_succ_edge(u, 0);
							u = e->Dest;
						} while (u->DegreeIn == 1 && u->degreeOut == 1);

						sLength += e->Length;
						maxPassCount = min(maxPassCount, e->MaxPassCount);
						seqBuffer[seqIndex] = kmer_get_base(u->KMer, kmer_graph_get_kmer_size(Graph) - 1);
						++seqIndex;
						seqBuffer[seqIndex] = '\0';
						ret = utils_malloc(sizeof(KMER_GRAPH_SHORTCUT), (void **)&tmpShortcut);
						if (ret == ERR_SUCCESS) {
							tmpShortcut->StartVertex = shortcutStart;
							tmpShortcut->EndVertex = u;
							tmpShortcut->Length = sLength;
							tmpShortcut->PassCount = 0;
							tmpShortcut->MaxPassCount = maxPassCount;
							tmpShortcut->Probability = shortcutEdge->Probability;
							ret = utils_copy_string(seqBuffer, &tmpShortcut->Sequence);
							if (ret == ERR_SUCCESS)
								shortcutEdge->Shortcut = tmpShortcut;

							if (ret != ERR_SUCCESS)
								utils_free(tmpShortcut);
						}
					}
				}
			}

			ret = kmer_table_next(Graph->VertexTable, iter, &iter);
		}

		utils_free(seqBuffer);
	}

	return;
}


ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(KMer, Type, &v);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_table_insert(Graph->VertexTable, KMer, v);
			if (ret == ERR_TABLE_FULL) {
				ret = kmer_table_extend(Graph->VertexTable);
				if (ret == ERR_SUCCESS)
					ret = ERR_TABLE_FULL;
			}
		} while (ret == ERR_TABLE_FULL);

		if (ret == ERR_SUCCESS) {
			Graph->NumberOfVertices++;
			if (Vertex != NULL)
				*Vertex = v;
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(v);
	}

	return ret;
}


ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type)
{
	return kmer_graph_add_vertex_ex(Graph, KMer, Type, NULL);
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight, const uint32_t Length, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX v = NULL;
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	u = kmer_table_get(Graph->VertexTable, Source);
	v = kmer_table_get(Graph->VertexTable, Dest);
	if (u != NULL && v != NULL) {
		ret = _edge_create(u, v, Type, weight, Length, &edge);
		if (ret == ERR_SUCCESS) {
			do {
				ret = kmer_edge_table_insert(Graph->EdgeTable, Source, Dest, edge);
				if (ret == ERR_SUCCESS) {
					ret = dym_array_prepare_for_insert(&u->Successors, 1);
					if (ret == ERR_SUCCESS) {
						ret = dym_array_prepare_for_insert(&v->Predecessors, 1);
						if (ret == ERR_SUCCESS) {
							u->degreeOut++;
							dym_array_push_back_no_alloc(&u->Successors, edge);
							v->DegreeIn++;
							dym_array_push_back_no_alloc(&v->Predecessors, edge);
							Graph->NumberOfEdges++;
							*Edge = edge;
						}
					}
				}

				if (ret == ERR_TABLE_FULL) {
					ret = kmer_edge_table_extend(Graph->EdgeTable);
					if (ret == ERR_SUCCESS)
						ret = ERR_TABLE_FULL;
				}
			} while (ret == ERR_TABLE_FULL);

			if (ret != ERR_SUCCESS)
				_edge_destroy(edge);

			if (ret == ERR_ALREADY_EXISTS)
				*Edge = (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source, Dest);
		}
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_graph_add_edge(PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest, const long weight)
{
	PKMER_EDGE dummy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_add_edge_ex(Graph, Source, Dest, weight, 1, kmetReference, &dummy);

	return ret;
}


PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest)
{
	return (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source, Dest);
}


PKMER_VERTEX kmer_graph_get_vertex(const struct _KMER_GRAPH *Graph, const struct _KMER *KMer)
{
	return (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
}


ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PDYM_ARRAY VertexArray)
{
	return kmer_table_get_multiple(Graph->VertexTable, KMer, VertexArray);
}


ERR_VALUE kmer_graph_delete_vertex(PKMER_GRAPH Graph, PKMER_VERTEX Vertex)
{
	PKMER_EDGE inEdge = NULL;
	PKMER_EDGE outEdge = NULL;
	PKMER_VERTEX predVertex = NULL;
	PKMER_VERTEX succVertex = NULL;
	PKMER_EDGE newEdge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Vertex->DegreeIn <= 1 && Vertex->degreeOut <= 1) {
		if (Vertex->degreeOut == 1) {
			outEdge = kmer_vertex_get_succ_edge(Vertex, 0);
			succVertex = outEdge->Dest;
		}

		if (Vertex->DegreeIn == 1) {
			inEdge = kmer_vertex_get_pred_edge(Vertex, 0);
			predVertex = inEdge->Source;
		}

		if (predVertex == NULL || (predVertex != succVertex)) {
			if (inEdge != NULL && outEdge != NULL) {
				ret = kmer_graph_merge_edges(Graph, inEdge, outEdge);
			} else if (inEdge != NULL) {
				kmer_graph_delete_edge(Graph, inEdge);
				ret = ERR_SUCCESS;
			} else if (outEdge != NULL) {
				kmer_graph_delete_edge(Graph, outEdge);
				ret = ERR_SUCCESS;
			} else ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				kmer_table_delete(Graph->VertexTable, Vertex->KMer);
				--Graph->NumberOfVertices;
			}
		} else ret = ERR_PRED_IS_SUCC;
	} else ret = ERR_TOO_MANY_EDGES;

	return ret;
}


void kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge)
{
	PKMER_VERTEX source = Edge->Source;
	PKMER_VERTEX dest = Edge->Dest;

	dym_array_remove_by_item_fast(&source->Successors, Edge);
	dym_array_remove_by_item_fast(&dest->Predecessors, Edge);
	--source->degreeOut;
	--dest->DegreeIn;

	kmer_edge_table_delete(Graph->EdgeTable, source->KMer, dest->KMer);
	--Graph->NumberOfEdges;

	return;
}


ERR_VALUE kmer_graph_merge_edges(PKMER_GRAPH Graph, PKMER_EDGE Source, PKMER_EDGE Dest)
{
	PKMER_VERTEX v = NULL;
	PKMER_VERTEX u = NULL;
	PKMER_VERTEX w = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	u = Source->Source;
	v = Source->Dest;
	w = Dest->Dest;
	if (v == Dest->Source && u != w) {
		if (kmer_graph_get_edge(Graph, u->KMer, w->KMer) == NULL) {
			long weight = max(Source->Weight, Dest->Weight);
			uint32_t Length = Source->Length + Dest->Length;
			uint32_t mpc = Source->MaxPassCount;
			EKMerEdgeType type = (Source->Type == kmetReference && Dest->Type == kmetReference) ? kmetReference : kmetRead;
			PKMER_EDGE newEdge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, u->KMer, w->KMer, weight, Length, type, &newEdge);
			if (ret == ERR_SUCCESS) {
				newEdge->MaxPassCount = mpc;
				newEdge->SeqLen = Source->SeqLen + 1 + Dest->SeqLen;
				ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
				if (ret == ERR_SUCCESS) {
					memcpy(newEdge->Seq, Source->Seq, Source->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
					memcpy(newEdge->Seq + Source->SeqLen + 1, Dest->Seq, Dest->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen + 1 + Dest->SeqLen] = '\0';
					kmer_graph_delete_edge(Graph, Source);
					kmer_graph_delete_edge(Graph, Dest);					
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, newEdge);
			}
		} else ret = ERR_TRIANGLE;
	} else ret = (u == w) ? ERR_PRED_IS_SUCC : ERR_NOT_ADJACENT;

	return ret;
}


ERR_VALUE kmer_graph_get_seqs(PKMER_GRAPH Graph, PDYM_ARRAY SeqArray)
{
	size_t l = 0;
	char s[16384];
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX v = Graph->StartingVertex;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	size_t edgeIndex = 0;
	PKMER_EDGE edge = NULL;
	DYM_ARRAY edgeIndices;
	DYM_ARRAY edges;

	dym_array_create(&edgeIndices, 140);
	dym_array_create(&edges, 140);
	memset(s, 0, sizeof(s));
	memcpy(s, v->KMer->Bases + 1, kmerSize - 1);
	l += (kmerSize - 1);
	do {
		if (edgeIndex == v->degreeOut) {
			v->Finished = FALSE;
			if (v == Graph->StartingVertex)
				break;

			edge = dym_array_remove_back(&edges);
			edgeIndex = (size_t)dym_array_remove_back(&edgeIndices);
			v = edge->Source;
			l -= (edge->SeqLen + 1);
			memset(s + l, 0, (edge->SeqLen + 1)*sizeof(char));
		} else {
			edge = kmer_vertex_get_succ_edge(v, edgeIndex);
			memcpy(s + l, edge->Seq, edge->SeqLen*sizeof(char));
			l += edge->SeqLen;
			if (edge->Dest == Graph->EndingVertex) {
				char *tmp = NULL;

				ret = utils_copy_string(s, &tmp);
				if (ret == ERR_SUCCESS) {
					ret = dym_array_push_back(SeqArray, tmp);
					if (ret == ERR_SUCCESS) {
						l -= edge->SeqLen;
						memset(s + l, 0, edge->SeqLen*sizeof(char));
						++edgeIndex;
					}

					if (ret != ERR_SUCCESS)
						utils_free(tmp);
				}
			} else {
				v = edge->Dest;
				if (!v->Finished) {
					edge->Source->Finished = TRUE;
					s[l] = kmer_get_base(v->KMer, kmerSize - 1);
					++l;
					ret = dym_array_push_back(&edgeIndices, (void *)(edgeIndex + 1));
					if (ret != ERR_SUCCESS)
						break;

					ret = dym_array_push_back(&edges, edge);
					if (ret != ERR_SUCCESS)
						break;

					edgeIndex = 0;
				} else {
					v = edge->Source;
					l -= edge->SeqLen;
					memset(s + l, 0, edge->SeqLen*sizeof(char));
					++edgeIndex;
				}
			}
		}
	} while (ret == ERR_SUCCESS && v != Graph->StartingVertex);

	dym_array_destroy(&edges);
	dym_array_destroy(&edgeIndices);
	if (ret != ERR_SUCCESS) {
		for (size_t i = 0; i < dym_array_size(SeqArray); ++i)
			utils_free(dym_array_get(SeqArray, i));

		dym_array_clear(SeqArray);
	}

	if (ret == ERR_SUCCESS)
		printf("%u\n", (uint32_t)dym_array_size(SeqArray));

	return ret;
}


