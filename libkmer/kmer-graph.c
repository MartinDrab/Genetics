
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

static ERR_VALUE _vertex_create(const PKMER KMer, EKMerVertexType Type, PKMER_VERTEX *Result)
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
			tmp->PassCount = 0;
			tmp->CurrentPass = NULL;
			tmp->Passes = NULL;
			tmp->StartIndex = 0;
			*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmp);
	}

	return ret;
}


static void _vertex_destroy(PKMER_VERTEX Vertex)
{
	if (Vertex->PassCount > 0)
		utils_free(Vertex->Passes);

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
			if (ret == ERR_SUCCESS) {
				tmp->PassCount = Vertex->PassCount;
				if (tmp->PassCount > 0) {
					ret = utils_calloc(tmp->PassCount, sizeof(KMER_VERTEX_PASS), &tmp->Passes);
					if (ret == ERR_SUCCESS) {
						memcpy(tmp->Passes, Vertex->Passes, tmp->PassCount*sizeof(KMER_VERTEX_PASS));
						tmp->CurrentPass = tmp->Passes + (Vertex->CurrentPass - Vertex->Passes);
					}

					if (ret != ERR_SUCCESS)
						tmp->PassCount = 0;
				}

				if (ret == ERR_SUCCESS)
					*Result = tmp;
			}			
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(tmp);
	}

	return ret;
}

/************************************************************************/
/*                        EDGE BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _edge_create(const PKMER_VERTEX Source, const PKMER_VERTEX Dest, EKMerEdgeType Type, const long Weight, const uint32_t Length, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc_KMER_EDGE(&tmp);
	if (ret == ERR_SUCCESS) {
		dym_array_create(&tmp->ReadIndices, 140);
		dym_array_create(&tmp->ReadPositions, 140);
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

	dym_array_destroy(&Edge->ReadPositions);
	dym_array_destroy(&Edge->ReadIndices);
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
			ret = dym_array_copy(&tmp->ReadIndices, &Edge->ReadIndices);
			if (ret == ERR_SUCCESS) {
				ret = dym_array_copy(&tmp->ReadPositions, &Edge->ReadPositions);
				if (ret == ERR_SUCCESS) {
					tmp->SeqLen = Edge->SeqLen;
					tmp->PassCount = Edge->PassCount;
					tmp->MaxPassCount = Edge->MaxPassCount;
					tmp->Probability = Edge->Probability;
					tmp->Order = Edge->Order;
					*Result = tmp;
				}
			}
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


static void _correct_connecing_vertices_passes(PKMER_EDGE SourceOutgoing, PKMER_EDGE DestIncomming, PKMER_EDGE ConnectingEdge)
{
	PKMER_VERTEX_PASS pass = NULL;

	for (size_t k = 0; k < kmer_vertex_get_pass_count(ConnectingEdge->Source); ++k) {
		pass = kmer_vertex_get_pass(ConnectingEdge->Source, k);

		if (pass->Outgoing == SourceOutgoing)
			pass->Outgoing = ConnectingEdge;
	}

	for (size_t k = 0; k < kmer_vertex_get_pass_count(ConnectingEdge->Dest); ++k) {
		pass = kmer_vertex_get_pass(ConnectingEdge->Dest, k);

		if (pass->Incomming == DestIncomming)
			pass->Incomming = ConnectingEdge;
	}

	return;
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


void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph)
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

			if (deleted)
				kmer_graph_delete_vertex(Graph, v);
		}

		ret = (deleted) ?
			kmer_table_first(Graph->VertexTable, &iter) :
			kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
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


ERR_VALUE kmer_graph_add_vertex(PKMER_GRAPH Graph, const PKMER KMer, const EKMerVertexType Type)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(KMer, Type, &v);
	if (ret == ERR_SUCCESS) {
		do {
			ret = kmer_table_insert(Graph->VertexTable, KMer, v);
			if (ret == ERR_TABLE_FULL)
				kmer_table_extend(Graph->VertexTable);
		} while (ret == ERR_TABLE_FULL);

		if (ret == ERR_SUCCESS)
			Graph->NumberOfVertices++;

		if (ret != ERR_SUCCESS)
			_vertex_destroy(v);
	}

	return ret;
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

				if (ret == ERR_TABLE_FULL)
					kmer_edge_table_extend(Graph->EdgeTable);
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



static ERR_VALUE _kmer_graph_capture_refseq(const KMER_GRAPH *Graph, PKMER_VERTEX Start, const KMER_VERTEX *End, char **Seq, size_t *SeqLen)
{
	size_t l = 0;
	char *s = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(16384, sizeof(char), (void **)&s);
	if (ret == ERR_SUCCESS) {
		PKMER_VERTEX v = Start;

		if (v->Type == kmvtRefSeqStart) {
			l = Graph->KMerSize - 2;
			for (size_t i = 0; i < l; i++)
				s[i] = kmer_get_base(v->KMer, i + 1);
		} else {
			l = Graph->KMerSize - 1;
			for (size_t i = 0; i < l; i++)
				s[i] = kmer_get_base(v->KMer, i);
		}

		while (v != End) {
			const struct _KMER_EDGE *edge = v->CurrentPass->Outgoing;

			s[l] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
			++l;
			if (edge->Seq != NULL) {
				memcpy(s + l, edge->Seq, edge->SeqLen*sizeof(char));
				l += edge->SeqLen;
			}

			++v->CurrentPass;
			v = edge->Dest;
		}

		s[l] = '\0';
		*Seq = s;
		*SeqLen = l;
	}

	return ret;
}


ERR_VALUE kmer_graph_get_seq(const KMER_GRAPH *Graph, char **Seq, size_t *SeqLen)
{
	return _kmer_graph_capture_refseq(Graph, Graph->StartingVertex, Graph->EndingVertex, Seq, SeqLen);
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


ERR_VALUE kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge)
{
	size_t index = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX source = Edge->Source;
	PKMER_VERTEX dest = Edge->Dest;
	PKMER_VERTEX_PASS pass = NULL;

	dym_array_remove_by_item_fast(&source->Successors, Edge);
	dym_array_remove_by_item_fast(&dest->Predecessors, Edge);
	--source->degreeOut;
	--dest->DegreeIn;

	index = 0;
	while (index < kmer_vertex_get_pass_count(source)) {
		pass = kmer_vertex_get_pass(source, index);
		if (pass->Outgoing == Edge) {
			kmer_vertex_remove_pass(source, index);
			continue;
		}

		++index;
	}

	index = 0;
	while (index < kmer_vertex_get_pass_count(dest)) {
		pass = kmer_vertex_get_pass(dest, index);
		if (pass->Incomming == Edge) {
			kmer_vertex_remove_pass(dest, index);
			continue;
		}

		++index;
	}

	kmer_edge_table_delete(Graph->EdgeTable, source->KMer, dest->KMer);
	--Graph->NumberOfEdges;

	return ret;
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
					PKMER_VERTEX_PASS pass = NULL;
					size_t passCount = 0;

					memcpy(newEdge->Seq, Source->Seq, Source->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
					memcpy(newEdge->Seq + Source->SeqLen + 1, Dest->Seq, Dest->SeqLen*sizeof(char));
					newEdge->Seq[Source->SeqLen + 1 + Dest->SeqLen] = '\0';

					_correct_connecing_vertices_passes(Source, Dest, newEdge);
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
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char seq[16384];
	char *s = seq;
	size_t l = 0;
	PKMER_VERTEX v = NULL;
	DYM_ARRAY oldStartIndices;
	DYM_ARRAY lastEdges;
	DYM_ARRAY lastPasses;

	dym_array_create(SeqArray, 140);
	dym_array_create(&oldStartIndices, 140);
	dym_array_create(&lastEdges, 140);
	dym_array_create(&lastPasses, 140);
	memset(seq, 0, sizeof(seq));
	ret = ERR_SUCCESS;

	v = Graph->StartingVertex;
	l += (kmer_graph_get_kmer_size(Graph) - 1);
	for (size_t i = 0; i < kmer_graph_get_kmer_size(Graph) - 1; ++i)
		s[i] = kmer_get_base(v->KMer, i + 1);

	while (ret == ERR_SUCCESS && v != Graph->EndingVertex) {
		if ((size_t)(v->CurrentPass - v->Passes) >= v->PassCount) {
			if (v == Graph->StartingVertex)
				break;

			PKMER_EDGE lastEdge = (PKMER_EDGE)dym_array_remove_back(&lastEdges);
			size_t oldStartIndex = (size_t)dym_array_remove_back(&oldStartIndices);
			PKMER_VERTEX_PASS lastPass = (PKMER_VERTEX_PASS)dym_array_remove_back(&lastPasses);

			--l;
			s[l] = '\0';
			l -= lastEdge->SeqLen;
			memset(s + l, 0, lastEdge->SeqLen*sizeof(char));

			v->StartIndex = oldStartIndex;
			v->CurrentPass = lastPass;
			v = lastEdge->Source;
		} else {
			const KMER_EDGE *nextEdge = NULL;
			
			for (size_t i = v->StartIndex; i < v->degreeOut; ++i) {
				nextEdge = kmer_vertex_get_succ_edge(v, i);
				if (nextEdge->Type == kmetRead) {
					v->StartIndex = i + 1;
					break;
				}

				nextEdge = NULL;
			}

			if (nextEdge == NULL) {
				nextEdge = v->CurrentPass->Outgoing;
				++v->CurrentPass;
			}

			ret = dym_array_prepare_for_insert(&lastPasses, 1);
			if (ret == ERR_SUCCESS) {
				ret = dym_array_prepare_for_insert(&lastEdges, 1);
				if (ret == ERR_SUCCESS) {
					ret = dym_array_prepare_for_insert(&oldStartIndices, 1);
					if (ret == ERR_SUCCESS) {
						memcpy(s + l, nextEdge->Seq, nextEdge->SeqLen*sizeof(char));
						l += nextEdge->SeqLen;
						v = nextEdge->Dest;
						if (v == Graph->EndingVertex) {
							char *s2 = NULL;

							ret = utils_calloc(l + 1, sizeof(char), &s2);
							if (ret == ERR_SUCCESS) {
								memcpy(s2, s, l*sizeof(char));
								s2[l] = '\0';
								ret = dym_array_push_back(SeqArray, s2);
								if (ret == ERR_SUCCESS) {
									l -= nextEdge->SeqLen;
									memset(s + l, 0, nextEdge->SeqLen*sizeof(char));

									v = nextEdge->Source;
								}

								if (ret != ERR_SUCCESS)
									utils_free(s2);
							}
						} else  {
							s[l] = kmer_get_base(v->KMer, kmer_graph_get_kmer_size(Graph) - 1);
							++l;
							dym_array_push_back_no_alloc(&oldStartIndices, (void *)v->StartIndex);
							dym_array_push_back_no_alloc(&lastEdges, nextEdge);
							dym_array_push_back_no_alloc(&lastPasses, v->CurrentPass);
						}
					}
				}
			}
		}
	}

	dym_array_destroy(&lastPasses);
	dym_array_destroy(&lastEdges);
	dym_array_destroy(&oldStartIndices);
	if (ret != ERR_SUCCESS) {
		for (size_t i = 0; i < dym_array_size(SeqArray); ++i)
			utils_free(dym_array_get(SeqArray, i));
	
		dym_array_destroy(SeqArray);
	}

	return ret;
}


ERR_VALUE kmer_vertex_add_pass(PKMER_VERTEX Vertex, const struct _KMER_EDGE *Incomming, const struct _KMER_EDGE *Outgoing, const EVertexPassType PassType)
{
	PKMER_VERTEX_PASS newPasses = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(Vertex->PassCount + 1, sizeof(KMER_VERTEX_PASS), &newPasses);
	if (ret == ERR_SUCCESS) {
		memcpy(newPasses, Vertex->Passes, Vertex->PassCount*sizeof(KMER_VERTEX_PASS));
		if (Vertex->PassCount > 0)
			utils_free(Vertex->Passes);
		
		Vertex->CurrentPass = newPasses;
		Vertex->Passes = newPasses;
		newPasses[Vertex->PassCount].Incomming = Incomming;
		newPasses[Vertex->PassCount].Outgoing = Outgoing;
		newPasses[Vertex->PassCount].PassType = PassType;
		++Vertex->PassCount;
	}

	return ret;
}


void kmer_vertex_remove_pass(PKMER_VERTEX Vertex, const size_t Index)
{
	memmove(Vertex->Passes + Index, Vertex->Passes + Index + 1, (Vertex->PassCount - Index - 1)*sizeof(KMER_VERTEX_PASS));
	--Vertex->PassCount;
	if (Vertex->PassCount == 0) {
		utils_free(Vertex->Passes);
		Vertex->Passes = NULL;
	}

	return;
}


void kmer_vertex_remove_passes(PKMER_VERTEX Vertex, const struct _KMER_EDGE *Source, const struct _KMER_EDGE *Dest)
{
	size_t index = 0;
	PKMER_VERTEX_PASS p = NULL;

	while (index < kmer_vertex_get_pass_count(Vertex)) {
		p = kmer_vertex_get_pass(Vertex, index);
		if (p->Incomming == Source && p->Outgoing == Dest) {
			kmer_vertex_remove_pass(Vertex, index);
			continue;
		}

		++index;
	}

	return;
}


ERR_VALUE kmer_edge_add_read(PKMER_EDGE Edge, const size_t ReadIndex, const size_t ReadPosition)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = dym_array_push_back(&Edge->ReadIndices, (void *)ReadIndex);
	if (ret == ERR_SUCCESS) {
		ret = dym_array_push_back(&Edge->ReadPositions, (void *)ReadPosition);
		if (ret != ERR_SUCCESS)
			dym_array_remove_back(&Edge->ReadIndices);
	}

	return ret;
}


ERR_VALUE kmer_edge_connecting_reads(const struct _KMER_EDGE **EdgeSet, const size_t EdgeSetSize, PDYM_ARRAY ReadIndices)
{
	size_t readIndex = 0;
	size_t readIndex2 = 0;
	boolean found = FALSE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const struct _KMER_EDGE *e = NULL;
	const struct _KMER_EDGE *firstEdge = *EdgeSet;

	ret = ERR_SUCCESS;
	if (EdgeSetSize > 0) {
		for (size_t i = 0; i < dym_array_size(&firstEdge->ReadIndices); ++i) {
			readIndex = (size_t)dym_array_get(&firstEdge->ReadIndices, i);
			found = FALSE;
			for (size_t j = 1; j < EdgeSetSize; ++j) {
				e = EdgeSet[j];
				for (size_t k = 0; k < dym_array_size(&e->ReadIndices); ++k) {
					readIndex2 = (size_t)dym_array_get(&e->ReadIndices, k);
					found = (readIndex2 == readIndex);
					if (found || readIndex2 > readIndex)
						break;
				}

				if (!found)
					break;
			}

			if (found) {
				ret = dym_array_push_back(ReadIndices, (void *)readIndex);
				if (ret != ERR_SUCCESS)
					break;
			}
		}
	}

	return ret;
}


void kmer_edge_get_read(const KMER_EDGE *Edge, const size_t Index, size_t *ReadIndex, size_t *ReadPosition)
{
	assert(Index < kmer_edge_get_read_count(Edge));
	*ReadIndex = (size_t)dym_array_get(&Edge->ReadIndices, Index);
	*ReadPosition = (size_t)dym_array_get(&Edge->ReadPositions, Index);

	return;
}


ERR_VALUE kmer_edge_find_read(const KMER_EDGE *Edge, const size_t StartIndex, size_t *EndIndex, const size_t ReadIndex, size_t *ReadPosition)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_NOT_FOUND;
	for (size_t i = StartIndex; i < kmer_edge_get_read_count(Edge); ++i) {
		if ((size_t)dym_array_get(&Edge->ReadIndices, i) == ReadIndex) {
			if (EndIndex != NULL)
				*EndIndex = i + 1;

			*ReadPosition = (size_t)dym_array_get(&Edge->ReadPositions, i);
			ret = ERR_SUCCESS;
			break;
		}
	}

	return ret;
}


void kmer_edge_remove_read(PKMER_EDGE Edge, const size_t StartIndex, const size_t ReadIndex)
{
	size_t dummy = 0;
	size_t index = 0;

	while (kmer_edge_find_read(Edge, StartIndex, &index, ReadIndex, &dummy) == ERR_SUCCESS) {
		dym_array_remove(&Edge->ReadIndices, index - 1);
		dym_array_remove(&Edge->ReadPositions, index - 1);
	}

	return;
}


void kmer_graph_separate_distinct_passes(PKMER_GRAPH Graph)
{
	PKMER_VERTEX_PASS p = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_table_first(Graph->VertexTable, &iter);
	while (err == ERR_SUCCESS) {
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;
		boolean found = FALSE;
			
		do {
			size_t count = kmer_vertex_get_pass_count(v);
				
			found = FALSE;
			if (v->degreeOut > 1 || v->DegreeIn > 1) {
				for (size_t i = 0; i < count; ++i) {
					p = kmer_vertex_get_pass(v, i);
					for (size_t j = 0; j < count; ++j) {
						PKMER_VERTEX_PASS p2 = kmer_vertex_get_pass(v, j);

						if ((p->Incomming == p2->Incomming && p->Outgoing != p2->Outgoing) ||
							(p->Incomming != p2->Incomming && p->Outgoing == p2->Outgoing))
							break;

						if (j == count - 1)
							found = (kmer_graph_merge_edges(Graph, (PKMER_EDGE)p->Incomming, (PKMER_EDGE)p->Outgoing)) == ERR_SUCCESS;
					}

					if (found)
						break;
				}
			}
		} while (found);

		err = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}


static void _solve_one_repetition(PKMER_GRAPH Graph, PKMER_VERTEX StartVertex, PKMER_VERTEX EndVertex, const size_t PassIndex)
{


	return;
}


void kmer_graph_resolve_repetitions(PKMER_GRAPH Graph)
{
	PKMER_VERTEX_PASS p = NULL;
	PKMER_TABLE_ENTRY iter = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_table_first(Graph->VertexTable, &iter);
	while (err == ERR_SUCCESS) {
		PKMER_VERTEX v = (PKMER_VERTEX)iter->Data;

		if (v->DegreeIn > 1 && v->degreeOut == 1) {
			PKMER_VERTEX endV = kmer_vertex_get_successor(v, 0);

			while (endV->degreeOut == 1)
				endV = kmer_vertex_get_successor(endV, 0);

			if (endV->Type == kmvtRefSeqMiddle) {
				boolean found = FALSE;

				do {
					PKMER_VERTEX_PASS inPass = v->Passes;
					uint32_t inPassCount = v->PassCount;
					PKMER_VERTEX_PASS outPass = endV->Passes;
					uint32_t outPassCount = endV->PassCount;
					PKMER_VERTEX_PASS inp = inPass;

					found = FALSE;
					if (inPassCount == outPassCount) {
						for (size_t i = 0; i < inPassCount; ++i) {
							PKMER_VERTEX_PASS outp = outPass;

							for (size_t j = 0; j < outPassCount; ++j) {

								if ((inp->Incomming == outp->Incomming && inp->Outgoing != outp->Outgoing) ||
									(inp->Incomming != outp->Incomming && inp->Outgoing == outp->Outgoing))
									break;

								if (j == outPassCount - 1) {
									_solve_one_repetition(Graph, v, endV, i);
									found = TRUE;
								}

								++outp;
							}

							++inp;
						}
					}
				} while (found);
			}
		}

		err = kmer_table_next(Graph->VertexTable, iter, &iter);
	}

	return;
}

static void _connecting_read_get_seq(const KMER_GRAPH *Graph, const ONE_READ *Read, const size_t ReadIndex, const KMER_EDGE *StartingEdge, const KMER_EDGE *EndingEdge, char **StartPos, char **EndPos)
{
	size_t pos = 0;
	ERR_VALUE err = ERR_INTERNAL_ERROR;

	err = kmer_edge_find_read(StartingEdge, 0, NULL, ReadIndex, &pos);
	assert(err == ERR_SUCCESS);
	*StartPos = Read->ReadSequence + pos + kmer_graph_get_kmer_size(Graph);
	err = kmer_edge_find_read(EndingEdge, 0, NULL, ReadIndex, &pos) + kmer_graph_get_kmer_size(Graph);
	assert(err == ERR_SUCCESS);
	*EndPos = Read->ReadSequence + pos;

	return;
}


ERR_VALUE kmer_graph_connect_reads(PKMER_GRAPH Graph, const size_t Threshold, const ONE_READ *Reads)
{
	boolean connectingPerformed = FALSE;
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE_ENTRY iter = NULL;
	DYM_ARRAY divergenceStarts;
	DYM_ARRAY divergenceEnds;

	dym_array_create(&divergenceStarts, 140);
	dym_array_create(&divergenceEnds, 140);
	ret = kmer_edge_table_first(Graph->EdgeTable, &iter);
	while (ret == ERR_SUCCESS) {
		edge = (PKMER_EDGE)iter->Data;
		if (edge->Source->Type == kmvtRefSeqMiddle && edge->Dest->Type == kmvtRead)
			ret = dym_array_push_back(&divergenceStarts, edge);
		else if (edge->Source->Type == kmvtRead && edge->Dest->Type == kmvtRefSeqMiddle)
			ret = dym_array_push_back(&divergenceEnds, edge);

		if (ret == ERR_SUCCESS)
			ret = kmer_edge_table_next(Graph->EdgeTable, iter, &iter);
	}

	if (ret == ERR_NO_MORE_ENTRIES)
		ret = ERR_SUCCESS;

	if (ret == ERR_SUCCESS) {
		DYM_ARRAY connectingReads;
		PKMER_EDGE endingEdge = NULL;
		PKMER_EDGE startingEdge = NULL;
		const KMER_EDGE *examinedEdges[2];

		do {
			connectingPerformed = FALSE;
			for (size_t i = 0; i < dym_array_size(&divergenceEnds); ++i) {
				endingEdge = (PKMER_EDGE)dym_array_get(&divergenceEnds, i);
				for (size_t j = 0; j < dym_array_size(&divergenceStarts); ++j) {
					startingEdge = (PKMER_EDGE)dym_array_get(&divergenceStarts, j);
					dym_array_create(&connectingReads, 140);
					examinedEdges[0] = endingEdge;
					examinedEdges[1] = startingEdge;
					ret = kmer_edge_connecting_reads(examinedEdges, sizeof(examinedEdges) / sizeof(PKMER_EDGE), &connectingReads);
					if (ret == ERR_SUCCESS) {
						if (dym_array_size(&connectingReads) >= ((Threshold > 0) ? Threshold : 1)) {
							PKMER_EDGE newEdge = NULL;
							size_t readIndex = (size_t)dym_array_get(&connectingReads, 0);
							const ONE_READ *read = Reads + readIndex;
							const char *startingPos = NULL;
							const char *endingPos = NULL;

							_connecting_read_get_seq(Graph, read, readIndex, endingEdge, startingEdge, &startingPos, &endingPos);
							if (startingPos < endingPos) {
								ret = kmer_graph_add_edge_ex(Graph, endingEdge->Source->KMer, startingEdge->Dest->KMer, dym_array_size(&connectingReads), 1, kmetRead, &newEdge);
								if (ret == ERR_SUCCESS) {
									newEdge->Length = (endingPos - startingPos) + 1;
									newEdge->MaxPassCount = 0;
									newEdge->SeqLen = edge->Length - 1;
									ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
									if (ret == ERR_SUCCESS) {
										memcpy(newEdge->Seq, startingPos, newEdge->SeqLen);
										newEdge->Seq[newEdge->SeqLen / sizeof(char)] = '\0';
										for (size_t k = 0; k < dym_array_size(&connectingReads); ++k) {
											// TODO: Add only reads that have the same sequence
											// The given vertices may be contained also by reads that follow different paths
											// between them (the vertices).
											ret = kmer_edge_add_read(newEdge, (size_t)dym_array_get(&connectingReads, k), (size_t)-1);
											if (ret != ERR_SUCCESS)
												break;
										}

										if (ret == ERR_SUCCESS) {
											// Since we are going to create a detour for certain reads, it is required to
											// remove them from edges we are going to bypass.
											for (size_t k = 0; k < kmer_edge_get_read_count(newEdge); ++k) {
												size_t index = 0;
												size_t pos = 0;
												const ONE_READ *r = NULL;
												
												kmer_edge_get_read(newEdge, k, &index, &pos);
												r = Reads + index;
												const char *start = NULL;
												const char *stop = NULL;

												_connecting_read_get_seq(Graph, r, index, endingEdge, startingEdge, &start, &stop);
												start -= kmer_graph_get_kmer_size(Graph);
												stop -= (kmer_graph_get_kmer_size(Graph) - 1);
												if (start < stop) {
													PKMER sourceKMer = NULL;
													PKMER destKMer = NULL;

													KMER_STACK_ALLOC(sourceKMer, kmer_graph_get_kmer_size(Graph), start);
													KMER_STACK_ALLOC(destKMer, kmer_graph_get_kmer_size(Graph), start + 1);
													while (start != stop) {
														edge = kmer_graph_get_edge(Graph, sourceKMer, destKMer);
														if (edge != NULL)
															kmer_edge_remove_read(edge, 0, index);

														kmer_advance(sourceKMer, start[kmer_graph_get_kmer_size(Graph)]);
														kmer_advance(destKMer, start[kmer_graph_get_kmer_size(Graph) + 1]);
														++start;
													}
												}
											}

											_correct_connecing_vertices_passes(endingEdge, startingEdge, newEdge);
											kmer_graph_delete_edge(Graph, endingEdge);
											kmer_graph_delete_edge(Graph, startingEdge);
											dym_array_remove_by_item_fast(&divergenceEnds, endingEdge);
											dym_array_remove_by_item_fast(&divergenceStarts, startingEdge);
											connectingPerformed = TRUE;
										}
									}

									if (ret != ERR_SUCCESS)
										kmer_graph_delete_edge(Graph, newEdge);
								}

								if (ret == ERR_ALREADY_EXISTS)
									ret = ERR_SUCCESS;
							}
						}
					}

					dym_array_destroy(&connectingReads);
					if (ret != ERR_SUCCESS || connectingPerformed)
						break;
				}

				if (ret != ERR_SUCCESS || connectingPerformed)
					break;
			}
		} while (ret == ERR_SUCCESS && connectingPerformed);
	}

	dym_array_destroy(&divergenceEnds);
	dym_array_destroy(&divergenceStarts);

	return ret;
}


static boolean _kmer_vertex_no_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = FALSE;
	const size_t edgeCount = Vertex->degreeOut;

	for (size_t i = 0; i < Vertex->degreeOut; ++i) {
		ret = (kmer_vertex_get_succ_edge(Vertex, i))->Type != kmetRead;
		if (!ret)
			break;
	}

	return ret;
}


static void _kmer_vertex_remove_previous_passes(PKMER_VERTEX Vertex)
{
	size_t previousCount = (Vertex->PassCount > 0) ? (Vertex->CurrentPass - Vertex->Passes) : 0;

	if (previousCount > 0) {
		memmove(Vertex->Passes, Vertex->CurrentPass, (Vertex->PassCount - previousCount)*sizeof(KMER_VERTEX_PASS));
		Vertex->CurrentPass -= previousCount;
		Vertex->PassCount -= previousCount;
	}

	return;
}

ERR_VALUE kmer_graph_merge_unbranched_refseq(PKMER_GRAPH Graph)
{
	size_t l = 0;
	char *s = NULL;
	boolean noReads = TRUE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX u = Graph->StartingVertex;
	PKMER_VERTEX v = NULL;

	do {
		long weight = 0;
		uint32_t length = 0;
		DYM_ARRAY vertexArray;

		dym_array_create(&vertexArray, 140);
		v = u->CurrentPass->Outgoing->Dest;
		++u->CurrentPass;
		while (ret == ERR_SUCCESS && noReads && v != Graph->EndingVertex) {
			noReads = _kmer_vertex_no_read_edges(v);
			if (noReads) {
				const KMER_EDGE *tmp = v->CurrentPass->Outgoing;

				weight = max(tmp->Weight, weight);
				ret = dym_array_push_back(&vertexArray, v);
				if (ret == ERR_SUCCESS) {
					++v->CurrentPass;
					v = tmp->Dest;
				}
			}
		}

		if (ret == ERR_SUCCESS) {
			--u->CurrentPass;
			for (size_t i = 0; i < dym_array_size(&vertexArray); ++i)
				--((PKMER_VERTEX)dym_array_get(&vertexArray, i))->CurrentPass;

			ret = _kmer_graph_capture_refseq(Graph, u, v, &s, &l);
			--u->CurrentPass;
			for (size_t i = 0; i < dym_array_size(&vertexArray); ++i)
				--((PKMER_VERTEX)dym_array_get(&vertexArray, i))->CurrentPass;

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE newEdge = NULL;

				ret = kmer_graph_add_edge_ex(Graph, u->KMer, v->KMer, weight, length, kmetReference, &newEdge);
				if (ret == ERR_SUCCESS) {
					const char *seqStart = s + ((u->Type == kmvtRefSeqStart) ? (kmer_graph_get_kmer_size(Graph) - 1) : (kmer_graph_get_kmer_size(Graph)));
					const char *seqEnd = s + l - ((v->Type == kmvtRefSeqEnd) ? 0 : 1);

					newEdge->SeqLen = seqEnd - seqStart;
					ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
					if (ret == ERR_SUCCESS) {
						memcpy(newEdge->Seq, seqStart, newEdge->SeqLen*sizeof(char));
						newEdge->Seq[newEdge->SeqLen / sizeof(char)] = '\0';
						u->CurrentPass->Outgoing = newEdge;
						v->CurrentPass->Incomming = newEdge;
						for (size_t i = 0; i < dym_array_size(&vertexArray); ++i)
							_kmer_vertex_remove_previous_passes((PKMER_VERTEX)dym_array_get(&vertexArray, i));
					}
				}

				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

				if (ret == ERR_SUCCESS)
					u = v;

				utils_free(s);
			}
		}

		dym_array_destroy(&vertexArray);
	} while (ret == ERR_SUCCESS && u != Graph->EndingVertex);

	return ret;
}
