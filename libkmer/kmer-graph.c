
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
			utils_free(tmp);
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
							if (v->Order <= u->Order)
								Graph->NumberOfBackwardEdges++;

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


PKMER_EDGE kmer_graph_get_edge(const PKMER_GRAPH Graph, const PKMER Source, const PKMER Dest)
{
	return (PKMER_EDGE)kmer_edge_table_get_data(Graph->EdgeTable, Source, Dest);
}


ERR_VALUE kmer_graph_get_seq(const KMER_GRAPH *Graph, char **Seq, size_t *SeqLen)
{
	size_t l = 0;
	char *s = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(16384, sizeof(char), (void **)&s);
	if (ret == ERR_SUCCESS) {
		PKMER_VERTEX v = Graph->StartingVertex;

		l = Graph->KMerSize - 2;
		for (size_t i = 0; i < l; i++)
			s[i] = kmer_get_base(v->KMer, i + 1);
	
		while (v != Graph->EndingVertex) {
			if (v->PassCount > 0) {
				PKMER_VERTEX old = v;
				PKMER_EDGE edge = v->CurrentPass->Outgoing;

				s[l] = kmer_get_base(v->KMer, Graph->KMerSize - 1);
				++l;
				if (edge->Seq != NULL) {
					memcpy(s + l, edge->Seq, edge->SeqLen*sizeof(char));
					l += edge->SeqLen;
				}

				v = edge->Dest;
				old->CurrentPass++;
			} else break;
		}

		s[l] = '\0';
		*Seq = s;
		*SeqLen = l;
	}

	return ret;
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
				if (kmer_graph_get_edge(Graph, predVertex->KMer, succVertex->KMer) == NULL) {
					uint32_t length = inEdge->Length + outEdge->Length;
					long weight = min(inEdge->Weight, outEdge->Weight);
					EKMerEdgeType edgeType;

					if (inEdge->Type == kmetReference && outEdge->Type == kmetReference)
						edgeType = kmetReference;
					else edgeType = kmetRead;

					ret = kmer_graph_add_edge_ex(Graph, predVertex->KMer, succVertex->KMer, weight, length, edgeType, &newEdge);
					if (ret == ERR_SUCCESS) {
						newEdge->SeqLen = inEdge->SeqLen + 1 + outEdge->SeqLen;
						ret = utils_calloc(newEdge->SeqLen + 1, sizeof(char), &newEdge->Seq);
						if (ret == ERR_SUCCESS) {
							size_t passCount = 0;

							newEdge->MaxPassCount = inEdge->MaxPassCount;
							memcpy(newEdge->Seq, inEdge->Seq, inEdge->SeqLen*sizeof(char));
							newEdge->Seq[inEdge->SeqLen] = kmer_get_base(Vertex->KMer, Graph->KMerSize - 1);
							memcpy(newEdge->Seq + inEdge->SeqLen + 1, outEdge->Seq, outEdge->SeqLen*sizeof(char));
							newEdge->Seq[newEdge->SeqLen] = '\0';

							passCount = kmer_vertex_get_pass_count(predVertex);
							for (size_t i = 0; i < passCount; ++i) {
								PKMER_VERTEX_PASS pass = kmer_vertex_get_pass(predVertex, i);

								if (pass->Outgoing == inEdge)
									pass->Outgoing = newEdge;
							}

							passCount = kmer_vertex_get_pass_count(succVertex);
							for (size_t i = 0; i < passCount; ++i) {
								PKMER_VERTEX_PASS pass = kmer_vertex_get_pass(succVertex, i);

								if (pass->Incomming == outEdge)
									pass->Incomming = newEdge;
							}

							kmer_graph_delete_edge(Graph, inEdge);
							kmer_graph_delete_edge(Graph, outEdge);
						}

						if (ret != ERR_SUCCESS)
							kmer_graph_delete_edge(Graph, newEdge);
					}
				} else ret = ERR_TRIANGLE;
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


ERR_VALUE kmer_vertex_add_pass(PKMER_VERTEX Vertex, const struct _KMER_EDGE *Incomming, const struct _KMER_EDGE *Outgoing)
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
