
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "err.h"
#include "utils.h"
#include "khash.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "reads.h"
#include "found-sequence.h"
#include "refseq-storage.h"
#include "pointer_array.h"
#include "gen_dym_array.h"
#include "ssw.h"
#include "kmer-graph.h"


/************************************************************************/
/*                      HELPER FUNCTIONS                                */
/************************************************************************/


UTILS_TYPED_MALLOC_FUNCTION(KMER_GRAPH)
UTILS_TYPED_MALLOC_FUNCTION(KMER_VERTEX)
UTILS_TYPED_MALLOC_FUNCTION(KMER_EDGE)

/************************************************************************/
/*                        VERTEX BASIC ROUTINES                         */
/************************************************************************/

static PKMER_VERTEX _default_vertex_allocator(struct _KMER_GRAPH *Graph, void *Context)
{
	PKMER_VERTEX ret = NULL;

	utils_malloc(sizeof(KMER_VERTEX) + kmer_graph_get_kmer_size(Graph)*sizeof(char), &ret);

	return ret;
}


static void _default_vertex_freer(struct _KMER_GRAPH *Graph, PKMER_VERTEX Vertex, void *Context)
{
	utils_free(Vertex);

	return;
}


static PKMER_EDGE _default_edge_allocator(struct _KMER_GRAPH *Graph, void *Context)
{
	PKMER_EDGE ret = NULL;

	utils_malloc(sizeof(KMER_EDGE), &ret);

	return ret;
}


static void _default_edge_freer(struct _KMER_GRAPH *Graph, PKMER_EDGE Edge, void *Context)
{
	utils_free(Edge);

	return;
}


static ERR_VALUE _vertex_create(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	tmp = Graph->Allocator.VertexAllocator(Graph, Graph->Allocator.VertexAllocatorContext);
	if (tmp != NULL) {
		tmp->Unique = TRUE;
		kmer_init_from_kmer(&tmp->KMer, KMer);;
		tmp->Type = Type;
		tmp->LongEdgeAllowed = FALSE;
		tmp->ReadStartAllowed = TRUE;
		tmp->Helper = FALSE;
		tmp->ShortVariant = FALSE;
		pointer_array_init_KMER_EDGE(&tmp->Successors, 140);
		pointer_array_init_KMER_EDGE(&tmp->Predecessors, 140);
		tmp->RefSeqPosition = 0;
		tmp->AbsPos = 0;
		tmp->Lists.Next = NULL;
		tmp->Lists.Graph = NULL;
		*Result = tmp;
		ret = ERR_SUCCESS;
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}


static void _vertex_destroy(PKMER_GRAPH Graph, PKMER_VERTEX Vertex)
{
	PKMER_GRAPH g = Vertex->Lists.Graph;

	if (g != NULL) {
		Vertex->Lists.Next = g->VerticesToDeleteList;
		g->VerticesToDeleteList = Vertex;
	} else {
		pointer_array_finit_KMER_EDGE(&Vertex->Predecessors);
		pointer_array_finit_KMER_EDGE(&Vertex->Successors);
		Graph->Allocator.VertexFreer(Graph, Vertex, Graph->Allocator.VertexAllocatorContext);
	}

	return;
}


static ERR_VALUE _vertex_copy(PKMER_GRAPH Graph, const KMER_VERTEX *Vertex, PKMER_VERTEX *Result)
{
	PKMER_VERTEX tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _vertex_create(Graph, &Vertex->KMer, Vertex->Type, &tmp);
	if (ret == ERR_SUCCESS) {
		tmp->Helper = Vertex->Helper;
		tmp->Order = Vertex->Order;
		tmp->Lists.Graph = Vertex->Lists.Graph;
		pointer_array_init_KMER_EDGE(&tmp->Successors, 140);
		ret = pointer_array_clean_copy_KMER_EDGE(&tmp->Successors, &Vertex->Successors);
		if (ret == ERR_SUCCESS) {
			pointer_array_init_KMER_EDGE(&tmp->Predecessors, 140);
			ret = pointer_array_clean_copy_KMER_EDGE(&tmp->Predecessors, &Vertex->Predecessors);
			if (ret == ERR_SUCCESS)
				*Result = tmp;
		}

		if (ret != ERR_SUCCESS)
			_vertex_destroy(Graph, tmp);
	}

	return ret;
}

/************************************************************************/
/*                        EDGE BASIC ROUTINES                         */
/************************************************************************/

static ERR_VALUE _edge_create(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_EDGE tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

		tmp = Graph->Allocator.EdgeAllocator(Graph, Graph->Allocator.EdgeAllocatorContext);
		if (tmp != NULL) {
			tmp->Source = Source;
			tmp->Dest = Dest;
			tmp->Type = Type;
			tmp->Order = 0;
			tmp->Seq = NULL;
			tmp->SeqLen = 0;
			tmp->Seq1Weight = 0;
			tmp->SeqType = kmetNone;
			read_info_init(&tmp->ReadInfo);
			tmp->MarkedForDelete = FALSE;
			pointer_array_init_VARIANT_CALL(&tmp->VCs, 140);
			tmp->LongData.LongEdge = FALSE;
			tmp->LongData.RefSeqEnd = 0;
			tmp->LongData.RefSeqStart = 0;
			dym_array_init_size_t(&tmp->Weights, 140);
			pointer_array_init_READ_INFO(&tmp->ReadIndices, 140);
			*Edge = tmp;
			ret = ERR_SUCCESS;
		} else ret = ERR_OUT_OF_MEMORY;
	
	return ret;
}


static void _edge_destroy(PKMER_GRAPH Graph, PKMER_EDGE Edge)
{
	for (size_t i = 0; i < pointer_array_size(&Edge->ReadIndices); ++i) {
		PREAD_INFO ri = Edge->ReadIndices.Data[i];

		read_info_finit(ri);
		utils_free(ri);
	}

	pointer_array_finit_READ_INFO(&Edge->ReadIndices);
	dym_array_finit_size_t(&Edge->Weights);
	pointer_array_finit_VARIANT_CALL(&Edge->VCs);
	read_info_finit(&Edge->ReadInfo);
	if (Edge->Seq != NULL)
		utils_free((void *)Edge->Seq);

	Graph->Allocator.EdgeFreer(Graph, Edge, Graph->Allocator.EdgeAllocatorContext);

	return;
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

static void _vertex_table_on_delete(struct _KMER_TABLE *Table, void *ItemData, void *Context)
{
	PKMER_GRAPH g = (PKMER_GRAPH)Context;
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;

	_vertex_destroy(g, v);

	return;
}


static ERR_VALUE _vertex_table_on_copy(struct _KMER_TABLE *Table, void *ItemData, void **Copy, void *Context)
{
	PKMER_GRAPH g = (PKMER_GRAPH)Context;
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;

	return _vertex_copy(g, v, (PKMER_VERTEX *)Copy);
}


static void _vertex_table_on_print(struct _KMER_TABLE *Table, void *ItemData, FILE *Stream)
{
	char *colors[] = { "yellow", "lightgreen", "blue", "red", "white" };
	PKMER_VERTEX v = (PKMER_VERTEX)ItemData;

	fprintf(Stream, "\t");
	kmer_print(Stream, &v->KMer);
	fprintf(Stream, "[label=\"");
	if (!v->Helper)
		kmer_print(Stream, &v->KMer);
	else fprintf(Stream, "Helper#%u", kmer_get_number(&v->KMer));
	
	if (v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart ||
		v->Type == kmvtRefSeqEnd) {
		fprintf(Stream, "\\nRef");
		fprintf(Stream, "\\nPOS: %" PRId64 " (%u)", v->AbsPos + 1, v->RefSeqPosition);
		fprintf(Stream, "\\n%s", v->Unique ? "unique" : "repeat");
	} else fprintf(Stream, "\\nRead");

	fprintf(Stream, "\",style=filled,color=%s]", colors[v->Type]);
	fprintf(Stream, ";\n");

	return;
}

/************************************************************************/
/*                KMER LIST TABLE                                       */
/************************************************************************/

static void _kmerlist_table_on_insert(struct _KMER_TABLE *Table, void *ItemData, const uint32_t Order)
{
	return;
}

static void _kmerlist_table_on_delete(struct _KMER_TABLE *Table, void *ItemData, void *Context)
{
	PKMER_LIST l = (PKMER_LIST)ItemData;

	pointer_array_finit_KMER_VERTEX(&l->Vertices);
	utils_free(l);

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

static void _edge_table_on_delete(struct _KMER_EDGE_TABLE *Table, void *ItemData, void *Context)
{
	PKMER_GRAPH g = (PKMER_GRAPH)Context;
	PKMER_EDGE e = (PKMER_EDGE)ItemData;

	_edge_destroy(g, e);

	return;
}


static void _edge_table_on_print(struct _KMER_EDGE_TABLE *Table, void *ItemData, void *Context, FILE *Stream)
{
	PKMER_EDGE e = (PKMER_EDGE)ItemData;
	const KMER_GRAPH *g = (const KMER_GRAPH *)Context;
	const char *typeStrings[] = {
		"Reference",
		"Read",
		"Variant",
		"None",
		"Max",
	};
	/*
	fprintf(Stream, "\t/**Edge(\n");
	fprintf(Stream, "\t\tSource(");
	kmer_print(Stream, &e->Source->KMer);
	fprintf(Stream, ")\n");
	fprintf(Stream, "\t\tDest(");
	kmer_print(Stream, &e->Dest->KMer);
	fprintf(Stream, ")\n");
	fprintf(Stream, "\t\tType(%s)\n", typeStrings[e->Type]);
	fprintf(Stream, "\t\tReads(\n");
	for (size_t i = 0; i < read_info_get_count(&e->ReadInfo); ++i) {
		const READ_INFO_ENTRY *entry = read_info_get_entry(&e->ReadInfo, i);
	
		fprintf(Stream, "\t\t\t%zu:%zu:%u", entry->ReadIndex, entry->ReadPosition, entry->Quality);
		if (i != read_info_get_count(&e->ReadInfo) - 1)
			fputs(",\n", Stream);
	}

	fprintf(Stream, "\t\t)\n");
*/
//	fprintf(Stream, "\t)**/\n");
	
	fprintf(Stream, "\t");
	kmer_print(Stream, &e->Source->KMer);
	fprintf(Stream, " -> ");
	kmer_print(Stream, &e->Dest->KMer);
	fprintf(Stream, " [color=");
	switch (e->Type) {
		case kmetReference:
			fprintf(Stream, "green");
			fprintf(Stream, ",label=\"W: %Iu (%Iu); L: %Iu", e->Seq1Weight, read_info_get_count(&e->ReadInfo), e->SeqLen);
			if (e->LongData.LongEdge || e->SeqLen == 0) {
				fprintf(Stream, "\\n(");
				for (size_t i = 0; i < read_info_get_count(&e->ReadInfo); ++i)
					fprintf(Stream, "%zu:%zu:%u ", e->ReadInfo.Array.Data[i].ReadIndex, e->ReadInfo.Array.Data[i].ReadPosition, e->ReadInfo.Array.Data[i].Quality);
			
				fprintf(Stream, ")");;
			}

			fprintf(Stream, "\"");
			break;
		case kmetRead:
			if (!e->LongData.LongEdge) {
				fprintf(Stream, "red");
				fprintf(Stream, ",label=\"W: %Iu (%Iu); L: %Iu", e->Seq1Weight, read_info_get_count(&e->ReadInfo), e->SeqLen);
			} else {
				fprintf(Stream, "purple");
				fprintf(Stream, ",label=\"W: %Iu (%Iu); L: %Iu\n%u-%u", e->Seq1Weight, read_info_get_count(&e->ReadInfo), e->SeqLen, e->LongData.RefSeqEnd, e->LongData.RefSeqStart);
			}

			if (e->LongData.LongEdge || e->SeqLen == 0) {
				fprintf(Stream, "\\n(");
				for (size_t i = 0; i < read_info_get_count(&e->ReadInfo); ++i)
					fprintf(Stream, "%zu:%zu:%u ", e->ReadInfo.Array.Data[i].ReadIndex, e->ReadInfo.Array.Data[i].ReadPosition, e->ReadInfo.Array.Data[i].Quality);

				fprintf(Stream, ")");;
			}

			fprintf(Stream, "\"");
			break;
		case kmetVariant: {
			const VARIANT_CALL **var = e->VCs.Data;

			fprintf(Stream, "blue");
			fprintf(Stream, ",label=\"W: %Iu; L: %Iu\\n", e->Seq1Weight, e->SeqLen);
			for (size_t i = 0; i < pointer_array_size(&e->VCs); ++i) {
				fprintf(Stream, "%" PRId64 "\\t%s\\t%s\\t%Iu-%Iu\\n", (uint64_t)(*var)->Pos, (*var)->Ref, (*var)->Alt, (*var)->RefWeight, (*var)->AltWeight);
				++var;
			}

			if (e->LongData.LongEdge || e->SeqLen == 0) {
				fprintf(Stream, "\\n(");
				for (size_t i = 0; i < read_info_get_count(&e->ReadInfo); ++i)
					fprintf(Stream, "%zu:%u ", e->ReadInfo.Array.Data[i].ReadIndex, e->ReadInfo.Array.Data[i].Quality);

				fprintf(Stream, ")");;
			}

			fprintf(Stream, "\"");
		} break;
	}

	fprintf(Stream, "];\n");

	return;
}


static boolean _kmer_vertex_no_read_edges(const KMER_VERTEX *Vertex)
{
	boolean ret = TRUE;

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		ret = (kmer_vertex_get_succ_edge(Vertex, i))->Type != kmetRead;
		if (!ret)
			break;
	}

	if (ret) {
		for (size_t i = 0; i < kmer_vertex_in_degree(Vertex); ++i) {
			ret = (kmer_vertex_get_pred_edge(Vertex, i))->Type != kmetRead;
			if (!ret)
				break;
		}
	}

	return ret;
}


static PKMER_EDGE _get_in_refseq_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < kmer_vertex_in_degree(Vertex); ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_pred_edge(Vertex, i);

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

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		PKMER_EDGE e = kmer_vertex_get_succ_edge(Vertex, i);

		ret = (e->Type == kmetRead);
		if (ret)
			break;
	}

	return ret;
}


static void _init_quality_table(uint8_t *Table)
{
	memset(Table, 100, 256 * sizeof(uint8_t));
	Table[0] = 0;
	memset(Table + 1, 0, 9 * sizeof(char));
	memset(Table + 10, 25, 7 * sizeof(char));
	memset(Table + 17, 50, 8 * sizeof(char));
	memset(Table + 25, 75, 5 * sizeof(char));

	return;
}


static void kmer_graph_check(const KMER_GRAPH *Graph)
{
	void *it = NULL;
	PKMER_LIST l = NULL;

	if (kmer_table_first(Graph->KmerListTable, &it, &l) == ERR_SUCCESS) {
		do {
			for (size_t i = 0; i < pointer_array_size(&l->Vertices); ++i) {
				assert(kmer_seq_equal(&l->Kmer, &(l->Vertices.Data[0]->KMer)));
			}

		} while (kmer_table_next(Graph->KmerListTable, it, &it, &l) == ERR_SUCCESS);
	}


	return;
}



/************************************************************************/
/*                     PUBLIC FUNCTIONS                                 */
/************************************************************************/


PKMER_EDGE _get_typed_edge(const KMER_VERTEX *Vertex, EKMerEdgeType Type)
{
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < kmer_vertex_out_degree(Vertex); ++i) {
		PKMER_EDGE tmp = kmer_vertex_get_succ_edge(Vertex, i);

		if (tmp->Type == Type) {
			ret = tmp;
			break;
		}
	}

	return ret;
}


PKMER_EDGE _get_refseq_edge(const KMER_VERTEX *Vertex)
{
	return _get_typed_edge(Vertex, kmetReference);
}


PKMER_EDGE _get_refseq_or_variant_edge(const KMER_VERTEX *Vertex)
{
	PKMER_EDGE ret = NULL;

	ret = _get_typed_edge(Vertex, kmetReference);
	if (ret == NULL)
		ret = _get_typed_edge(Vertex, kmetVariant);

	return ret;
}


static ERR_VALUE _capture_edge_sequence(const KMER_EDGE *Start, const POINTER_ARRAY_KMER_EDGE *RSEdges, const KMER_EDGE *End, char **Seq, size_t *SeqLen)
{
	REFSEQ_STORAGE rsStorage;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t kmerSize = kmer_get_size(&Start->Dest->KMer);
	const KMER_VERTEX *lastVertex = NULL;

	ret = ERR_SUCCESS;
	rs_storage_init(&rsStorage);
	if (Start != NULL) {
		ret = rs_storage_add_edge(&rsStorage, Start);
		lastVertex = Start->Dest;
	}

	if (ret == ERR_SUCCESS) {
		if (RSEdges != NULL) {
			for (size_t i = 0; i < pointer_array_size(RSEdges); ++i) {
				ret = rs_storage_add_edge(&rsStorage, RSEdges->Data[i]);
				lastVertex = RSEdges->Data[i]->Dest;
				if (ret != ERR_SUCCESS)
					break;
			}
		}

	}

	if (ret == ERR_SUCCESS) {
		if (End != NULL) {
			ret = rs_storage_add_edge(&rsStorage, End);
			lastVertex = End->Dest;
		}
	}

	if (ret == ERR_SUCCESS) {
		if (lastVertex != NULL && !lastVertex->Helper &&
			lastVertex->Type != kmvtRefSeqEnd)
			rs_storage_remove(&rsStorage, 1);

		ret = rs_storage_create_string(&rsStorage, Seq);
		*SeqLen = rsStorage.ValidLength;
	}

	rs_storage_finit(&rsStorage);

	return ret;
}


ERR_VALUE kmer_graph_create(const uint32_t KMerSize, const size_t VerticesHint, const size_t EdgesHint, PKMER_GRAPH *Graph)
{
	PKMER_GRAPH tmpGraph = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	KMER_TABLE_CALLBACKS vCallbacks;
	KMER_EDGE_TABLE_CALLBACKS eCallbacks;
	KMER_TABLE_CALLBACKS lCallbacks;

	ret = utils_malloc_KMER_GRAPH(&tmpGraph);
	if (ret == ERR_SUCCESS) {
		memset(tmpGraph, 0, sizeof(KMER_GRAPH));
		_init_quality_table(tmpGraph->QualityTable);
		tmpGraph->Allocator.VertexAllocator = _default_vertex_allocator;
		tmpGraph->Allocator.VertexFreer = _default_vertex_freer;
		tmpGraph->Allocator.EdgeAllocator = _default_edge_allocator;
		tmpGraph->Allocator.EdgeFreer = _default_edge_freer;
		tmpGraph->Allocator.VertexAllocatorContext = NULL;
		tmpGraph->Allocator.EdgeAllocatorContext = NULL;
		tmpGraph->NumberOfEdges = 0;
		tmpGraph->NumberOfVertices = 0;
		tmpGraph->KMerSize = KMerSize;
		tmpGraph->StartingVertex = NULL;
		tmpGraph->EndingVertex = NULL;
		tmpGraph->VerticesToDeleteList = NULL;
		vCallbacks.Context = tmpGraph;
		vCallbacks.OnCopy = _vertex_table_on_copy;
		vCallbacks.OnDelete = _vertex_table_on_delete;
		vCallbacks.OnInsert = _vertex_table_on_insert;
		vCallbacks.OnPrint = _vertex_table_on_print;
		ret = kmer_table_create(KMerSize, utils_next_prime(VerticesHint), &vCallbacks, &tmpGraph->VertexTable);
		if (ret == ERR_SUCCESS) {
			eCallbacks.Context = tmpGraph;
			eCallbacks.OnDelete = _edge_table_on_delete;
			eCallbacks.OnInsert = _edge_table_on_insert;
			eCallbacks.OnPrint = _edge_table_on_print;
			ret = kmer_edge_table_create(KMerSize, 4096, &eCallbacks, &tmpGraph->EdgeTable);
			if (ret == ERR_SUCCESS) {
				ret = kmer_edge_table_create(KMerSize, 47, NULL, &tmpGraph->DummyVertices);
				if (ret == ERR_SUCCESS) {
					lCallbacks.Context = Graph;
					lCallbacks.OnCopy = NULL;
					lCallbacks.OnDelete = _kmerlist_table_on_delete;
					lCallbacks.OnInsert = _kmerlist_table_on_insert;
					lCallbacks.OnPrint = NULL;
					ret = kmer_table_create(KMerSize, 37, &lCallbacks, &tmpGraph->KmerListTable);
					if (ret == ERR_SUCCESS) {
						pointer_array_init_KMER_VERTEX(&tmpGraph->RefVertices, 140);
						*Graph = tmpGraph;
					}

					if (ret != ERR_SUCCESS)
						kmer_edge_table_destroy(tmpGraph->DummyVertices);
				}

				if (ret != ERR_SUCCESS)
					kmer_edge_table_destroy(tmpGraph->EdgeTable);
			}

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
	pointer_array_finit_KMER_VERTEX(&Graph->RefVertices);
	kmer_table_destroy(Graph->KmerListTable);
	kmer_edge_table_destroy(Graph->DummyVertices);
	kmer_edge_table_destroy(Graph->EdgeTable);
	kmer_table_destroy(Graph->VertexTable);

	PKMER_VERTEX del = Graph->VerticesToDeleteList;
	PKMER_VERTEX old = Graph->VerticesToDeleteList;

	while (del != NULL) {
		old = del;
		del = del->Lists.Next;
		old->Lists.Graph = NULL;
		_vertex_destroy(Graph, old);
	}

	utils_free(Graph);

	return;
}


void kmer_graph_print(FILE *Stream, const KMER_GRAPH *Graph)
{
	void *it = NULL;
	KMER_VERTEX *v = NULL;

	fprintf(Stream, "digraph G {\n");
	fprintf(Stream, "\t/* number of vertices: %u */\n", Graph->NumberOfVertices);
	fprintf(Stream, "\t/* number of edges: %u */\n", Graph->NumberOfEdges);
	fprintf(Stream, "\t/* number of reference edges: %u */\n", Graph->TypedEdgeCount[kmetReference]);
	fprintf(Stream, "\t/* number of read edges: %u */\n", Graph->TypedEdgeCount[kmetRead]);
	fprintf(Stream, "\t/* number of variant edges: %u */\n", Graph->TypedEdgeCount[kmetVariant]);
	if (kmer_table_first(Graph->VertexTable, &it, &v) == ERR_SUCCESS) {
		do {
			_vertex_table_on_print(Graph->VertexTable, v, Stream);
			for (size_t i = 0; i < pointer_array_size(&v->Successors); ++i)
				_edge_table_on_print(Graph->EdgeTable, v->Successors.Data[i], Graph, Stream);
		} while (kmer_table_next(Graph->VertexTable, it, &it, &v) == ERR_SUCCESS);
	}

	fprintf(Stream, "}\n");
	fflush(Stream);

	return;
}


typedef struct _PATH_COMPARE_CONTEXT {
	const KMER_EDGE **Path;
	size_t PathLength;
	size_t CurrentIndex;
	const KMER_EDGE *CurrentEdge;
	const char *Seq;
	size_t SeqLen;
	boolean End;
} PATH_COMPARE_CONTEXT, *PPATH_COMPARE_CONTEXT;


static void _path_context_init(PPATH_COMPARE_CONTEXT Context, const KMER_EDGE **Path, const size_t PathLength)
{
	Context->Path = Path;
	Context->PathLength = PathLength;
	Context->CurrentIndex = 0;
	Context->CurrentEdge = Context->Path[Context->CurrentIndex];
	Context->Seq = Context->CurrentEdge->Seq;
	Context->SeqLen = Context->CurrentEdge->SeqLen;
	Context->End = FALSE;

	return;
}


static char _path_context_produce_base(PPATH_COMPARE_CONTEXT Context)
{
	char ret = '\0';

	if (Context->SeqLen == 0) {
		if (!Context->CurrentEdge->Dest->Helper && Context->CurrentIndex < Context->PathLength - 1)
			ret = kmer_get_last_base(&Context->CurrentEdge->Dest->KMer);

		++Context->CurrentIndex;
		Context->End = (Context->CurrentIndex == Context->PathLength);
		if (!Context->End) {
			Context->CurrentEdge = Context->Path[Context->CurrentIndex];
			Context->Seq = Context->CurrentEdge->Seq;
			Context->SeqLen = Context->CurrentEdge->SeqLen;
		}
	} else {
		ret = *Context->Seq;
		++Context->Seq;
		--Context->SeqLen;
	}

	return ret;
}



static boolean _paths_equal_by_seq(const KMER_EDGE **Path1, const size_t Path1Length, const KMER_EDGE **Path2, const size_t Path2Length)
{
	boolean ret = TRUE;
	PATH_COMPARE_CONTEXT ctx1;
	PATH_COMPARE_CONTEXT ctx2;

	_path_context_init(&ctx1, Path1, Path1Length);
	_path_context_init(&ctx2, Path2, Path2Length);
	while (ret && !ctx1.End && !ctx2.End) {
		char b1 = '\0';
		while (!ctx1.End && b1 == '\0')
			b1 = _path_context_produce_base(&ctx1);

		char b2 = '\0';
		while (!ctx2.End && b2 == '\0')
			b2 = _path_context_produce_base(&ctx2);
	
		ret = (b1 == b2);
	}

	if (ret && !ctx1.End)
		ret = (_path_context_produce_base(&ctx1) == '\0');

	if (ret && !ctx2.End)
		ret = (_path_context_produce_base(&ctx2) == '\0');

	return ret;
}


void kmer_graph_range(PKMER_GRAPH Graph, uint64_t Start, uint64_t End)
{
	PKMER_EDGE e = NULL;

	e = _get_refseq_or_variant_edge(Graph->StartingVertex);
	while (e->Dest != Graph->EndingVertex) {
		PKMER_VERTEX v = e->Source;

		e = _get_refseq_or_variant_edge(e->Dest);
		if (v->Type == kmvtRefSeqMiddle &&
			(v->AbsPos + 1 < Start || v->AbsPos + 1 > End)) {
			while (pointer_array_size(&v->Successors) > 0)
				kmer_graph_delete_edge(Graph, v->Successors.Data[0]);

			while (pointer_array_size(&v->Predecessors) > 0)
				kmer_graph_delete_edge(Graph, v->Predecessors.Data[0]);

			kmer_graph_delete_vertex(Graph, v);
		} else if (v->AbsPos == Start) {
			v->Type = kmvtRefSeqStart;
			Graph->StartingVertex = v;
		} else if (v->AbsPos == End) {
			v->Type = kmvtRefSeqEnd;
			Graph->EndingVertex = v;
		}
	}

//	size_t dummy = 0;

//	kmer_graph_delete_trailing_things(Graph, &dummy);

	return;
}


void kmer_graph_delete_1to1_vertices(PKMER_GRAPH Graph)
{
	void *iter = NULL;
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	EKMerVertexType vTypes[] = {kmvtRead, kmvtRefSeqMiddle};
	
	for (size_t j = 0; j < sizeof(vTypes) / sizeof(vTypes[0]); ++j) {
		ret = kmer_table_first(Graph->VertexTable, &iter, (void **)&v);
		while (ret == ERR_SUCCESS) {
			if (v->Type == vTypes[j] && kmer_vertex_in_degree(v) == 1 && kmer_vertex_out_degree(v) == 1) {
				const KMER_EDGE *inEdge = kmer_vertex_get_pred_edge(v, 0);
				const KMER_EDGE *outEdge = kmer_vertex_get_succ_edge(v, 0);

				if (!inEdge->LongData.LongEdge && inEdge->Type != kmetVariant && 
					!outEdge->LongData.LongEdge && outEdge->Type != kmetVariant)
					kmer_graph_delete_vertex(Graph, v);
			}

			ret = kmer_table_next(Graph->VertexTable, iter, &iter, (void **)&v);
		}

		if (ret == ERR_NO_MORE_ENTRIES)
			ret = ERR_SUCCESS;
	}

	return;
}


void kmer_graph_delete_edges_under_threshold(PKMER_GRAPH Graph, const size_t Threshold)
{
	PKMER_EDGE e = NULL;
	void *iter = NULL;
	ERR_VALUE err = ERR_INTERNAL_ERROR;
	const size_t realThreshold = Threshold * 100;

	err = kmer_edge_table_first(Graph->EdgeTable, &iter, (void **)&e);
	while (err == ERR_SUCCESS) {
		if (e->Type == kmetRead && read_info_get_count(&e->ReadInfo) <= Threshold)
			kmer_graph_delete_edge(Graph, e);

		err = kmer_edge_table_next(Graph->EdgeTable, iter, &iter, (void **)&e);
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


void kmer_graph_delete_trailing_things(PKMER_GRAPH Graph, size_t *DeletedThings)
{
	void *iter = NULL;
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	POINTER_ARRAY_KMER_VERTEX stack;
	PKMER_EDGE e = NULL;

	pointer_array_init_KMER_VERTEX(&stack, 140);
	ret = kmer_table_first(Graph->VertexTable, &iter, (void **)&v);
	while (ret == ERR_SUCCESS) {
		if (v->Type != kmvtRefSeqStart && v->Type != kmvtRefSeqEnd) {
			if (kmer_vertex_in_degree(v) == 0 ||
				kmer_vertex_out_degree(v) == 0) {
				ret = pointer_array_push_back_KMER_VERTEX(&stack, v);
			}
		}

		ret = kmer_table_next(Graph->VertexTable, iter, &iter, (void **)&v);
	}

	while (pointer_array_size(&stack) > 0) {
		v = *pointer_array_pop_back_KMER_VERTEX(&stack);
		while (kmer_vertex_in_degree(v) > 0) {
			e = kmer_vertex_get_pred_edge(v, 0);
			if (kmer_vertex_out_degree(e->Source) == 1 &&
				kmer_vertex_in_degree(e->Source) > 0 &&
				(e->Source->Type == kmvtRefSeqMiddle || e->Source->Type == kmvtRead))
				ret = pointer_array_push_back_KMER_VERTEX(&stack, e->Source);

			kmer_graph_delete_edge(Graph, e);
		}

		while (kmer_vertex_out_degree(v) > 0) {
			e = kmer_vertex_get_succ_edge(v, 0);
			if (kmer_vertex_in_degree(e->Dest) == 1 &&
				kmer_vertex_out_degree(e->Dest) > 0 &&
				(e->Dest->Type == kmvtRefSeqMiddle || e->Dest->Type == kmvtRead))
				ret = pointer_array_push_back_KMER_VERTEX(&stack, e->Dest);

			kmer_graph_delete_edge(Graph, e);
		}

		kmer_graph_delete_vertex(Graph, v);
		*DeletedThings += 1;
	}

	pointer_array_finit_KMER_VERTEX(&stack);
	/*
	ret = kmer_table_first(Graph->VertexTable, &iter, (void **)&v);
	while (ret != ERR_SUCCESS) {
		boolean deleted = FALSE;

		if (v->Type != kmvtRefSeqStart && v->Type != kmvtRefSeqEnd) {
			if (kmer_vertex_in_degree(v) == 0) {
				deleted = TRUE;
				while (kmer_vertex_out_degree(v) > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_succ_edge(v, 0));
			} else if (kmer_vertex_out_degree(v) == 0) {
				deleted = TRUE;
				while (kmer_vertex_in_degree(v) > 0)
					kmer_graph_delete_edge(Graph, kmer_vertex_get_pred_edge(v, 0));
			}

			if (deleted) {
				kmer_graph_delete_vertex(Graph, v);
				*DeletedThings += 1;
			}
		}

		ret = (deleted) ?
			kmer_table_first(Graph->VertexTable, &iter, (void **)&v) :
			kmer_table_next(Graph->VertexTable, iter, &iter, (void **)&v);
	}
	*/
	return;
}


ERR_VALUE kmer_graph_add_vertex_ex(PKMER_GRAPH Graph, const KMER *KMer, const EKMerVertexType Type, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	v = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
	if (v == NULL) {
		ret = _vertex_create(Graph, KMer, Type, &v);
		if (ret == ERR_SUCCESS) {
			ret = kmer_table_insert(Graph->VertexTable, &v->KMer, v);
			if (ret == ERR_SUCCESS) {
				PKMER_LIST list;
				PKMER lk = NULL;

				KMER_STACK_ALLOC(lk, 0, kmer_get_size(KMer), KMer->Bases);
				list = (PKMER_LIST)kmer_table_get(Graph->KmerListTable, lk);
				if (list == NULL) {
					ret = utils_malloc(sizeof(KMER_LIST) + kmer_get_size(lk)*sizeof(char), &list);
					if (ret == ERR_SUCCESS) {
						pointer_array_init_KMER_VERTEX(&list->Vertices, 140);
						kmer_init_from_kmer(&list->Kmer, lk);
						ret = kmer_table_insert(Graph->KmerListTable, &list->Kmer, list);
						if (ret != ERR_SUCCESS)
							utils_free(list);
					}
				}
				
				if (ret == ERR_SUCCESS)
					ret = pointer_array_push_back_KMER_VERTEX(&list->Vertices, v);
				
				assert(kmer_seq_equal(&list->Kmer, &v->KMer));
				if (ret == ERR_SUCCESS) {
					Graph->NumberOfVertices++;
					v->Lists.Graph = Graph;
					*Vertex = v;
				}

				if (ret != ERR_SUCCESS)
					kmer_table_delete(Graph->VertexTable, KMer);
			}

			if (ret != ERR_SUCCESS)
				_vertex_destroy(Graph, v);
		}
	} else {
		*Vertex = v;
		ret = ERR_ALREADY_EXISTS;
	}

	return ret;
}


ERR_VALUE kmer_graph_add_helper_vertex(PKMER_GRAPH Graph, const KMER *KMer1, const KMER *KMer2, PKMER_VERTEX *Vertex)
{
	PKMER_VERTEX v = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *dummySeq = alloca(kmerSize*sizeof(char));
	PKMER kmer = NULL;

	ret = ERR_SUCCESS;
	memset(dummySeq, 'D', kmerSize*sizeof(char));
	KMER_STACK_ALLOC(kmer, 0, kmerSize, dummySeq);
	v = (PKMER_VERTEX)kmer_edge_table_get(Graph->DummyVertices, KMer1, KMer2);
	if (v == NULL) {
		do {
			ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
			if (ret == ERR_ALREADY_EXISTS)
				kmer_get_number(kmer) += 1;
		} while (ret == ERR_ALREADY_EXISTS);
	
		if (ret == ERR_SUCCESS) {
			v->Helper = TRUE;
			ret = kmer_edge_table_insert(Graph->DummyVertices, KMer1, KMer2, v);
			if (ret == ERR_SUCCESS)
				*Vertex = v;

			if (ret != ERR_SUCCESS)
				kmer_graph_delete_vertex(Graph, v);
		}
	} else {
		ret = ERR_ALREADY_EXISTS;
		*Vertex = v;
	}


	return ret;
}


ERR_VALUE kmer_graph_add_edge_ex(PKMER_GRAPH Graph, PKMER_VERTEX Source, PKMER_VERTEX Dest, const EKMerEdgeType Type, PKMER_EDGE *Edge)
{
	PKMER_EDGE edge = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	edge = (PKMER_EDGE)kmer_edge_table_get(Graph->EdgeTable, &Source->KMer, &Dest->KMer);
	if (edge == NULL) {
		ret = _edge_create(Graph, Source, Dest, Type, &edge);
		if (ret == ERR_SUCCESS) {
			ret = kmer_edge_table_insert(Graph->EdgeTable, &Source->KMer, &Dest->KMer, edge);
			if (ret == ERR_SUCCESS) {
				ret = pointer_array_reserve_KMER_EDGE(&Source->Successors, pointer_array_size(&Source->Successors) + 1);
				if (ret == ERR_SUCCESS) {
					ret = pointer_array_reserve_KMER_EDGE(&Dest->Predecessors, pointer_array_size(&Dest->Predecessors) + 1);
					if (ret == ERR_SUCCESS) {
						pointer_array_push_back_no_alloc_KMER_EDGE(&Source->Successors, edge);
						pointer_array_push_back_no_alloc_KMER_EDGE(&Dest->Predecessors, edge);
						Graph->NumberOfEdges++;
						Graph->TypedEdgeCount[Type]++;
						*Edge = edge;
					}
				}
			}

			if (ret != ERR_SUCCESS)
				_edge_destroy(Graph, edge);
		}
	} else {
		*Edge = edge;
		ret = ERR_ALREADY_EXISTS;
	}

	return ret;
}


PKMER_EDGE kmer_graph_get_edge(const struct _KMER_GRAPH *Graph, const struct _KMER *Source, const struct _KMER *Dest)
{
	return (PKMER_EDGE)kmer_edge_table_get(Graph->EdgeTable, Source, Dest);
}


PKMER_VERTEX kmer_graph_get_vertex(const struct _KMER_GRAPH *Graph, const struct _KMER *KMer)
{
	return (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, KMer);
}


ERR_VALUE kmer_graph_get_vertices(const KMER_GRAPH *Graph, const KMER *KMer, PPOINTER_ARRAY_KMER_VERTEX *VertexArray)
{
	PKMER lk = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const KMER_LIST *l = NULL;

	ret = ERR_NOT_FOUND;
	KMER_STACK_ALLOC(lk, 0, kmer_get_size(KMer), KMer->Bases);
	l = (PKMER_LIST)kmer_table_get(Graph->KmerListTable, lk);
	if (l != NULL) {
		for (size_t i = 0; i < pointer_array_size(&l->Vertices); ++i) {
			assert(kmer_seq_equal(KMer, &l->Vertices.Data[i]->KMer));
			assert(kmer_seq_equal(&l->Kmer, &l->Vertices.Data[i]->KMer));
		}

		*VertexArray = &l->Vertices;
		ret = ERR_SUCCESS;
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

	if (kmer_vertex_in_degree(Vertex) <= 1 && kmer_vertex_out_degree(Vertex) <= 1) {
		if (kmer_vertex_out_degree(Vertex) == 1) {
			outEdge = kmer_vertex_get_succ_edge(Vertex, 0);
			succVertex = outEdge->Dest;
		}

		if (kmer_vertex_in_degree(Vertex) == 1) {
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
				PKMER_LIST list = NULL;
				PKMER lk = NULL;

				KMER_STACK_ALLOC(lk, 0, kmer_get_size(&Vertex->KMer), Vertex->KMer.Bases);
				list = (PKMER_LIST)kmer_table_get(Graph->KmerListTable, lk);
				pointer_array_remove_KMER_VERTEX(&list->Vertices, Vertex);
				ret = kmer_table_delete(Graph->VertexTable, &Vertex->KMer);
				--Graph->NumberOfVertices;
			}
		} else ret = ERR_PRED_IS_SUCC;
	} else ret = ERR_TOO_MANY_EDGES;

	return ret;
}


void kmer_graph_delete_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge)
{
	ERR_VALUE err = ERR_INTERNAL_ERROR;
	PKMER_VERTEX source = Edge->Source;
	PKMER_VERTEX dest = Edge->Dest;
	EKMerEdgeType edgeType = Edge->Type;

	if (Graph->DeleteEdgeCallback != NULL)
		Graph->DeleteEdgeCallback(Graph, Edge, Graph->DeleteEdgeCallbackContext);

	err = kmer_edge_table_delete(Graph->EdgeTable, &source->KMer, &dest->KMer);
	if (err == ERR_SUCCESS) {
		pointer_array_remove_by_item_fast_KMER_EDGE(&source->Successors, Edge);
		pointer_array_remove_by_item_fast_KMER_EDGE(&dest->Predecessors, Edge);
		--Graph->TypedEdgeCount[edgeType];
		--Graph->NumberOfEdges;
	}

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
		if (kmer_graph_get_edge(Graph, &u->KMer, &w->KMer) == NULL) {
			boolean mfd = (Source->MarkedForDelete || Dest->MarkedForDelete);
			EKMerEdgeType type = (Source->Type == kmetReference && Dest->Type == kmetReference) ? kmetReference : kmetRead;
			PKMER_EDGE newEdge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, u, w, type, &newEdge);
			if (ret == ERR_SUCCESS) {
				char *newSeq = NULL;

				newEdge->MarkedForDelete = mfd;
				ret = _capture_edge_sequence(Source, NULL, Dest, &newSeq, &newEdge->SeqLen);
				if (ret == ERR_SUCCESS) {
					newEdge->Seq = newSeq;
					ret = dym_array_reserve_size_t(&newEdge->Weights, newEdge->SeqLen + 1);
					if (ret == ERR_SUCCESS) {						
						dym_array_push_back_array_no_alloc_size_t(&newEdge->Weights, &Source->Weights);
						dym_array_push_back_array_no_alloc_size_t(&newEdge->Weights, &Dest->Weights);
						newEdge->SeqType = newEdge->SeqType;
						newEdge->Seq1Weight = max(Source->Seq1Weight, Dest->Seq1Weight);
						ret = pointer_array_reserve_READ_INFO(&newEdge->ReadIndices, gen_array_size(&newEdge->Weights));
						if (ret == ERR_SUCCESS) {
							pointer_array_push_back_array_no_alloc_READ_INFO(&newEdge->ReadIndices, &Source->ReadIndices);
							pointer_array_push_back_array_no_alloc_READ_INFO(&newEdge->ReadIndices, &Dest->ReadIndices);
							pointer_array_clear_READ_INFO(&Source->ReadIndices);
							pointer_array_clear_READ_INFO(&Dest->ReadIndices);
						}

						if (ret == ERR_SUCCESS)
							ret = read_info_merge(&newEdge->ReadInfo, &Source->ReadInfo, &Dest->ReadInfo);
						
						if (ret == ERR_SUCCESS) {
							kmer_graph_delete_edge(Graph, Source);
							kmer_graph_delete_edge(Graph, Dest);
						}
					}
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, newEdge);
			}
		} else ret = ERR_TRIANGLE;
	} else ret = (u == w) ? ERR_PRED_IS_SUCC : ERR_NOT_ADJACENT;

	return ret;
}


ERR_VALUE kmer_graph_get_splitted_edge(PKMER_GRAPH Graph, const KMER_VERTEX *Source, const KMER_VERTEX *Dest, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX helperVertex = NULL;
	PKMER_EDGE es = NULL;
	PKMER_EDGE ed = NULL;

	helperVertex = (PKMER_VERTEX)kmer_edge_table_get(Graph->DummyVertices, &Source->KMer, &Dest->KMer);
	if (helperVertex != NULL) {
		es = kmer_graph_get_edge(Graph, &Source->KMer, &helperVertex->KMer);
		ed = kmer_graph_get_edge(Graph, &helperVertex->KMer, &Dest->KMer);
		ret = ERR_SUCCESS;
		if (SplitVertex != NULL)
			*SplitVertex = helperVertex;

		if (SourceEdge != NULL)
			*SourceEdge = es;

		if (DestEdge != NULL)
			*DestEdge = ed;
	} else ret = ERR_NOT_FOUND;

	return ret;
}


ERR_VALUE kmer_graph_split_edge(PKMER_GRAPH Graph, PKMER_EDGE Edge, PKMER_EDGE *SourceEdge, PKMER_EDGE *DestEdge, PKMER_VERTEX *SplitVertex)
{
	PKMER_EDGE es = NULL;
	PKMER_EDGE ed = NULL;
	PKMER_VERTEX helperVertex = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_add_helper_vertex(Graph, &Edge->Source->KMer, &Edge->Dest->KMer, &helperVertex);
	if (ret == ERR_SUCCESS) {
		ret = kmer_graph_add_edge_ex(Graph, Edge->Source, helperVertex, Edge->Type, &es);
		if (ret == ERR_SUCCESS) {
			helperVertex->RefSeqPosition = Edge->Dest->RefSeqPosition;
			if (Edge->Source->Type == kmvtRefSeqMiddle && Edge->Source->RefSeqPosition > Edge->Dest->RefSeqPosition)
				helperVertex->RefSeqPosition = Edge->Source->RefSeqPosition + 1;

			ret = kmer_graph_add_edge_ex(Graph, helperVertex, Edge->Dest, Edge->Type, &ed);
			if (ret == ERR_SUCCESS) {
				char *seq = NULL;
				size_t seqLen = Edge->SeqLen;

				if (seqLen > 0) {
					ret = utils_calloc(seqLen + 1, sizeof(char), &seq);
					if (ret == ERR_SUCCESS) {
						memcpy(seq, Edge->Seq, seqLen*sizeof(char));
						seq[seqLen] = '\0';
						kmer_edge_add_seq(es, Edge->SeqType, seq, seqLen);
						seq = NULL;
					}
				}

				es->Seq1Weight = Edge->Seq1Weight;
				ed->Seq1Weight = Edge->Seq1Weight;
				if (ret == ERR_SUCCESS) {
					ret = read_info_copy(&es->ReadInfo, &Edge->ReadInfo);
					if (ret == ERR_SUCCESS) {
						ret = read_info_copy(&ed->ReadInfo, &Edge->ReadInfo);
						if (ret == ERR_SUCCESS) {
							kmer_graph_delete_edge(Graph, Edge);
							if (SourceEdge != NULL)
								*SourceEdge = es;

							if (DestEdge != NULL)
								*DestEdge = ed;

							if (SplitVertex != NULL)
								*SplitVertex = helperVertex;
						}
					}
				}

				if (ret != ERR_SUCCESS)
					kmer_graph_delete_edge(Graph, ed);
			}

			if (ret != ERR_SUCCESS)
				kmer_graph_delete_edge(Graph, es);
		}

		if (ret != ERR_SUCCESS)
			kmer_graph_delete_vertex(Graph, helperVertex);
	} else if (ret == ERR_ALREADY_EXISTS) {
		if (SplitVertex != NULL)
			*SplitVertex = helperVertex;

		if (SourceEdge != NULL)
			*SourceEdge = kmer_graph_get_edge(Graph, &Edge->Source->KMer, &helperVertex->KMer);

		if (DestEdge != NULL)
			*DestEdge = kmer_graph_get_edge(Graph, &helperVertex->KMer, &Edge->Dest->KMer);
	}

	assert(*SourceEdge != NULL);
	assert(*DestEdge != NULL);
	assert(*SplitVertex != NULL);

	return ret;
}


ERR_VALUE kmer_vertex_get_certain_edges(const KMER_VERTEX *Vertex, const EKMerEdgeType EdgeType, const boolean Incomming, PPOINTER_ARRAY_KMER_EDGE Array)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const PKMER_EDGE *edgeArray = (Incomming) ? (Vertex->Predecessors.Data) : (Vertex->Successors.Data);
	const size_t count = pointer_array_size((Incomming) ? (&Vertex->Predecessors) : (&Vertex->Successors));

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < count; ++i) {
		const KMER_EDGE *e = *edgeArray;
	
		if (e->Type == EdgeType)
			ret = pointer_array_push_back_KMER_EDGE(Array, e);
		
		if (ret != ERR_SUCCESS)
			break;

		++edgeArray;
	}

	return ret;
}


void kmer_edge_add_seq(PKMER_EDGE Edge, EKMerEdgeType Type, const char *Seq, const size_t Length)
{
	if (Edge->SeqType == kmetNone) {
		Edge->SeqType = Type;
		Edge->Seq = Seq;
		Edge->SeqLen = Length;
	} else {
		printf("Attempt to add second sequence to a non-variant edge\n");
		exit(0);
	}

	return;
}


static void _remove_read_info_from_edge_array(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY ToRemove, size_t Distance, size_t *NewDistance)
{
	for (size_t k = 0; k < pointer_array_size(Array); ++k) {
		PKMER_EDGE tmp = *pointer_array_item_KMER_EDGE(Array, k);

		read_info_subtract(&tmp->ReadInfo, ToRemove, Distance);
		Distance += tmp->SeqLen;
		if (!tmp->Dest->Helper)
			Distance++;
	}

	if (NewDistance != NULL)
		*NewDistance = Distance;

	return;
}


static void _remove_read_info_from_edges(PKMER_GRAPH Graph, PKMER_EDGE Edge1, PKMER_EDGE Edge2, PPOINTER_ARRAY_KMER_EDGE Array, PGEN_ARRAY_READ_INFO_ENTRY Intersection)
{
	size_t distance = 0;

	read_info_subtract(&Edge1->ReadInfo, Intersection, distance);
	distance += Edge1->SeqLen;
	if (!Edge1->Dest->Helper)
		++distance;

	_remove_read_info_from_edge_array(Graph, Array, Intersection, distance, &distance);
	read_info_subtract(&Edge2->ReadInfo, Intersection, distance);

	return;
}


typedef struct _EDGE_REMOVE_CONTEXT {
	POINTER_ARRAY_KMER_EDGE RSEdges;
	PKMER_EDGE TargetEdge;
	size_t ReadDistance;
} EDGE_REMOVE_CONTEXT, *PEDGE_REMOVE_CONTEXT;

POINTER_ARRAY_TYPEDEF(EDGE_REMOVE_CONTEXT);
POINTER_ARRAY_IMPLEMENTATION(EDGE_REMOVE_CONTEXT)


static ERR_VALUE _remove_context_create(PPOINTER_ARRAY_KMER_EDGE RSEdges, PKMER_EDGE TargetEdge, const size_t ReadDistance, PEDGE_REMOVE_CONTEXT *Context)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PEDGE_REMOVE_CONTEXT tmpContext = NULL;

	ret = utils_malloc(sizeof(EDGE_REMOVE_CONTEXT), &tmpContext);
	if (ret == ERR_SUCCESS) {
		tmpContext->ReadDistance = ReadDistance;
		tmpContext->TargetEdge = TargetEdge;
		pointer_array_init_KMER_EDGE(&tmpContext->RSEdges, 140);
		ret = pointer_array_push_back_array_KMER_EDGE(&tmpContext->RSEdges, RSEdges);
		if (ret == ERR_SUCCESS)
			*Context = tmpContext;
		
		if (ret != ERR_SUCCESS) {
			pointer_array_finit_KMER_EDGE(&tmpContext->RSEdges);
			utils_free(tmpContext);
		}
	}

	return ret;
}


static void _remove_context_destroy(PEDGE_REMOVE_CONTEXT Context)
{
	pointer_array_finit_KMER_EDGE(&Context->RSEdges);
	utils_free(Context);

	return;
}

static void _remove_context_apply(PKMER_GRAPH Graph, PEDGE_REMOVE_CONTEXT Context)
{
	_remove_read_info_from_edge_array(Graph, &Context->RSEdges, &Context->TargetEdge->ReadInfo.Array, Context->ReadDistance, NULL);
	_remove_context_destroy(Context);

	return;
}

ERR_VALUE kmer_graph_connect_reads_by_pairs(PKMER_GRAPH Graph, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, size_t *ChangeCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_READ_INFO_ENTRY intersection;
	PKMER_EDGE rsLastEdge = NULL;
	PKMER_EDGE rsNextEdge = NULL;
	const KMER_EDGE *nextRsLastEdge = NULL;
	POINTER_ARRAY_KMER_EDGE edgesToDelete;
	POINTER_ARRAY_KMER_EDGE rsEdges;
	size_t dummy = 0;
	boolean deleteEOut = FALSE;
	boolean deleteeIn = FALSE;
	POINTER_ARRAY_TYPE(EDGE_REMOVE_CONTEXT) removeContexts;

	pointer_array_init_EDGE_REMOVE_CONTEXT(&removeContexts, 140);
	ret = pointer_array_reserve_EDGE_REMOVE_CONTEXT(&removeContexts, 2 * gen_array_size(PairArray));
	if (ret == ERR_SUCCESS) {
		dym_array_init_READ_INFO_ENTRY(&intersection, 140);
		pointer_array_init_KMER_EDGE(&edgesToDelete, 140);
		ret = pointer_array_reserve_KMER_EDGE(&edgesToDelete, 2 * gen_array_size(PairArray));
		if (ret == ERR_SUCCESS) {
			for (size_t h = 0; h < gen_array_size(PairArray); ++h) {
				boolean edgeCreated = FALSE;
				KMER_EDGE_PAIR pair = *dym_array_item_KMER_EDGE_PAIR(PairArray, h);
				PKMER_EDGE eIn = pair.U;
				PKMER_EDGE eOut = pair.V;
				boolean allowed = TRUE;

				if (eIn == eOut)
					continue;

				if (pair.U == NULL || pair.V == NULL || pair.ConnectingEdge == NULL)
					continue;

				if (eIn->Dest->Type == kmvtRefSeqMiddle && !eIn->Source->LongEdgeAllowed)
					allowed = FALSE;

				rsLastEdge = _get_in_refseq_edge(eIn->Dest);
				rsNextEdge = _get_refseq_edge(eOut->Source);
				rsEdges.Data = pair.Edges;
				rsEdges.ValidLength = pair.EdgeCount;
				rsEdges.AllocLength = pair.EdgeCount;
				if (pointer_array_size(&rsEdges) != pair.ReadDistance)
					continue;

					ret = read_info_intersection(&eIn->ReadInfo, &eOut->ReadInfo, &intersection, eIn->SeqLen + (!eIn->Dest->Helper ? 1 : 0) + pair.ReadDistance);
					if (ret == ERR_SUCCESS && gen_array_size(&intersection) > Threshold && allowed) {
						if (kmer_equal(&eIn->Source->KMer, &eOut->Dest->KMer))
							ret = ERR_ALREADY_EXISTS;

						if (ret == ERR_SUCCESS) {
							*ChangeCount++;
							edgeCreated = TRUE;
							_remove_read_info_from_edges(Graph, eIn, eOut, &rsEdges, &intersection);
						}
						
						ret = ERR_SUCCESS;
					} else if (ret == ERR_SUCCESS) {
						PKMER_EDGE e = kmer_graph_get_edge(Graph, &eIn->Source->KMer, &eOut->Dest->KMer);

						if (e != NULL)
							kmer_graph_delete_edge(Graph, e);
					}

					if (edgeCreated) {
						deleteeIn = FALSE;
						deleteEOut = FALSE;
						if (rsNextEdge != NULL) {
							ret = read_info_intersection(&eIn->ReadInfo, &rsNextEdge->ReadInfo, &intersection, eIn->SeqLen + 1 + pair.ReadDistance);
							if (ret == ERR_SUCCESS)
								deleteeIn = (gen_array_size(&intersection) <= Threshold);
						} else deleteeIn = TRUE;

						if (ret == ERR_SUCCESS) {
							if (rsLastEdge != NULL) {
								ret = read_info_intersection(&rsLastEdge->ReadInfo, &eOut->ReadInfo, &intersection, rsLastEdge->SeqLen + 1 + pair.ReadDistance);
								if (ret == ERR_SUCCESS)
									deleteEOut = (gen_array_size(&intersection) <= Threshold);
							} else deleteEOut = TRUE;
						}

						if (deleteeIn) {
							PEDGE_REMOVE_CONTEXT removeContext = NULL;

							ret = _remove_context_create(&rsEdges, eIn, eIn->SeqLen + 1, &removeContext);
							if (ret == ERR_SUCCESS) {
								pointer_array_push_back_no_alloc_EDGE_REMOVE_CONTEXT(&removeContexts, removeContext);
								if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eIn))
									pointer_array_push_back_no_alloc_KMER_EDGE(&edgesToDelete, eIn);
							}
						}

						if (deleteEOut) {
							PEDGE_REMOVE_CONTEXT removeContext = NULL;

							ret = _remove_context_create(&rsEdges, eOut, -pair.ReadDistance, &removeContext);
							if (ret == ERR_SUCCESS) {
								pointer_array_push_back_no_alloc_EDGE_REMOVE_CONTEXT(&removeContexts, removeContext);
								if (!pointer_array_contains_KMER_EDGE(&edgesToDelete, eOut))
									pointer_array_push_back_no_alloc_KMER_EDGE(&edgesToDelete, eOut);
							}
						}
					}

				if (ret != ERR_SUCCESS)
					break;
			}

			dym_array_finit_READ_INFO_ENTRY(&intersection);

			if (ret == ERR_NO_MORE_ENTRIES)
				ret = ERR_SUCCESS;

			PEDGE_REMOVE_CONTEXT *pctx = removeContexts.Data;
			for (size_t i = 0; i < pointer_array_size(&removeContexts); ++i) {
				_remove_context_apply(Graph, *pctx);
				++pctx;
			}

			for (size_t i = 0; i < pointer_array_size(&edgesToDelete); ++i) {
				kmer_graph_delete_edge(Graph, *pointer_array_item_KMER_EDGE(&edgesToDelete, i));
				*ChangeCount++;
			}
		}

		pointer_array_finit_KMER_EDGE(&edgesToDelete);
	}

	pointer_array_finit_EDGE_REMOVE_CONTEXT(&removeContexts);

	kmer_graph_delete_trailing_things(Graph, &dummy);

	return ret;
}


KHASH_SET_INIT_INT(es);


static ERR_VALUE _follow_line(PKMER_EDGE Start, char **Seq, size_t *SeqLen, PKMER_EDGE *LastEdge, PGEN_ARRAY_size_t Weights, PPOINTER_ARRAY_READ_INFO ReadPath, PPOINTER_ARRAY_KMER_EDGE Edges)
{
	REFSEQ_STORAGE rs;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	rs_storage_init(&rs);
	while (ret == ERR_SUCCESS && kmer_vertex_in_degree(Start->Dest) == 1 &&
		kmer_vertex_out_degree(Start->Dest) == 1) {
		ret = rs_storage_add_edge(&rs, Start);
		if (ret == ERR_SUCCESS)
			ret = dym_array_push_back_array_size_t(Weights, &Start->Weights);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_array_READ_INFO(ReadPath, &Start->ReadIndices);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_KMER_EDGE(Edges, Start);

		Start = kmer_vertex_get_succ_edge(Start->Dest, 0);
	}

	if (ret == ERR_SUCCESS) {
		ret = rs_storage_add_edge(&rs, Start);
		if (ret == ERR_SUCCESS)
			ret = dym_array_push_back_array_size_t(Weights, &Start->Weights);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_array_READ_INFO(ReadPath, &Start->ReadIndices);

		if (ret == ERR_SUCCESS)
			ret = pointer_array_push_back_KMER_EDGE(Edges, Start);

		if (ret == ERR_SUCCESS)
			ret = rs_storage_create_string(&rs, Seq);
	
		if (ret == ERR_SUCCESS) {
			*SeqLen = strlen(*Seq);
			*LastEdge = Start;
		}
	}

	rs_storage_finit(&rs);

	return ret;
}


static uint32_t _binomic_probability(const size_t RefReads, const size_t AltReads)
{
	uint32_t ret = 0;
	double res = 1.0;
	const size_t n = (RefReads + AltReads) / 100;
	const double kn = (double)AltReads / (n*100);
	const size_t k = AltReads / 100;

	if (AltReads > 0 && RefReads > 0) {
		for (size_t i = 0; i < k; ++i) {
			res *= (n - i);
			res /= (i + 1);
			res *= kn;
		}

		for (size_t i = 0; i < n - k; ++i)
			res *= (1 - kn);
	} else if (RefReads == 0)
		res = 1.0;
	else res = 0;

	ret = (uint32_t)round(res * 100);

	return ret;
}


static ERR_VALUE _create_variants(const uint32_t KMerSize, const char *Chrom, uint64_t Pos, const char *Ref, size_t RefLen, const char *Alt, size_t AltLen, const GEN_ARRAY_size_t *RSWeights, const GEN_ARRAY_size_t *ReadWeights, const POINTER_ARRAY_READ_INFO *RefReads, const POINTER_ARRAY_READ_INFO *AltReads, void *Context, const PARSE_OPTIONS *Options, PGEN_ARRAY_VARIANT_CALL VCArray)
{
	VARIANT_CALL vc;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int rfwStartIndex = 0;
	int rfwEndIndex = 0;
	int rewStartIndex = 0;
	int rewEndIndex = 0;
	char *opString = NULL;
	size_t opStringLen = 0;

	++Alt;
	++Ref;
	--AltLen;
	--RefLen;
	++Pos;

	assert(gen_array_size(RSWeights) >= RefLen);
	assert(gen_array_size(ReadWeights) >= AltLen);
	assert(pointer_array_size(RefReads) == gen_array_size(RSWeights));
	assert(pointer_array_size(AltReads) == gen_array_size(ReadWeights));

	while (RefLen != 0 || AltLen != 0) {
		ret = ssw_clever(Ref, RefLen, Alt, AltLen, 2, -1, -1, &opString, &opStringLen);;
		if (ret == ERR_SUCCESS) {
			const char *opIt = opString;
			boolean nothing = TRUE;
			const char *tmpRS = Ref;
			const char *tmpAltS = Alt;
			boolean found = FALSE;

			while (!found) {
				switch (*opIt) {
				case 'X':
					++tmpRS;
					++rfwEndIndex;
					++rewEndIndex;
					++tmpAltS;
					nothing = FALSE;
					break;
				case '\0':
				case 'M':
					if (!nothing) {
						size_t rLen = tmpRS - Ref;
						size_t aLen = tmpAltS - Alt;
						uint32_t offset = 0;

						if (rLen == 0 || aLen == 0 ||
							((rLen > 1 || aLen > 1) && *Ref != *Alt))
							offset++;

						GEN_ARRAY_size_t refIndices;
						GEN_ARRAY_size_t altIndices;

						dym_array_init_size_t(&refIndices, 140);
						for (int i = rfwStartIndex; i < rfwEndIndex + 1; ++i)
							read_info_to_indices(RefReads->Data[i], &refIndices);

						dym_array_init_size_t(&altIndices, 140);
						for (int i = rewStartIndex; i < rewEndIndex + 1; ++i)
							read_info_to_indices(AltReads->Data[i], &altIndices);

						ret = variant_call_init(Chrom, Pos + 1 - offset, ".", Ref - offset, rLen + offset, Alt - offset, aLen + offset, 60, &refIndices, &altIndices, &vc);
						if (ret == ERR_SUCCESS) {
							size_t total = 0;

							vc.KMerSize = KMerSize;
							vc.Context = Context;
							for (int i = rfwStartIndex; i < rfwEndIndex + 1; ++i)
								total += RSWeights->Data[i];

							vc.RefWeight = (total / (rfwEndIndex + 1 - rfwStartIndex));
							total = 0;
							vc.AltWeight = 0;
							for (int i = rewStartIndex; i < rewEndIndex + 1; ++i)
								total += ReadWeights->Data[i];

							vc.AltWeight = (total / (rewEndIndex + 1 - rewStartIndex));
							vc.ProbByCounts = _binomic_probability(100*gen_array_size(&vc.RefReads), 100*gen_array_size(&vc.AltReads));
							vc.ProbByWeights = _binomic_probability(vc.RefWeight, vc.AltWeight);
							
							if (Options->ReadCoverage > 0) {
								vc.AltWeight /= Options->ReadCoverage;
								vc.RefWeight /= Options->ReadCoverage;
							}

							ret = vc_array_add(VCArray, &vc, NULL);
							if (ret != ERR_SUCCESS) {
								variant_call_finit(&vc);
								if (ret == ERR_ALREADY_EXISTS)
									ret = ERR_SUCCESS;
							}
						}

						dym_array_finit_size_t(&altIndices);
						dym_array_finit_size_t(&refIndices);

						Pos += (tmpRS - Ref);
						RefLen -= rLen;
						AltLen -= aLen;
						rfwStartIndex = rfwEndIndex;
						rewStartIndex = rewEndIndex;
						Ref = tmpRS;
						Alt = tmpAltS;
						nothing = TRUE;
					}

					while ((RefLen > 0 && AltLen > 0) && *Alt == *Ref) {
						Pos++;
						++rfwStartIndex;
						++rewStartIndex;
						Ref++;
						Alt++;
						--RefLen;
						--AltLen;
						++rfwEndIndex;
						++rewEndIndex;
					}

					found = TRUE;
					break;
				case 'I':
					++tmpAltS;
					++rewEndIndex;
					nothing = FALSE;
					break;
				case 'D':
					++tmpRS;
					++rfwEndIndex;
					nothing = FALSE;
					break;
				}

				if (*opIt == '\0')
					break;

				++opIt;
			}

			utils_free(opString);
		}
	}

	return ret;
}


ERR_VALUE kmer_graph_detect_uncertainities(PKMER_GRAPH Graph, PGEN_ARRAY_VARIANT_CALL VCArray, const char *CHrom, const PARSE_OPTIONS *Options, boolean *Changed)
{
	boolean edgeCreated = FALSE;
	REFSEQ_STORAGE s1;
	REFSEQ_STORAGE s2;
	GEN_ARRAY_size_t w1;
	GEN_ARRAY_size_t w2;
	POINTER_ARRAY_READ_INFO rp1;
	POINTER_ARRAY_READ_INFO rp2;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE e = NULL;
	PKMER_VERTEX v = NULL;
	POINTER_ARRAY_KMER_EDGE es1;
	POINTER_ARRAY_KMER_EDGE es2;

	ret = ERR_SUCCESS;
	pointer_array_init_KMER_EDGE(&es1, 140);
	pointer_array_init_KMER_EDGE(&es2, 140);
	rs_storage_init(&s1);
	rs_storage_init(&s2);
	dym_array_init_size_t(&w1, 140);
	dym_array_init_size_t(&w2, 140);
	pointer_array_init_READ_INFO(&rp1, 140);
	pointer_array_init_READ_INFO(&rp2, 140);
	v = _get_refseq_edge(Graph->StartingVertex)->Dest;
	while (v != Graph->EndingVertex) {
		edgeCreated = FALSE;
		if (v->Type == kmvtRefSeqMiddle && kmer_vertex_out_degree(v) >= 2 &&
			kmer_vertex_in_degree(v) == 1) {
			for (size_t m = 0; m < kmer_vertex_out_degree(v); ++m) {
				PKMER_EDGE path1Start = _get_refseq_or_variant_edge(v);
				PKMER_EDGE path2Start = kmer_vertex_get_succ_edge(v, m);

				if (path1Start == path2Start)
					continue;

				size_t weight1 = path1Start->Seq1Weight;
				size_t weight2 = path2Start->Seq1Weight;


				PKMER_VERTEX path1Vertex = path1Start->Dest;
				PKMER_VERTEX path2Vertex = path2Start->Dest;

				dym_array_clear_size_t(&w1);
				pointer_array_clear_KMER_EDGE(&es1);
				pointer_array_clear_READ_INFO(&rp1);
				rs_storage_reset(&s1);
				rs_storage_add_vertex(&s1, path1Start->Source);
				rs_storage_add_edge(&s1, path1Start);
				dym_array_push_back_array_size_t(&w1, &path1Start->Weights);
				pointer_array_push_back_array_READ_INFO(&rp1, &path1Start->ReadIndices);
				pointer_array_push_back_KMER_EDGE(&es1, path1Start);
				while (ret == ERR_SUCCESS && kmer_vertex_in_degree(path1Vertex) == 1 && kmer_vertex_out_degree(path1Vertex) == 1 && path1Vertex->Type == kmvtRefSeqMiddle ) {
					PKMER_EDGE e = NULL;

					e = kmer_vertex_get_succ_edge(path1Vertex, 0);
					weight1 = max(weight1, e->Seq1Weight);
					path1Vertex = e->Dest;
					rs_storage_add_edge(&s1, e);
					dym_array_push_back_array_size_t(&w1, &e->Weights);
					pointer_array_push_back_array_READ_INFO(&rp1, &e->ReadIndices);
					pointer_array_push_back_KMER_EDGE(&es1, e);
				}

				if (!path1Vertex->Helper)
					rs_storage_remove(&s1, 1);

				{
					khash_t(es)	*table = kh_init(es);

					pointer_array_clear_KMER_EDGE(&es2);
					pointer_array_clear_READ_INFO(&rp2);
					dym_array_clear_size_t(&w2);
					rs_storage_reset(&s2);
					rs_storage_add_vertex(&s2, path2Start->Source);
					rs_storage_add_edge(&s2, path2Start);
					dym_array_push_back_array_size_t(&w2, &path2Start->Weights);
					pointer_array_push_back_array_READ_INFO(&rp2, &path2Start->ReadIndices);
					pointer_array_push_back_KMER_EDGE(&es2, path2Start);
					while (ret == ERR_SUCCESS && kmer_vertex_out_degree(path2Vertex) == 1 && path2Vertex->Type == kmvtRead) {
						int r;
						PKMER_EDGE e = NULL;

						e = kmer_vertex_get_succ_edge(path2Vertex, 0);
						if (kh_get(es, table, e->Order) != kh_end(table))
							break;

						dym_array_push_back_array_size_t(&w2, &e->Weights);
						pointer_array_push_back_array_READ_INFO(&rp2, &e->ReadIndices);
						pointer_array_push_back_KMER_EDGE(&es2, e);
						ret = rs_storage_add_edge(&s2, e);
						weight2 = max(weight2, e->Seq1Weight);
						path2Vertex = e->Dest;
						kh_put(es, table, e->Order, &r);
					}
					
					if (ret == ERR_SUCCESS && path2Vertex->Type == kmvtRead && kmer_vertex_out_degree(path2Vertex) > 1 && kmer_vertex_in_degree(path2Vertex) == 1) {
						GEN_ARRAY_size_t tmpw2;
						POINTER_ARRAY_READ_INFO tmpRP;
						POINTER_ARRAY_KMER_EDGE tmpE;

						pointer_array_init_KMER_EDGE(&tmpE, 140);
						pointer_array_init_READ_INFO(&tmpRP, 140);
						dym_array_init_size_t(&tmpw2, 140);
						for (size_t i = 0; i < kmer_vertex_out_degree(path2Vertex); ++i) {
							char *tmpSeq = NULL;
							size_t tmpSeqLen = 0;
							PKMER_EDGE succEdge = kmer_vertex_get_succ_edge(path2Vertex, i);

							pointer_array_clear_KMER_EDGE(&tmpE);
							pointer_array_clear_READ_INFO(&tmpRP);
							dym_array_clear_size_t(&tmpw2);
							ret = _follow_line(succEdge, &tmpSeq, &tmpSeqLen, &succEdge, &tmpw2, &tmpRP, &tmpE);
							if (ret == ERR_SUCCESS) {
								if (path1Vertex == succEdge->Dest) {
									path2Start = succEdge;
									path2Vertex = succEdge->Dest;
									rs_storage_add_seq(&s2, tmpSeq, tmpSeqLen);
									utils_free(tmpSeq);
									dym_array_push_back_array_size_t(&w2, &tmpw2);
									pointer_array_push_back_array_READ_INFO(&rp2, &tmpRP);
									pointer_array_push_back_array_KMER_EDGE(&es2, &tmpE);
									break;
								} else utils_free(tmpSeq);
							}

							if (ret != ERR_SUCCESS)
								break;
						}

						dym_array_finit_size_t(&tmpw2);
						pointer_array_finit_READ_INFO(&tmpRP);
						pointer_array_finit_KMER_EDGE(&tmpE);
					}
					
					if (path2Vertex->Type == kmvtRefSeqMiddle && !path2Vertex->Helper)
						rs_storage_remove(&s2, 1);

					kh_destroy(es, table);
				}

				if (path2Vertex->Type == kmvtRefSeqMiddle) {
					if (ret == ERR_SUCCESS && path1Vertex == path2Vertex) {
						PKMER_EDGE e = kmer_graph_get_edge(Graph, &v->KMer, &path1Vertex->KMer);

						if (e == NULL || e == path1Start || e == path2Start) {
							for (size_t i = 0; i < pointer_array_size(&es1); ++i)
								pointer_array_clear_READ_INFO(&(es1.Data[i]->ReadIndices));
							
							_create_variants(kmer_graph_get_kmer_size(Graph), CHrom, v->AbsPos, s1.Sequence, s1.ValidLength, s2.Sequence, s2.ValidLength, &w1, &w2, &rp1, &rp2, NULL, Options, VCArray);
							kmer_graph_delete_edge(Graph, path1Start);
							kmer_graph_delete_edge(Graph, path2Start);
							ret = kmer_graph_add_edge_ex(Graph, v, path1Vertex, kmetVariant, &e);
							if (ret == ERR_SUCCESS) {
								char *tmpSeq = NULL;

								v = _get_refseq_or_variant_edge(path1Vertex)->Dest;
								edgeCreated = TRUE;
								ret = rs_storage_create_string_with_offset(&s1, 1, &tmpSeq);
								if (ret == ERR_SUCCESS) {
									kmer_edge_add_seq(e, kmetReference, tmpSeq, s1.ValidLength - 1);
									e->Seq1Weight = weight1;
									ret = dym_array_push_back_array_size_t(&e->Weights,  &w1);
									if (ret == ERR_SUCCESS)
										ret = pointer_array_push_back_array_READ_INFO(&e->ReadIndices, &rp1);
								}
							}
						}

						*Changed = TRUE;
					}

					if (ret == ERR_TWO_READ_SEQUENCES)
						ret = ERR_SUCCESS;
				}

				rs_storage_reset(&s2);
				dym_array_clear_size_t(&w2);
				pointer_array_clear_READ_INFO(&rp2);
				if (ret != ERR_SUCCESS || edgeCreated)
					break;
			}
		}

		if (ret != ERR_SUCCESS)
			break;

		if (!edgeCreated)
			v = _get_refseq_or_variant_edge(v)->Dest;
	}

	size_t dummy = 0;
	pointer_array_finit_READ_INFO(&rp2);
	pointer_array_finit_READ_INFO(&rp1);
	dym_array_finit_size_t(&w2);
	dym_array_finit_size_t(&w1);
	rs_storage_finit(&s2);
	rs_storage_finit(&s1);
	pointer_array_finit_KMER_EDGE(&es2);
	pointer_array_finit_KMER_EDGE(&es1);

	kmer_graph_delete_trailing_things(Graph, &dummy);
//	kmer_graph_delete_1to1_vertices(Graph);

	return ret;
}


PKMER_EDGE kmer_vertex_get_edge_by_base(PKMER_VERTEX Vertex, const char Base)
{
	boolean found = FALSE;
	PKMER_EDGE ret = NULL;

	for (size_t i = 0; i < gen_array_size(&Vertex->Successors); ++i) {
		const KMER *kmer = NULL;
		
		ret = Vertex->Successors.Data[i];
		kmer = &ret->Dest->KMer;
		found = (kmer_get_base(kmer, kmer_get_size(kmer) - 1) == Base);
		if (found)
			break;
	}

	if (!found)
		ret = NULL;

	return ret;
}


void kmer_graph_check_weights(PKMER_GRAPH Graph)
{
	void *iter = NULL;
	PKMER_EDGE e = NULL;

	if (kmer_edge_table_first(Graph->EdgeTable, &iter, &e) == ERR_SUCCESS) {
		do {
			assert(gen_array_size(&e->Weights) == pointer_array_size(&e->ReadIndices));
			assert(pointer_array_size(&e->ReadIndices) == e->SeqLen + (!e->Dest->Helper ? 1 : 0));
		} while (kmer_edge_table_next(Graph->EdgeTable, iter, &iter, &e) == ERR_SUCCESS);
	}

	return;
}

ERR_VALUE kmer_graph_compute_weights(PKMER_GRAPH Graph)
{
	void *iter = NULL;
	PKMER_EDGE e = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (kmer_edge_table_first(Graph->EdgeTable, &iter, &e) == ERR_SUCCESS) {
		do {
			size_t wLen = e->SeqLen + 1;

			if (e->Dest->Helper)
				wLen -= 1;

			assert(e->Seq1Weight == 0);
			assert(e->SeqLen == 0 || e->LongData.LongEdge);
			e->Seq1Weight = read_info_weight(&e->ReadInfo, Graph->QualityTable);
			ret = dym_array_reserve_size_t(&e->Weights, wLen);			
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < wLen; ++i)
					dym_array_push_back_no_alloc_size_t(&e->Weights, e->Seq1Weight);

				ret = pointer_array_reserve_READ_INFO(&e->ReadIndices, wLen);
				if (ret == ERR_SUCCESS) {
					for (size_t i = 0; i < wLen; ++i) {
						PREAD_INFO ri = NULL;

						ret = utils_malloc(sizeof(READ_INFO), &ri);
						if (ret == ERR_SUCCESS) {
							read_info_init(ri);
							pointer_array_push_back_no_alloc_READ_INFO(&e->ReadIndices, ri);
							ret = read_info_assign(ri, &e->ReadInfo.Array);
						}

						if (ret != ERR_SUCCESS)
							break;
					}
				}
			}
		} while (ret == ERR_SUCCESS && kmer_edge_table_next(Graph->EdgeTable, iter, &iter, &e) == ERR_SUCCESS);
	}

	return ret;
}


ERR_VALUE kmer_edge_add_read(PKMER_EDGE Edge, size_t ReadIndex, size_t ReadPosition, uint8_t Quality)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = read_info_add(&Edge->ReadInfo, ReadIndex, ReadPosition, Quality);

	return ret;
}
