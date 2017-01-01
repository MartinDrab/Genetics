
#include <malloc.h>
#include <stdlib.h>
#include "khash.h"
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "refseq-storage.h"
#include "assembly.h"




/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

typedef struct _DISTANCE_RECORD {
	size_t Distance;
	size_t BackIndex;
	PKMER_VERTEX Vertex;
	size_t Index;
} DISTANCE_RECORD, *PDISTANCE_RECORD;

GEN_ARRAY_TYPEDEF(DISTANCE_RECORD);
GEN_ARRAY_IMPLEMENTATION(DISTANCE_RECORD)


static ERR_VALUE _init_distance_arrays(PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfVertices, const size_t NumberOfRSVertices, PGEN_ARRAY_DISTANCE_RECORD *Distances)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PGEN_ARRAY_DISTANCE_RECORD tmpDistances = NULL;

	ret = utils_calloc(NumberOfRSVertices, sizeof(GEN_ARRAY_DISTANCE_RECORD), &tmpDistances);
	if (ret == ERR_SUCCESS) {
		PGEN_ARRAY_DISTANCE_RECORD currentD = tmpDistances;
		PPOINTER_ARRAY_KMER_VERTEX currentV = Vertices[0];
		DISTANCE_RECORD dr;

		for (size_t i = 0; i < NumberOfVertices; ++i) {
			currentV = Vertices[i];
			assert(pointer_array_size(currentV) > 0);
			if (currentV->Data[0]->Type == kmvtRefSeqMiddle) {
				dym_array_init_DISTANCE_RECORD(currentD, 140);
				dr.BackIndex = (size_t)-1;
				dr.Distance = 0xffffffff;
				dr.Index = i;
				dr.Vertex = NULL;
				ret = dym_array_reserve_DISTANCE_RECORD(currentD, pointer_array_size(currentV));
				if (ret == ERR_SUCCESS) {
					PKMER_VERTEX *pv = currentV->Data;

					for (size_t j = 0; j < pointer_array_size(currentV); ++j) {
						dr.Vertex = *pv;
						dym_array_push_back_no_alloc_DISTANCE_RECORD(currentD, dr);
						++pv;
					}
				}

				if (ret != ERR_SUCCESS)
					break;
				
				++currentD;
			}
		}

		if (ret == ERR_SUCCESS) {
			currentD = tmpDistances;
			for (size_t i = 0; i < gen_array_size(currentD); ++i)
				currentD->Data[i].Distance = 0;
		
			*Distances = tmpDistances;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpDistances);
	}

	return ret;
}


static void _finit_distance_arrays(PGEN_ARRAY_DISTANCE_RECORD Distances, const size_t DistancesCount)
{
	PGEN_ARRAY_DISTANCE_RECORD currentD = Distances;

	for (size_t i = 0; i < DistancesCount; ++i) {
		dym_array_finit_DISTANCE_RECORD(currentD);
		++currentD;
	}

	utils_free(Distances);

	return;
}


static void _build_distance_graph(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, PGEN_ARRAY_DISTANCE_RECORD Distances, const size_t DistanceCount, PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfVertices)
{
	PGEN_ARRAY_DISTANCE_RECORD currentD = Distances;
	PGEN_ARRAY_DISTANCE_RECORD nextD = Distances + 1;

	for (size_t i = 0; i + 1 < DistanceCount; ++i) {
		const DISTANCE_RECORD *cdr = currentD->Data;

		for (size_t j = 0; j < gen_array_size(currentD); ++j) {
			PDISTANCE_RECORD ndr = nextD->Data;

			for (size_t k = 0; k < gen_array_size(nextD); ++k) {
				size_t distance = 0;
				PKMER_VERTEX u = cdr->Vertex;
				PKMER_VERTEX v = ndr->Vertex;
				PKMER_VERTEX x = NULL;
				PKMER_VERTEX y = NULL;

				if (ndr->Index - cdr->Index > 1) {
					x = Vertices[cdr->Index + 1]->Data[0];
					y = Vertices[ndr->Index - 1]->Data[0];
					if (kmer_graph_get_edge(Graph, &u->KMer, &x->KMer) == NULL)
						distance += Options->MissingEdgePenalty;

					if (kmer_graph_get_edge(Graph, &y->KMer, &v->KMer) == NULL)
						distance += Options->MissingEdgePenalty;
				} else {
					if (kmer_graph_get_edge(Graph, &u->KMer, &v->KMer) == NULL)
						distance += Options->MissingEdgePenalty;
				}

				if (u->RefSeqPosition > v->RefSeqPosition)
					distance += Options->BackwardRefseqPenalty;

				if (cdr->Distance + distance < ndr->Distance) {
					ndr->Distance = cdr->Distance + distance;
					ndr->BackIndex = j;
				}

				++ndr;
			}

			++cdr;
		}

		++currentD;
		++nextD;
	}

	return;
}


static void _shortest_path(PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfVertices, PGEN_ARRAY_DISTANCE_RECORD Distances, const size_t DistanceCount, PKMER_VERTEX *Path)
{
	const DISTANCE_RECORD *tmpRec = NULL;
	PGEN_ARRAY_DISTANCE_RECORD currentD = Distances + DistanceCount - 1;
	const DISTANCE_RECORD *cdr = currentD->Data;
	size_t minDistance = 0xffffffff;

	for (size_t j = 0; j < gen_array_size(currentD); ++j) {
		if (minDistance > cdr->Distance) {
			minDistance = cdr->Distance;
			tmpRec = cdr;
		}

		++cdr;
	}

	Path[tmpRec->Index] = tmpRec->Vertex;
	for (size_t j = tmpRec->Index + 1; j < NumberOfVertices; ++j)
		Path[j] = Vertices[j]->Data[0];

	size_t lastRSIndex = tmpRec->Index;

	currentD = Distances + DistanceCount - 2;
	for (size_t i = 0; i + 1 < DistanceCount; ++i) {
		tmpRec = &currentD->Data[tmpRec->BackIndex];

		Path[tmpRec->Index] = tmpRec->Vertex;
		for (size_t j = tmpRec->Index + 1; j < lastRSIndex; ++j)
			Path[j] = Vertices[j]->Data[0];

		lastRSIndex = tmpRec->Index;
		--currentD;
	}

	for (size_t j = 0; j < lastRSIndex; ++j)
		Path[j] = Vertices[j]->Data[0];

	return;
}


static ERR_VALUE _find_best_path(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfVertices, const boolean Linear, const boolean CreateDummyVertices, PKMER_VERTEX **Result, size_t *ResultLength)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_DISTANCE_RECORD *distances = NULL;
	size_t rsVertexCount = 0;

	if (!Linear) {
		for (size_t i = 0; i < NumberOfVertices; ++i) {
			if (Vertices[i]->Data[0]->Type == kmvtRefSeqMiddle)
				++rsVertexCount;
		}

		ret = _init_distance_arrays(Vertices, NumberOfVertices, rsVertexCount, &distances);
		if (ret == ERR_SUCCESS) {
			PKMER_VERTEX *tmpResult = NULL;

			_build_distance_graph(Options, Graph, distances, rsVertexCount, Vertices, NumberOfVertices);
			ret = utils_calloc(NumberOfVertices, sizeof(PKMER_VERTEX), &tmpResult);
			if (ret == ERR_SUCCESS) {
				_shortest_path(Vertices, NumberOfVertices, distances, rsVertexCount, tmpResult);
				*Result = tmpResult;
				*ResultLength = NumberOfVertices;
			}

			_finit_distance_arrays(distances, rsVertexCount);
		}
	} else {
		PKMER_VERTEX *tmpResult = NULL;

		ret = utils_calloc(NumberOfVertices, sizeof(PKMER_VERTEX), &tmpResult);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < NumberOfVertices; ++i)
				tmpResult[i] = Vertices[i]->Data[0];

			*Result = tmpResult;
			*ResultLength = NumberOfVertices;
		}
	}

	return ret;
}


static ERR_VALUE _assign_vertice_sets_to_kmers(PKMER_GRAPH Graph, const char *ReadSequence, PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfSets, boolean *Linear)
{
	size_t count = 0;
	PKMER kmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const uint32_t KMerSize = kmer_graph_get_kmer_size(Graph);

	KMER_STACK_ALLOC(kmer, 0, KMerSize, ReadSequence);
	for (size_t i = 0; i < NumberOfSets; ++i) {
		Vertices[i] = NULL;
		ret = kmer_graph_get_vertices(Graph, kmer, Vertices + i);
		if (ret == ERR_SUCCESS) {
			count += pointer_array_size(Vertices[i]);
		} else if (ret == ERR_NOT_FOUND) {
			PKMER_VERTEX v = NULL;

			ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
			if (ret == ERR_SUCCESS) {
				kmer_graph_get_vertices(Graph, kmer, Vertices + i);
				count += 1;
			}
		}

		if (ret != ERR_SUCCESS)
			break;

		kmer_advance(kmer, ReadSequence[i + KMerSize]);
	}

	if (ret == ERR_SUCCESS)
		*Linear = (count == NumberOfSets);

	return ret;
}


static ERR_VALUE _create_short_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t ReadIndex, const size_t ReadPosition, const uint8_t ReadQuality, PKMER_EDGE *Edge)
{
	PKMER_EDGE e = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	ret = kmer_graph_add_edge_ex(Graph, U, V, kmetRead, &e);
	if (ret == ERR_ALREADY_EXISTS)
		ret = ERR_SUCCESS;

	if (ret == ERR_SUCCESS)
		ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPosition, ReadQuality);

	*Edge = e;

	return ret;
}


static ERR_VALUE _create_long_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t StartIndex, const size_t EndIndex, const READ_PART *ReadPart, const size_t ReadIndex, const EKMerEdgeType Type, PKMER_EDGE *NewEdge)
{
	PKMER_EDGE e = NULL;
	size_t rsLen = EndIndex - StartIndex - 1;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
		
	assert(U->Helper);
	assert(V->Helper);
	ret = kmer_graph_add_edge_ex(Graph, U, V, Type, &e);
	if (ret == ERR_SUCCESS) {
		char *rs = NULL;

		ret = utils_calloc(rsLen + 1, sizeof(char), &rs);
		if (ret == ERR_SUCCESS) {
			memcpy(rs, ReadPart->ReadSequence + StartIndex, rsLen*sizeof(char));
			rs[rsLen] = '\0';
			kmer_edge_add_seq(e, Type, rs, rsLen);
		}
	} else if (ret == ERR_ALREADY_EXISTS) {
		if (rsLen == e->SeqLen && memcmp(ReadPart->ReadSequence + StartIndex, e->Seq, rsLen*sizeof(char)) == 0) {
			ret = ERR_SUCCESS;
		} else {
			printf("ONE EDGE - TWO SEQS\n");
//			exit(0);
			ret = ERR_SUCCESS;
		}
	}

	if (ret == ERR_SUCCESS) {
		ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + StartIndex, 0);
		*NewEdge = e;
	}

	return ret;
}


static ERR_VALUE _create_short_read_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, const size_t NumberOfVertices, const READ_PART *ReadPart, const size_t ReadIndex, PKMER_EDGE **EdgePath)
{
	PKMER_EDGE *tmpEdgePath = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t seqIndex = kmerSize;

	ret = ERR_SUCCESS;
	if (NumberOfVertices > 1) {
		ret = utils_calloc(NumberOfVertices - 1, sizeof(PKMER_EDGE), &tmpEdgePath);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
				PKMER_VERTEX v = Vertices[i];
				PKMER_VERTEX w = Vertices[i + 1];

				ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex, ReadPart->Quality[seqIndex], tmpEdgePath + i);
				if (ret != ERR_SUCCESS)
					break;

				++seqIndex;
			}

			if (ret == ERR_SUCCESS)
				*EdgePath = tmpEdgePath;

			if (ret != ERR_SUCCESS)
				utils_free(tmpEdgePath);
		}
	}

	return ret;
}


static boolean _edge_pair_exists(const GEN_ARRAY_KMER_EDGE_PAIR *Array, const KMER_EDGE_PAIR *Pair)
{
	boolean ret = FALSE;
	const KMER_EDGE_PAIR *item = Array->Data;
		
	for (size_t i = 0; i < gen_array_size(Array); ++i) {
		if (memcmp(item, Pair, sizeof(KMER_EDGE_PAIR_KEY)) == 0) {
			ret = TRUE;
			break;
		}

		++item;
	}

	return ret;
}

#define READ_EDGE_FLAG_LONG_START		0x1
#define READ_EDGE_FLAG_LONG_END			0x2
#define READ_EDGE_FLAG_LONG_REFSEQ		0x4

static ERR_VALUE _mark_long_edge_flags(const PARSE_OPTIONS *Options, const PKMER_VERTEX *Vertices, const size_t NumberOfVertices, uint8_t **Flags)
{
	size_t edgeStartIndex = (size_t)-1;
	uint8_t *tmpFlags = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (NumberOfVertices >= 2) {
		ret = utils_calloc(NumberOfVertices, sizeof(uint8_t), &tmpFlags);
		if (ret == ERR_SUCCESS) {
			memset(tmpFlags, 0, NumberOfVertices*sizeof(uint8_t));
			for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
				const KMER_VERTEX *u = Vertices[i];
				const KMER_VERTEX *v = Vertices[i + 1];
			
				if (u->Type == kmvtRead && v->Type == kmvtRefSeqMiddle) {
					if (edgeStartIndex != (size_t)-1) {
						tmpFlags[i] &= (~READ_EDGE_FLAG_LONG_START);
						edgeStartIndex = (size_t)-1;
					}

					if (Options->ConnectRefSeq) {
						tmpFlags[i] |= (READ_EDGE_FLAG_LONG_START | READ_EDGE_FLAG_LONG_REFSEQ);
						edgeStartIndex = i;
					}
				} else if (u->Type == kmvtRefSeqMiddle && v->Type == kmvtRead) {
					if (edgeStartIndex != (size_t)-1) {
						tmpFlags[i + 1] |= READ_EDGE_FLAG_LONG_END;
						edgeStartIndex = (size_t)-1;
					}
				} else if (u->Type == kmvtRead && v->Type == kmvtRead) {
					if (kmer_vertex_out_degree(u) > 1) {
						if (edgeStartIndex != (size_t)-1) {
							tmpFlags[i + 1] |= READ_EDGE_FLAG_LONG_END;
							edgeStartIndex = (size_t)-1;
						}
					}
					
					if (kmer_vertex_in_degree(v) > 1) {
						if (edgeStartIndex != (size_t)-1) {
							tmpFlags[i] |= READ_EDGE_FLAG_LONG_END;
							edgeStartIndex = (size_t)-1;
						}

						if (Options->ConnectReads) {
							tmpFlags[i] |= READ_EDGE_FLAG_LONG_START;
							edgeStartIndex = i;
						}
					}
				}
			}

			if (ret == ERR_SUCCESS) {
				if (edgeStartIndex != (size_t)-1)
					tmpFlags[edgeStartIndex] &= ~(READ_EDGE_FLAG_LONG_START);

				*Flags = tmpFlags;
			}

			if (ret != ERR_SUCCESS)
				utils_free(tmpFlags);
		}
	}

	return ret;
}


static ERR_VALUE _create_long_read_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, PKMER_EDGE *Edges, const uint8_t *LongFlags, const size_t NumberOfVertices, const READ_PART *ReadPart, const size_t ReadIndex, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	PKMER_EDGE e = NULL;
	size_t readGapStart = (size_t)-1;
	size_t readGapSeqStart = (size_t)-1;
	PKMER_EDGE gapStartE = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR; 
	size_t seqIndex = kmerSize;
	EKMerEdgeType longEdgeType;
	uint32_t lastRefSeqPos = 0;
	uint32_t maxRefSeqPos = 0;
	boolean permitLongEdge = FALSE;

	ret = ERR_SUCCESS;
	if (NumberOfVertices > 0) {
		for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
			PKMER_VERTEX v = Vertices[i];
			PKMER_VERTEX w = Vertices[i + 1];

			if (LongFlags[i + 1] & READ_EDGE_FLAG_LONG_END) {
				if (readGapStart != (size_t)-1) {
					PKMER_EDGE connectingEdge = NULL;

					ret = _create_long_edge(Graph, gapStartE->Source, w, readGapSeqStart, seqIndex + 1, ReadPart, ReadIndex, longEdgeType, &connectingEdge);
					if (ret == ERR_SUCCESS) {
						KMER_EDGE_PAIR p;

						p.U = gapStartE;
						p.V = Edges[i];
						p.ConnectingEdge = connectingEdge;
						if (p.U->Dest->Type == kmvtRefSeqMiddle) {
							p.ConnectingEdge->Source->LongEdgeAllowed |= permitLongEdge;
						}

						p.ReadDistance = seqIndex - readGapSeqStart - 1;
						if (!_edge_pair_exists(PairArray, &p)) {
							p.EdgeCount = i - readGapStart - 1;
							ret = utils_calloc(p.EdgeCount, sizeof(PKMER_EDGE), &p.Edges);
							if (ret == ERR_SUCCESS) {
								memcpy(p.Edges, Edges + readGapStart, p.EdgeCount*sizeof(PKMER_EDGE));
								ret = dym_array_push_back_KMER_EDGE_PAIR(PairArray, p);
								if (ret != ERR_SUCCESS)
									utils_free(p.Edges);
							}
						}
					}

				}

				readGapStart = (size_t)-1;
				readGapSeqStart = (size_t)-1;
			}

			if (LongFlags[i] & READ_EDGE_FLAG_LONG_START) {
				readGapStart = i;
				readGapSeqStart = seqIndex;
				gapStartE = Edges[i];
				longEdgeType = kmetRead;
				if (w->Type == kmvtRefSeqMiddle && (lastRefSeqPos > w->RefSeqPosition))
					permitLongEdge = TRUE;
			}

			if (ret != ERR_SUCCESS)
				break;

			if (v->Type == kmvtRefSeqMiddle) {
				lastRefSeqPos = v->RefSeqPosition;
				if (lastRefSeqPos > maxRefSeqPos)
					maxRefSeqPos = lastRefSeqPos;
			}

			if (!w->Helper)
				++seqIndex;
		}
	}

	return ret;
}


static ERR_VALUE _add_read_helper_vertices(PKMER_GRAPH Graph, PKMER_VERTEX **pVertices, PKMER_EDGE **pEdges, size_t *pNumberOfVertices)
{
	boolean split = FALSE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX *pathVertices = *pVertices;
	POINTER_ARRAY_KMER_VERTEX va;
	POINTER_ARRAY_KMER_EDGE ea;
	const size_t numberOfVertices = *pNumberOfVertices;

	ret = ERR_SUCCESS;
	if (numberOfVertices > 0) {
		pointer_array_init_KMER_VERTEX(&va, 140);
		pointer_array_init_KMER_EDGE(&ea, 140);
		ret = pointer_array_reserve_KMER_VERTEX(&va, 20);
		if (ret == ERR_SUCCESS) {
			ret = pointer_array_reserve_KMER_EDGE(&ea, 20);
			if (ret == ERR_SUCCESS) {
				pointer_array_push_back_no_alloc_KMER_VERTEX(&va, pathVertices[0]);
				if (ret == ERR_SUCCESS) {
					for (size_t i = 0; i < numberOfVertices - 1; ++i) {
						PKMER_VERTEX u = pathVertices[i];
						PKMER_VERTEX v = pathVertices[i + 1];
						PKMER_EDGE e = kmer_graph_get_edge(Graph, &u->KMer, &v->KMer);

						if (!u->Helper && !v->Helper) {
							if (u->Type == kmvtRead && v->Type == kmvtRead) {
								split = (
									kmer_vertex_out_degree(u) > 1 ||
									(kmer_vertex_in_degree(v) > 1)
									);
							}
							else split = (e == NULL || e->Type == kmetRead);
						}

						if (split) {
							PKMER_VERTEX sv = NULL;
							PKMER_EDGE se = NULL;
							PKMER_EDGE de = NULL;

							if (e != NULL)
								ret = kmer_graph_split_edge(Graph, e, &se, &de, &sv);
							else ret = kmer_graph_get_splitted_edge(Graph, u, v, &se, &de, &sv);

							if (ret != ERR_SUCCESS) {
								printf("_add_read_helper_vertices(): %u\n", ret);
								exit(0);
								ret = ERR_SUCCESS;
							}

							if (ret == ERR_SUCCESS)
								ret = pointer_array_push_back_KMER_VERTEX(&va, sv);

							if (ret == ERR_SUCCESS)
								ret = pointer_array_push_back_KMER_VERTEX(&va, v);

							if (ret == ERR_SUCCESS)
								ret = pointer_array_push_back_KMER_EDGE(&ea, se);

							if (ret == ERR_SUCCESS)
								ret = pointer_array_push_back_KMER_EDGE(&ea, de);

							split = FALSE;
						}
						else {
							ret = pointer_array_push_back_KMER_VERTEX(&va, v);
							if (ret == ERR_SUCCESS)
								ret = pointer_array_push_back_KMER_EDGE(&ea, e);
						}

						if (ret != ERR_SUCCESS)
							break;
					}

					if (ret == ERR_SUCCESS) {
						*pNumberOfVertices = pointer_array_size(&va);
						utils_free(*pVertices);
						*pVertices = va.Data;
						utils_free(*pEdges);
						*pEdges = ea.Data;
					}
				}
			}
		}

		if (ret != ERR_SUCCESS) {
			pointer_array_finit_KMER_EDGE(&ea);
			pointer_array_finit_KMER_VERTEX(&va);
		}
	}

	return ret;
}


static ERR_VALUE _produce_single_path(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, const char *ReadSequence, const size_t MaxNumberOfSets, const boolean CreateDummyVertices, PKMER_VERTEX **Path, size_t *PathLength)
{
	boolean linear = FALSE;
	PPOINTER_ARRAY_KMER_VERTEX *vertices = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(MaxNumberOfSets, sizeof(PPOINTER_ARRAY_KMER_VERTEX), &vertices);
	if (ret == ERR_SUCCESS) {
//		for (size_t i = 0; i < MaxNumberOfSets; ++i)
//			pointer_array_init_KMER_VERTEX(vertices + i, 140);

		ret = _assign_vertice_sets_to_kmers(Graph, ReadSequence, vertices, MaxNumberOfSets, &linear);
		if (ret == ERR_SUCCESS)
			ret = _find_best_path(Options, Graph, vertices, MaxNumberOfSets, linear, CreateDummyVertices, Path, PathLength);
	
//		for (size_t i = 0; i < MaxNumberOfSets; ++i)
//			pointer_array_finit_KMER_VERTEX(vertices + i);

		utils_free(vertices);
	}

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read_v2(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex, PKMER_VERTEX **Path, size_t *PathLength, PKMER_EDGE **EdgePath)
{
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	*Path = NULL;
	*PathLength = 0;
	ret = ERR_SUCCESS;
	const READ_PART *part = NULL;

	part = &Read->Part;
	if (part->ReadSequenceLength > kmerSize) {
		const size_t maxNumberOfVertices = part->ReadSequenceLength - (kmerSize - 1);
		size_t tmpPathLength = 0;
		PKMER_VERTEX *tmpPath = NULL;

		ret = _produce_single_path(Options, Graph, part->ReadSequence, maxNumberOfVertices, TRUE, &tmpPath, &tmpPathLength);
		if (ret == ERR_SUCCESS) {
			ret = _create_short_read_edges(Graph, tmpPath, tmpPathLength, part, ReadIndex, EdgePath);
			if (ret == ERR_SUCCESS) {
				*Path = tmpPath;
				*PathLength = tmpPathLength;
			}
		}
	}

	return ret;
}


static void _fix_read(const KMER_GRAPH *Graph, PONE_READ Reads, KMER_VERTEX **Vertices, const size_t NumberOfVertices, const size_t ReadIndex, PKMER_EDGE *Edges, const PARSE_OPTIONS *Options, PGEN_ARRAY_size_t FixedReads)
{
	boolean fixed = FALSE;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	PONE_READ read = Reads + ReadIndex;

	if (read->Part.ReadSequenceLength > kmerSize) {
		const KMER_VERTEX *u = Vertices[0];
		const KMER_VERTEX *v = Vertices[1];
		char *pBase = read->Part.ReadSequence + kmerSize;
		const size_t realThreshold = Options->ReadThreshold * 100;

		fixed = FALSE;
		for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
			const KMER_EDGE *e = Edges[i];
			const size_t origWeight = read_info_weight(&e->ReadInfo, Graph->QualityTable);

			if ((u->Type != kmvtRefSeqMiddle || u->RefSeqPosition > kmerSize) &&
				(v->Type != kmvtRefSeqMiddle || v->RefSeqPosition > kmerSize)) {
				if (read->Part.Quality[kmerSize + i] < 255 && (origWeight <= realThreshold || *pBase == 'N')) {
					size_t maxValue = 0;
					size_t maxIndex = 0;

					for (size_t j = 0; j < kmer_vertex_out_degree(u); ++j) {
						const KMER_EDGE *f = kmer_vertex_get_succ_edge(u, j);

						if (f == e)
							continue;

						const size_t weight = read_info_weight(&f->ReadInfo, Graph->QualityTable);

						if (maxValue <= weight) {
							maxValue = weight;
							maxIndex = j;
						}
					}

					if (maxValue > realThreshold /*|| maxValue > origWeight*/) {
						PKMER_EDGE newPath = kmer_vertex_get_succ_edge(u, maxIndex);

						char newBase = kmer_get_base(&(newPath->Dest->KMer), kmerSize - 1);
						if (*pBase != newBase) {
							const READ_INFO *ri = &e->ReadInfo;
							const READ_INFO_ENTRY *re = read_info_get_entry(ri, 0);
							const char oldBase = *pBase;
							
							if (read->Parent->Quality[kmerSize + i] > 20)
								read->NumberOfFixes++;

							read->Part.Quality[kmerSize + i] = 255;
							read->Part.ReadSequence[kmerSize + i] = newBase;
							dym_array_push_back_size_t(FixedReads, ReadIndex);
							/*
							for (size_t j = 0; j < read_info_get_count(ri); ++j) {
								PONE_READ r = Reads + re->ReadIndex;

								assert(r->Part.ReadSequence[re->ReadPosition] == oldBase);
								r->Part.Quality[re->ReadPosition] = 255;
								r->Part.ReadSequence[re->ReadPosition] = newBase;
								if (!dym_array_contains_size_t(FixedReads, re->ReadIndex))
									dym_array_push_back_size_t(FixedReads, re->ReadIndex);

								++re;
							}
							*/
							break;
						}
					}
				}
			}

			++pBase;
			u = v;
			v = Vertices[i + 2];
		}
	}

	return;
}


static void _remove_read_from_graph(PKMER_GRAPH Graph, const size_t ReadIndex, PKMER_VERTEX *Vertices, const size_t NumberOfVertices, PKMER_EDGE *EdgePath)
{
	PKMER_EDGE e = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);

	for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
		e = EdgePath[i];
		read_info_remove(&e->ReadInfo, ReadIndex, kmerSize + i);
		if (e->Type == kmetRead && read_info_get_count(&e->ReadInfo) == 0)
			kmer_graph_delete_edge(Graph, e);
	}

	utils_free(EdgePath);
	utils_free(Vertices);

	return;
}


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen)
{
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX sourceVertex = NULL;
	PKMER_VERTEX destVertex = NULL;

	ret = kmer_alloc(0, kmerSize, RefSeq, &sourceKMer);
	if (ret == ERR_SUCCESS) {
		kmer_back(sourceKMer, 'B');
		ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRefSeqStart, &sourceVertex);
		if (ret == ERR_SUCCESS) {
			kmer_graph_set_starting_vertex(Graph, sourceKMer);
			ret = kmer_alloc(0, kmerSize, RefSeq, &destKMer);
			if (ret == ERR_SUCCESS) {
				kmer_back(destKMer, 'B');
				for (size_t i = kmerSize - 1; i < RefSeqLen; ++i) {
					kmer_advance(destKMer, RefSeq[i]);
					kmer_set_number(destKMer, 0);
					ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqMiddle, &destVertex);
					if (ret == ERR_ALREADY_EXISTS) {
						do {
							kmer_set_number(destKMer, kmer_get_number(destKMer) + 1);
							ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqMiddle, &destVertex);
						} while (ret == ERR_ALREADY_EXISTS);
					}

					if (ret == ERR_SUCCESS) {
						PKMER_EDGE edge = NULL;

						destVertex->RefSeqPosition = i;
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
						kmer_advance(sourceKMer, RefSeq[i]);
						kmer_set_number(sourceKMer, kmer_get_number(destKMer));
						sourceVertex = destVertex;
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				if (ret == ERR_SUCCESS) {
					kmer_advance(destKMer, 'E');
					kmer_set_number(destKMer, 0);
					ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqEnd, &destVertex);
					if (ret == ERR_SUCCESS) {
						PKMER_EDGE edge = NULL;

						destVertex->RefSeqPosition = RefSeqLen;
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
						if (ret == ERR_SUCCESS)
							kmer_graph_set_ending_vertex(Graph, destKMer);
					}
				}

				kmer_free(destKMer);
			}
		}
	
		kmer_free(sourceKMer);
	}

	return ret;
}


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, PONE_READ Reads, const size_t ReadCount, const size_t Threshold, const PARSE_OPTIONS *Options, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	PKMER_VERTEX **paths = NULL;
	PKMER_EDGE **edgePaths = NULL;
	uint8_t **flagPaths = NULL;
	size_t *pathLengths = 0;

	ret = utils_calloc(ReadCount, sizeof(PKMER_VERTEX *), &paths);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(ReadCount, sizeof(size_t), &pathLengths);
		if (ret == ERR_SUCCESS) {
			ret = utils_calloc(ReadCount, sizeof(PKMER_EDGE *), &edgePaths);
			if (ret == ERR_SUCCESS) {
				ret = utils_calloc(ReadCount, sizeof(uint8_t *), &flagPaths);
				if (ret == ERR_SUCCESS) {
					for (size_t i = 0; i < ReadCount; ++i) {
						paths[i] = NULL;
						pathLengths[i] = 0;
						edgePaths[i] = NULL;
						flagPaths[i] = NULL;
					}

					currentRead = Reads;
					for (size_t i = 0; i < ReadCount; ++i) {
						ret = _kmer_graph_parse_read_v2(Options, Graph, currentRead, currentRead->ReadIndex, paths + i, pathLengths + i, edgePaths + i);
						if (ret != ERR_SUCCESS)
							break;
						
						++currentRead;
					}

					if (ret == ERR_SUCCESS) {						
						currentRead = Reads;
						for (size_t i = 0; i < ReadCount; ++i) {
							if (Options->HelperVertices)
								ret = _add_read_helper_vertices(Graph, paths + i, edgePaths + i, pathLengths + i);
							
							if (ret == ERR_SUCCESS)
								ret = _mark_long_edge_flags(Options, paths[i], pathLengths[i], flagPaths + i);

							if (ret != ERR_SUCCESS)
								break;

							++currentRead;
						}

						if (ret == ERR_SUCCESS) {
							currentRead = Reads;
							for (size_t i = 0; i < ReadCount; ++i) {
								ret = _create_long_read_edges(Graph, paths[i], edgePaths[i], flagPaths[i], pathLengths[i], &currentRead->Part, currentRead->ReadIndex, PairArray);
								if (ret != ERR_SUCCESS)
									break;

								++currentRead;
							}
						}
					}

					for (size_t i = 0; i < ReadCount; ++i) {
						if (flagPaths[i] != NULL)
							utils_free(flagPaths[i]);
					}

					utils_free(flagPaths);
				}

				for (size_t i = 0; i < ReadCount; ++i) {
					if (edgePaths[i] != NULL)
						utils_free(edgePaths[i]);
				}

				utils_free(edgePaths);
			}

			utils_free(pathLengths);
		}

		for (size_t i = 0; i < ReadCount; ++i) {
			if (paths[i] != NULL)
				utils_free(paths[i]);
		}

		utils_free(paths);
	}

	return ret;
}


ERR_VALUE assembly_repair_reads(const KMER_GRAPH_ALLOCATOR *GraphAllocator, const uint32_t KMerSize, PONE_READ Reads, const size_t ReadCount, const char*RefSeq, const size_t RefSeqLen, const PARSE_OPTIONS *ParseOptions)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead = NULL;
	PKMER_VERTEX **paths = NULL;
	PKMER_EDGE **edgePaths = NULL;
	size_t *pathLengths = 0;

	ret = kmer_graph_create(KMerSize, 2500, 6000, &g);
	if (ret == ERR_SUCCESS) {
		if (GraphAllocator != NULL)
			g->Allocator = *GraphAllocator;

		ret = kmer_graph_parse_ref_sequence(g, RefSeq, RefSeqLen);
		if (ret == ERR_SUCCESS) {
			ret = utils_calloc(ReadCount, sizeof(PKMER_VERTEX *), &paths);
			if (ret == ERR_SUCCESS) {
				ret = utils_calloc(ReadCount, sizeof(size_t), &pathLengths);
				if (ret == ERR_SUCCESS) {
					ret = utils_calloc(ReadCount, sizeof(PKMER_EDGE *), &edgePaths);
					if (ret == ERR_SUCCESS) {
						currentRead = Reads;
						for (size_t i = 0; i < ReadCount; ++i) {
							ret = _kmer_graph_parse_read_v2(ParseOptions, g, currentRead, currentRead->ReadIndex, paths + i, pathLengths + i, edgePaths + i);
							if (ret != ERR_SUCCESS)
								break;

							++currentRead;
						}

						if (ret == ERR_SUCCESS) {
							GEN_ARRAY_size_t fixedReads;

							dym_array_init_size_t(&fixedReads, 140);
							do {
								if (gen_array_size(&fixedReads) > 0) {
									const size_t *preadIndex = dym_array_const_item_size_t(&fixedReads, 0);

									for (size_t i = 0; i < gen_array_size(&fixedReads); ++i) {
										PONE_READ r = Reads + *preadIndex;
										
										_remove_read_from_graph(g, r->ReadIndex, paths[*preadIndex], pathLengths[*preadIndex], edgePaths[*preadIndex]);
										paths[*preadIndex] = NULL;
										edgePaths[*preadIndex] = NULL;
										if (100*r->NumberOfFixes / r->Part.ReadSequenceLength < ParseOptions->ReadMaxErrorRate) {
											ret = _kmer_graph_parse_read_v2(ParseOptions, g, Reads + *preadIndex, (Reads + *preadIndex)->ReadIndex, paths + *preadIndex, pathLengths + *preadIndex, edgePaths + *preadIndex);
											++preadIndex;
											if (ret != ERR_SUCCESS)
												break;
										}
									}
								}

								dym_array_clear_size_t(&fixedReads);
								for (size_t i = 0; i < ReadCount; ++i) {
									PONE_READ r = Reads + i;

									if (100 * r->NumberOfFixes / r->Part.ReadSequenceLength < ParseOptions->ReadMaxErrorRate) {
										_fix_read(g, Reads, paths[i], pathLengths[i], i, edgePaths[i], ParseOptions, &fixedReads);
										if (gen_array_size(&fixedReads) > 0)
											break;
									}
								}
							} while (ret == ERR_SUCCESS && gen_array_size(&fixedReads) > 0);

							dym_array_finit_size_t(&fixedReads);
						}

						for (size_t i = 0; i < ReadCount; ++i) {
							if (edgePaths[i] != NULL)
								utils_free(edgePaths[i]);
						}

						utils_free(edgePaths);
					}

					utils_free(pathLengths);
				}

				for (size_t i = 0; i < ReadCount; ++i) {
					if (paths[i] != NULL)
						utils_free(paths[i]);
				}

				utils_free(paths);
			}
		}

		kmer_graph_destroy(g);
	}

	return ret;
}
