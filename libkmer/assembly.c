
#include <malloc.h>
#include <stdlib.h>
#include <inttypes.h>
#include "khash.h"
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "pointer_array.h"
#include "refseq-storage.h"
#include "ssw.h"
#include "assembly.h"




/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

/** Represents one vertex in the helper graph. */
typedef struct _DISTANCE_RECORD {
	/** Minimal distance from the first layer to this vertex. */
	size_t Distance;
	/** Previous vertex on the minimal path. */
	size_t BackIndex;
	/** Corresponding "main" graph vertex. */
	PKMER_VERTEX Vertex;
	/** Index within the layer. */
	size_t Index;
} DISTANCE_RECORD, *PDISTANCE_RECORD;

GEN_ARRAY_TYPEDEF(DISTANCE_RECORD);
GEN_ARRAY_IMPLEMENTATION(DISTANCE_RECORD)


UTILS_TYPED_CALLOC_FUNCTION(PPOINTER_ARRAY_KMER_VERTEX)
UTILS_TYPED_CALLOC_FUNCTION(GEN_ARRAY_DISTANCE_RECORD)


/** @brief
 *  Allocates and initializes storage for the helper graph.
 *
 *  @param VertexSets Vertex sets (M_i) assigned to individual k-mers of the read.
 *  @param NumberOfSets Number of the vertex sets.
 *  @param NumberOfRSSets Number of sets containing reference vertices.
 *  @param hGraph Receives the initialized helper graph. 
 */
static ERR_VALUE _helper_graph_init(PPOINTER_ARRAY_KMER_VERTEX *VertexSets, const size_t NumberOfSets, const size_t NumberOfRSSets, PGEN_ARRAY_DISTANCE_RECORD *hGraph)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PGEN_ARRAY_DISTANCE_RECORD tmpDistances = NULL;

	ret = utils_calloc_GEN_ARRAY_DISTANCE_RECORD(NumberOfRSSets, &tmpDistances);
	if (ret == ERR_SUCCESS) {
		PGEN_ARRAY_DISTANCE_RECORD currentD = tmpDistances;
		PPOINTER_ARRAY_KMER_VERTEX currentV = VertexSets[0];
		DISTANCE_RECORD dr;

		for (size_t i = 0; i < NumberOfSets; ++i) {
			currentV = VertexSets[i];
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
		
			*hGraph = tmpDistances;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpDistances);
	}

	return ret;
}


/** @brief
 *  Destroys a given helper graph.
 *
 *  @param hGraph The helper graph vertices.
 *  @param LayerCount Number of layers in the graph.
 */
static void _helper_graph_finit(PGEN_ARRAY_DISTANCE_RECORD hGraph, const size_t LayerCount)
{
	PGEN_ARRAY_DISTANCE_RECORD currentD = hGraph;

	for (size_t i = 0; i < LayerCount; ++i) {
		dym_array_finit_DISTANCE_RECORD(currentD);
		++currentD;
	}

	utils_free(hGraph);

	return;
}


static size_t _missing_edge_price(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, const KMER_VERTEX *Source, const KMER_VERTEX *Dest)
{
	size_t ret = Options->MissingEdgePenalty;
	const KMER_EDGE *e = NULL;

	e = kmer_graph_get_edge(Graph, &Source->KMer, &Dest->KMer);
	if (e != NULL) {
		if (e->Type == kmetReference || read_info_weight(&e->ReadInfo, Graph->QualityTable) > 100*Options->ReadThreshold)
			ret = 0;
	}


	return ret;
}


/** @brief
 *  Finds the shortest path from the first to the last layer of the given helper graph.
 *
 *  @param Options Algorithm parameters (useful for edge weight computation).
 *  @param Graph The main graph.
 *  @param hGraph The helper graph.
 *  @param LayerCount Number of layers within the helper graph.
 *  @param VertexSets Vertex sets assigned to individual k-mers of the read.
 *  @param VertexSetCount Number of the vertex sets.
 *
 *  @remark
 *  The main graph is used when computing weights of edges inside the helper graph,
 *  and to detect when the reads attempts to go backwards in the reference.
 */
static void _helper_graph_build(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, PGEN_ARRAY_DISTANCE_RECORD hGraph, const size_t LayerCount, PPOINTER_ARRAY_KMER_VERTEX *VertexSets, const size_t VertexSetCount)
{
	PGEN_ARRAY_DISTANCE_RECORD currentD = hGraph;
	PGEN_ARRAY_DISTANCE_RECORD nextD = hGraph + 1;

	for (size_t i = 0; i + 1 < LayerCount; ++i) {
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
					x = VertexSets[cdr->Index + 1]->Data[0];
					y = VertexSets[ndr->Index - 1]->Data[0];
					distance += _missing_edge_price(Options, Graph, u, x);
					distance += _missing_edge_price(Options, Graph, y, v);
					for (size_t l = cdr->Index + 1; l < ndr->Index - 1; ++l) {
						x = VertexSets[l]->Data[0];
						y = VertexSets[l + 1]->Data[0];
						distance += _missing_edge_price(Options, Graph, x, y);
					}
				} else distance += _missing_edge_price(Options, Graph, u, v);

				if (u->RefSeqPosition >= v->RefSeqPosition)
					distance += (ndr->Index - cdr->Index)*Options->BackwardRefseqPenalty;
				else if (v->RefSeqPosition - u->RefSeqPosition >= kmer_graph_get_kmer_size(Graph) + 30)
					distance += (v->RefSeqPosition - u->RefSeqPosition)*Options->BackwardRefseqPenalty;

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


/** @brief
 *  Retrieves the shortest path found previously in the helper graph. This is equivalent
 *  to selection of representatives in all the M_i sets.
 *
 *  @param VertexSets Vertex sets associated to individual k-mers of the read.
 *  @param VertexSetCount Number of the vertex sets.
 *  @param hGraph The helper graph.
 *  @param LayerCount Number of layers in the helper graph.
 *  @param Path Receives the shortest path as a sequence of vertices representing the
 *  read in the graph.
 */
static void _shortest_path(PPOINTER_ARRAY_KMER_VERTEX *VertexSets, const size_t VertexSetCount, PGEN_ARRAY_DISTANCE_RECORD hGraph, const size_t LayersCount, PKMER_VERTEX *Path)
{
	const DISTANCE_RECORD *tmpRec = NULL;
	PGEN_ARRAY_DISTANCE_RECORD currentD = hGraph + LayersCount - 1;
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
	for (size_t j = tmpRec->Index + 1; j < VertexSetCount; ++j)
		Path[j] = VertexSets[j]->Data[0];

	size_t lastRSIndex = tmpRec->Index;

	currentD = hGraph + LayersCount - 2;
	for (size_t i = 0; i + 1 < LayersCount; ++i) {
		tmpRec = &currentD->Data[tmpRec->BackIndex];

		Path[tmpRec->Index] = tmpRec->Vertex;
		for (size_t j = tmpRec->Index + 1; j < lastRSIndex; ++j)
			Path[j] = VertexSets[j]->Data[0];

		lastRSIndex = tmpRec->Index;
		--currentD;
	}

	for (size_t j = 0; j < lastRSIndex; ++j)
		Path[j] = VertexSets[j]->Data[0];

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

		ret = _helper_graph_init(Vertices, NumberOfVertices, rsVertexCount, &distances);
		if (ret == ERR_SUCCESS) {
			PKMER_VERTEX *tmpResult = NULL;

			_helper_graph_build(Options, Graph, distances, rsVertexCount, Vertices, NumberOfVertices);
			ret = utils_calloc_PKMER_VERTEX(NumberOfVertices, &tmpResult);
			if (ret == ERR_SUCCESS) {
				_shortest_path(Vertices, NumberOfVertices, distances, rsVertexCount, tmpResult);
				*Result = tmpResult;
				*ResultLength = NumberOfVertices;
			}

			_helper_graph_finit(distances, rsVertexCount);
		}
	} else {
		PKMER_VERTEX *tmpResult = NULL;

		ret = utils_calloc_PKMER_VERTEX(NumberOfVertices, &tmpResult);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < NumberOfVertices; ++i)
				tmpResult[i] = Vertices[i]->Data[0];

			*Result = tmpResult;
			*ResultLength = NumberOfVertices;
		}
	}

	return ret;
}


/** @brief
 *  Assign a vertex set to a given k-mer.
 *
 */
static ERR_VALUE _assign_vertice_set_to_kmer(PKMER_GRAPH Graph, const KMER *KMer, PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t Index, const PARSE_OPTIONS *Options, size_t *SetSize)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	Vertices[Index] = NULL;
	ret = kmer_graph_get_vertices(Graph, KMer, Vertices + Index);
	if (ret == ERR_SUCCESS) {
		*SetSize += pointer_array_size(Vertices[Index]);
	} else if (ret == ERR_NOT_FOUND) {
		PKMER_VERTEX v = NULL;

		ret = kmer_graph_add_vertex_ex(Graph, KMer, kmvtRead, &v);
		if (ret == ERR_SUCCESS) {
			kmer_graph_get_vertices(Graph, KMer, Vertices + Index);
			*SetSize += 1;
		}
	}

	return ret;
}


static ERR_VALUE _assign_vertice_sets_to_kmers(PKMER_GRAPH Graph, const ONE_READ *Read, PPOINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfSets, const PARSE_OPTIONS *Options, boolean *Linear)
{
	size_t count = 0;
	PKMER kmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const uint32_t KMerSize = kmer_graph_get_kmer_size(Graph);

	KMER_STACK_ALLOC(kmer, 0, KMerSize, Read->ReadSequence);
	ret = _assign_vertice_set_to_kmer(Graph, kmer, Vertices, 0, Options, &count);
	if (ret == ERR_SUCCESS) { {
		const KMER_VERTEX *tmp = Vertices[0]->Data[0];

		if (tmp->Type == kmvtRead && !tmp->ReadStartAllowed) {
			PKMER dummyKMer = NULL;

			count--;
			KMER_STACK_ALLOC(dummyKMer, 0, KMerSize, NULL);
			kmer_init_by_base(dummyKMer, KMerSize, 'H');
			ret = _assign_vertice_set_to_kmer(Graph, dummyKMer, Vertices, 0, Options, &count);
		}
	}
		
		if (ret == ERR_SUCCESS) {
			int i = 1;
			int optimizationIndex = -1;

			kmer_advance(KMerSize, kmer, Read->ReadSequence[KMerSize]);
			while (i < (int)NumberOfSets) {
				ret = _assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
				if (ret == ERR_SUCCESS && Options->OptimizeShortVariants && 
					NumberOfSets - i > 16 && optimizationIndex < i - 1) {
					PPOINTER_ARRAY_KMER_VERTEX currVertices = Vertices[i];
					PPOINTER_ARRAY_KMER_VERTEX prevVertices = Vertices[i - 1];

					const KMER_VERTEX *rsv = prevVertices->Data[0];
					PKMER_VERTEX rev = currVertices->Data[0];

					boolean isOk = TRUE;
					if (pointer_array_size(currVertices) > 0) {
						for (size_t i = 0; i < pointer_array_size(currVertices); ++i) {
							isOk = (currVertices->Data[i]->RefSeqPosition - 1 != rsv->RefSeqPosition);
							if (!isOk)
								break;
						}
					}

					if (isOk && pointer_array_size(prevVertices) == 1 &&
						rsv->Type == kmvtRefSeqMiddle &&
						rsv->RefSeqPosition < Options->RegionLength - 1 &&
						rev->RefSeqPosition - 1 != rsv->RefSeqPosition) {
						size_t opStringSize = 0;
						char *opString = NULL;
						const char *ref = Options->Reference + rsv->RefSeqPosition + 1;
						const char *alt = Read->ReadSequence + KMerSize + i - 1;

						ret = ssw_clever(ref, 12, alt, 12, 2, -1, -1, &opString, &opStringSize);
						if (ret == ERR_SUCCESS) {
							boolean oneType = TRUE;
							char typeChar = opString[0];
							size_t matchIndex = 0;
							size_t MCount = 0;

							while (matchIndex < opStringSize && opString[matchIndex] != 'M' && opString[matchIndex] == typeChar)
								++matchIndex;

							oneType = (matchIndex < opStringSize && opString[matchIndex] == 'M');
							while (opString[matchIndex + MCount] == 'M')
								++MCount;

							if (oneType && matchIndex <= 4 && MCount >= 5) {
								switch (typeChar) {
									case 'I': {
										for (uint32_t k = 0; k < matchIndex; ++k)
											kmer_set_base(KMerSize, kmer, k, 'H');

										_assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
										for (size_t j = 0; j < matchIndex - 1; ++j) {
											kmer_advance(KMerSize, kmer, Read->ReadSequence[i + KMerSize]);
											++i;
											_assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
										}

										kmer_init_from_kmer(kmer, KMerSize, &rsv->KMer);
										rev->ReadStartAllowed = FALSE;
										optimizationIndex = i;
									} break;
									case 'D': {
										if (rsv->RefSeqPosition + matchIndex + 1 < Options->RegionLength) {
											kmer_init_from_kmer(kmer, KMerSize, &Graph->RefVertices.Data[rsv->RefSeqPosition + matchIndex + 1]->KMer);
											--count;
											_assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
											optimizationIndex = i;
										}
									} break;
									case 'X': {
										if (rsv->RefSeqPosition + matchIndex < Options->RegionLength) {
											for (uint32_t k = 0; k < matchIndex; ++k)
												kmer_set_base(KMerSize, kmer, k, 'H');

											_assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
											for (size_t j = 0; j < matchIndex - 1; ++j) {
												kmer_advance(KMerSize, kmer, Read->ReadSequence[i + KMerSize]);
												++i;
												_assign_vertice_set_to_kmer(Graph, kmer, Vertices, i, Options, &count);
											}

											kmer_init_from_kmer(kmer, KMerSize, &Graph->RefVertices.Data[rsv->RefSeqPosition + matchIndex]->KMer);
											rev->ReadStartAllowed = FALSE;
											optimizationIndex = i;
										}
									} break;
									default:
										break;
								}
							}

							utils_free(opString);
						}
					}
				}

				if (ret != ERR_SUCCESS)
					break;

				if (i < (int)NumberOfSets - 1)
					kmer_advance(KMerSize, kmer, Read->ReadSequence[i + KMerSize]);
				
				++i;
			}
		}
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
		ret = kmer_edge_add_read(e, ReadIndex, ReadPosition, ReadQuality);

	*Edge = e;

	return ret;
}


static ERR_VALUE _create_long_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t StartIndex, const size_t EndIndex, const ONE_READ *Read, const size_t ReadIndex, const EKMerEdgeType Type, PKMER_EDGE *NewEdge)
{
	PKMER_EDGE e = NULL;
	size_t rsLen = EndIndex - StartIndex - 1;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
		
	assert(U->Helper);
	assert(V->Helper);
	ret = kmer_graph_add_edge_ex(Graph, U, V, Type, &e);
	if (ret == ERR_SUCCESS) {
		char *rs = NULL;

		ret = utils_calloc_char(rsLen + 1, &rs);
		if (ret == ERR_SUCCESS) {
			memcpy(rs, Read->ReadSequence + StartIndex, rsLen*sizeof(char));
			rs[rsLen] = '\0';
			kmer_edge_add_seq(e, Type, rs, rsLen);
		}
	} else if (ret == ERR_ALREADY_EXISTS) {
		if (rsLen == e->SeqLen && memcmp(Read->ReadSequence + StartIndex, e->Seq, rsLen*sizeof(char)) == 0) {
			ret = ERR_SUCCESS;
		} else {
//			printf("ONE EDGE - TWO SEQS\n");
//			exit(0);
			ret = ERR_SUCCESS;
		}
	}

	if (ret == ERR_SUCCESS) {
		ret = kmer_edge_add_read(e, ReadIndex, Read->Offset + StartIndex, Read->Quality[StartIndex]);
		*NewEdge = e;
	}

	return ret;
}


static ERR_VALUE _create_short_read_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, const size_t NumberOfVertices, const ONE_READ *Read, const size_t ReadIndex, PKMER_EDGE **EdgePath)
{
	PKMER_EDGE *tmpEdgePath = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t seqIndex = kmerSize;

	ret = ERR_SUCCESS;
	if (NumberOfVertices > 1) {
		ret = utils_calloc_PKMER_EDGE(NumberOfVertices - 1, &tmpEdgePath);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
				PKMER_VERTEX v = Vertices[i];
				PKMER_VERTEX w = Vertices[i + 1];

				ret = _create_short_edge(Graph, v, w, ReadIndex, Read->Offset + seqIndex, Read->Quality[seqIndex], tmpEdgePath + i);
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
		ret = utils_calloc_uint8_t(NumberOfVertices, &tmpFlags);
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


/** @brief
 *  Connect possible bubbles and attempts to bypass repeats.
 *
 *  @param Graph The de Bruijn graph.
 *  @param Vertices Vertices representing the read.
 *  @param Edges Edges connecting the read vertices.
 *  @param LongFlags For each read edge, defines a bit field indicating whether a
 *  connecting edge should start or end there.
 *  @param NumberOfvertices Number of read vertices.
 *  @param Read The read.
 *  @param ReadIndex The read index.
 *  @param PairArray Receives information about the created connecting edges.
 *
 *  @remark
 *  The LongFlags data are initialized by the _mark_long_edge_flags routine.
 */
static ERR_VALUE _create_long_read_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, PKMER_EDGE *Edges, const uint8_t *LongFlags, const size_t NumberOfVertices, const ONE_READ *Read, const size_t ReadIndex, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
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
	boolean permitLongEdge = FALSE;

	ret = ERR_SUCCESS;
	if (NumberOfVertices > 0) {
		for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
			PKMER_VERTEX v = Vertices[i];
			PKMER_VERTEX w = Vertices[i + 1];

			if (LongFlags[i + 1] & READ_EDGE_FLAG_LONG_END) {
				if (readGapStart != (size_t)-1) {
					PKMER_EDGE connectingEdge = NULL;

					ret = _create_long_edge(Graph, gapStartE->Source, w, readGapSeqStart, seqIndex + 1, Read, ReadIndex, longEdgeType, &connectingEdge);
					if (ret == ERR_SUCCESS) {
						KMER_EDGE_PAIR p;

						p.U = gapStartE;
						p.V = Edges[i];
						p.ReadDistance = seqIndex - readGapSeqStart - 1;
						p.ConnectingEdge = connectingEdge;
						connectingEdge->LongData.LongEdge = TRUE;
						if (p.U->Dest->Type == kmvtRefSeqMiddle) {
							p.ConnectingEdge->Source->LongEdgeAllowed |= (permitLongEdge);
							connectingEdge->LongData.RefSeqEnd = p.U->Dest->RefSeqPosition;
							connectingEdge->LongData.RefSeqStart = p.V->Source->RefSeqPosition + 1;
						}

						if (!_edge_pair_exists(PairArray, &p)) {
							p.EdgeCount = i - readGapStart - 1;
							ret = utils_calloc_PKMER_EDGE(p.EdgeCount, &p.Edges);
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
				if (w->Type == kmvtRefSeqMiddle && (lastRefSeqPos > w->RefSeqPosition || (lastRefSeqPos != 0 && lastRefSeqPos + kmer_graph_get_kmer_size(Graph) + 30 < w->RefSeqPosition)))
					permitLongEdge = TRUE;
			}

			if (ret != ERR_SUCCESS)
				break;

			if (v->Type == kmvtRefSeqMiddle)
				lastRefSeqPos = v->RefSeqPosition;

			if (!w->Helper)
				++seqIndex;
		}
	}

	return ret;
}


/** @brief
 *  Inserts helper vertices into the sequence of vertices and edges of a read.
 *
 *  @param Graph The main de Bruijn graph.
 *  @param pVertices Address of the read vertex array, will be updated by helper ones.
 *  @param pEdges Address of the read edge array, will be updated by edges connecting the
 *  helper vertices.
 *  @param pNumberOfVertices Number of vertices representing the read. Will be updated to reflect
 *  addition of the helper vertices.
 */
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


static ERR_VALUE _produce_single_path(const PARSE_OPTIONS *Options, PKMER_GRAPH Graph, const ONE_READ *Read, const size_t MaxNumberOfSets, const boolean CreateDummyVertices, PKMER_VERTEX **Path, size_t *PathLength)
{
	boolean linear = FALSE;
	PPOINTER_ARRAY_KMER_VERTEX *vertices = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc_PPOINTER_ARRAY_KMER_VERTEX(MaxNumberOfSets, &vertices);
	if (ret == ERR_SUCCESS) {
		ret = _assign_vertice_sets_to_kmers(Graph, Read, vertices, MaxNumberOfSets, Options, &linear);
		if (ret == ERR_SUCCESS)
			ret = _find_best_path(Options, Graph, vertices, MaxNumberOfSets, linear, CreateDummyVertices, Path, PathLength);
	
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

	if (Read->ReadSequenceLen > kmerSize) {
		const size_t maxNumberOfVertices = Read->ReadSequenceLen - (kmerSize - 1);
		size_t tmpPathLength = 0;
		PKMER_VERTEX *tmpPath = NULL;

		ret = _produce_single_path(Options, Graph, Read, maxNumberOfVertices, TRUE, &tmpPath, &tmpPathLength);
		if (ret == ERR_SUCCESS) {
			ret = _create_short_read_edges(Graph, tmpPath, tmpPathLength, Read, ReadIndex, EdgePath);
			if (ret == ERR_SUCCESS) {
				*Path = tmpPath;
				*PathLength = tmpPathLength;
			}
		}
	}

	return ret;
}


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE assembly_parse_reference(PASSEMBLY_STATE State)
{
	boolean refRepeats = FALSE;
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	PKMER_GRAPH Graph = State->Graph;
	const PARSE_OPTIONS *ParseOptions = &State->ParseOptions;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX sourceVertex = NULL;
	PKMER_VERTEX destVertex = NULL;
	const char *RefSeq = ParseOptions->Reference;
	const size_t RefSeqLen = ParseOptions->RegionLength;

	ret = pointer_array_reserve_KMER_VERTEX(&Graph->RefVertices, RefSeqLen + 1);
	if (ret == ERR_SUCCESS) {
		Graph->RefVertices.ValidLength = RefSeqLen + 1;
		for (size_t i = 0; i < RefSeqLen; ++i)
			Graph->RefVertices.Data[i] = NULL;

		ret = kmer_alloc(0, kmerSize, RefSeq, &sourceKMer);
	}

	if (ret == ERR_SUCCESS) {
		kmer_back(kmerSize, sourceKMer, 'B');
		ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRefSeqStart, &sourceVertex);
		if (ret == ERR_SUCCESS) {
			kmer_graph_set_starting_vertex(Graph, sourceKMer);
			ret = kmer_alloc(0, kmerSize, RefSeq, &destKMer);
			if (ret == ERR_SUCCESS) {
				kmer_back(kmerSize, destKMer, 'B');
				for (uint32_t i = kmerSize - 1; i < RefSeqLen; ++i) {
					kmer_advance(kmerSize, destKMer, RefSeq[i]);
					kmer_set_number(destKMer, 0);
					ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqMiddle, &destVertex);
					if (ret == ERR_ALREADY_EXISTS) {
						do {
							destVertex->Unique = FALSE;
							refRepeats = TRUE;
							kmer_set_number(destKMer, kmer_get_number(destKMer) + 1);
							ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqMiddle, &destVertex);
						} while (ret == ERR_ALREADY_EXISTS);
					}

					if (ret == ERR_SUCCESS) {
						PKMER_EDGE edge = NULL;

						destVertex->RefSeqPosition = i;
						Graph->RefVertices.Data[i] = destVertex;
						destVertex->AbsPos = ParseOptions->RegionStart + destVertex->RefSeqPosition;
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
						sourceVertex->RefEdge = edge;
						sourceVertex->RefVarEdge = edge;
						kmer_advance(kmerSize, sourceKMer, RefSeq[i]);
						kmer_set_number(sourceKMer, kmer_get_number(destKMer));
						sourceVertex = destVertex;
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				if (ret == ERR_SUCCESS) {
					kmer_advance(kmerSize, destKMer, 'B');
					kmer_set_number(destKMer, 0);
					ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRefSeqEnd, &destVertex);
					if (ret == ERR_SUCCESS) {
						PKMER_EDGE edge = NULL;

						destVertex->RefSeqPosition = (uint32_t)RefSeqLen;
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
						sourceVertex->RefEdge = edge;
						sourceVertex->RefVarEdge = edge;
						if (ret == ERR_SUCCESS)
							kmer_graph_set_ending_vertex(Graph, destKMer);
					}
				}

				kmer_free(destKMer);
			}
		}
	
		kmer_free(sourceKMer);
	}

	if (ret == ERR_SUCCESS && refRepeats)
		ret = ERR_REF_REPEATS;

	return ret;
}


ERR_VALUE assembly_parse_reads(PASSEMBLY_STATE State)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead = NULL;
	PKMER_GRAPH Graph = State->Graph;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	const PARSE_OPTIONS *Options = &State->ParseOptions;
	const size_t ReadCount = State->ReadCount;
	PONE_READ Reads = State->Reads;
	PKMER_VERTEX **paths = State->Paths;
	PKMER_EDGE **edgePaths = State->EdgePaths;
	uint8_t **flagPaths = State->FlagPaths;
	size_t *pathLengths = State->PathLengths;

//	_sort_reads(Graph, Reads, ReadCount);
	currentRead = Reads;
	for (size_t i = 0; i < ReadCount; ++i) {
		ret = _kmer_graph_parse_read_v2(Options, Graph, currentRead, currentRead->ReadIndex, paths + i, pathLengths + i, edgePaths + i);
		if (ret != ERR_SUCCESS)
			break;
						
		++currentRead;
	}

	return ret;
}


ERR_VALUE assembly_add_helper_vertices(PASSEMBLY_STATE State)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead;
	const size_t ReadCount = State->ReadCount;
	PKMER_GRAPH Graph = State->Graph;
	PKMER_VERTEX **paths = State->Paths;
	PKMER_EDGE **edgePaths = State->EdgePaths;
	size_t *pathLengths = State->PathLengths;

	if (State->ParseOptions.HelperVertices) {
		currentRead = State->Reads;
		for (size_t i = 0; i < ReadCount; ++i) {
			ret = _add_read_helper_vertices(Graph, paths + i, edgePaths + i, pathLengths + i);
			if (ret != ERR_SUCCESS)
				break;

			++currentRead;
		}
	}

	return ret;
}


ERR_VALUE assembly_create_long_edges(PASSEMBLY_STATE State, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead;
	const size_t ReadCount = State->ReadCount;
	PKMER_GRAPH Graph = State->Graph;
	PKMER_VERTEX **paths = State->Paths;
	PKMER_EDGE **edgePaths = State->EdgePaths;
	size_t *pathLengths = State->PathLengths;
	uint8_t **flagPaths = State->FlagPaths;

	currentRead = State->Reads;
	for (size_t i = 0; i < ReadCount; ++i) {
		ret = _mark_long_edge_flags(&State->ParseOptions, paths[i], pathLengths[i], flagPaths + i);
		if (ret != ERR_SUCCESS)
			break;

		++currentRead;
	}

	if (ret == ERR_SUCCESS) {
		currentRead = State->Reads;
		for (size_t i = 0; i < ReadCount; ++i) {
			ret = _create_long_read_edges(Graph, paths[i], edgePaths[i], flagPaths[i], pathLengths[i], currentRead, currentRead->ReadIndex, PairArray);
			if (ret != ERR_SUCCESS)
				break;

			++currentRead;
		}
	}

	return ret;
}


ERR_VALUE assembly_variants_to_edges(PASSEMBLY_STATE State, const GEN_ARRAY_VARIANT_CALL *VCArray)
{
	PKMER_EDGE e = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	e = _get_refseq_or_variant_edge(State->Graph->StartingVertex);
	while (ret == ERR_SUCCESS && e != NULL) {
		if (e->Type == kmetVariant) {
			const uint64_t startPos = e->Source->AbsPos + 1;
			const uint64_t endPos = e->Dest->AbsPos + 1;
			const VARIANT_CALL *vc = VCArray->Data;

			for (size_t i = 0; i < gen_array_size(VCArray); ++i) {
				if (startPos <= vc->Pos && vc->Pos < endPos) {
					ret = pointer_array_push_back_VARIANT_CALL(&e->VCs, vc);
					if (ret != ERR_SUCCESS)
						break;
				}

				++vc;
			}
		}

		e = _get_refseq_or_variant_edge(e->Dest);
	}

	return ret;
}



ERR_VALUE assembly_state_init(PKMER_GRAPH Graph, const PARSE_OPTIONS *ParseOptions, PONE_READ Reads, size_t ReadCount, PASSEMBLY_STATE State)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	State->Graph = Graph;
	State->Reads = Reads;
	State->ReadCount = ReadCount;
	State->ParseOptions = *ParseOptions;
	ret = utils_calloc_PPKMER_VERTEX(ReadCount, &State->Paths);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc_size_t(ReadCount, &State->PathLengths);
		if (ret == ERR_SUCCESS) {
			ret = utils_calloc_PPKMER_EDGE(ReadCount, &State->EdgePaths);
			if (ret == ERR_SUCCESS) {
				ret = utils_calloc_puint8_t(ReadCount, &State->FlagPaths);
				if (ret == ERR_SUCCESS) {
					for (size_t i = 0; i < State->ReadCount; ++i) {
						State->Paths[i] = NULL;
						State->PathLengths[i] = 0;
						State->EdgePaths[i] = NULL;
						State->FlagPaths[i] = NULL;
					}
				}

				if (ret != ERR_SUCCESS) {
					for (size_t i = 0; i < State->ReadCount; ++i) {
						if (State->EdgePaths[i] != NULL)
							utils_free(State->EdgePaths[i]);
					}

					utils_free(State->EdgePaths);
				}
			}

			if (ret != ERR_SUCCESS)
				utils_free(State->PathLengths);
		}

		if (ret != ERR_SUCCESS) {
			for (size_t i = 0; i < State->ReadCount; ++i) {
				if (State->Paths[i] != NULL)
					utils_free(State->Paths[i]);
			}

			utils_free(State->Paths);
		}
	}

	return ret;
}


void assembly_state_finit(PASSEMBLY_STATE State)
{

	for (size_t i = 0; i < State->ReadCount; ++i) {
		if (State->FlagPaths[i] != NULL)
			utils_free(State->FlagPaths[i]);
	}

	utils_free(State->FlagPaths);
	for (size_t i = 0; i < State->ReadCount; ++i) {
		if (State->EdgePaths[i] != NULL)
			utils_free(State->EdgePaths[i]);
	}

	utils_free(State->EdgePaths);
	utils_free(State->PathLengths);

	for (size_t i = 0; i < State->ReadCount; ++i) {
		if (State->Paths[i] != NULL)
			utils_free(State->Paths[i]);
	}

	utils_free(State->Paths);

	return;
}
