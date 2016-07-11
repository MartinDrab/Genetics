
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

static ERR_VALUE _capture_connect_candidates(const KMER_GRAPH *Graph, const KMER_VERTEX **Vertices, const size_t NumberOfEdges, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t firstReadPoint = (size_t)-1;
	const KMER_VERTEX *u = NULL;
	const KMER_VERTEX *v = NULL;
	PKMER_EDGE e = NULL;
	size_t seqIndex = 0;
	size_t firstSeqPoint = (size_t)-1;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < NumberOfEdges; ++i) {
		u = Vertices[i];
		v = Vertices[i + 1];
		e = kmer_graph_get_edge(Graph, &u->KMer, &v->KMer);
		if (e != NULL && e->Type == kmetRead) {
			if (firstReadPoint != (size_t)-1 && e->Source->Type == kmvtRefSeqMiddle) {
				PKMER_EDGE se = kmer_graph_get_edge(Graph, &(Vertices[firstReadPoint])->KMer, &(Vertices[firstReadPoint + 1])->KMer);
				PKMER_EDGE de = kmer_graph_get_edge(Graph, &(Vertices[i])->KMer, &(Vertices[i + 1])->KMer);
				KMER_EDGE_PAIR pair;

				if (se->Dest->Order <= de->Source->Order) {
					pair.U = se;
					pair.V = de;
					pair.ReadDistance = seqIndex - firstSeqPoint - 1;
					assert(pair.U->Type == kmetRead);
					assert(pair.V->Type == kmetRead);
					if (!dym_array_contains_KMER_EDGE_PAIR(PairArray, pair))
						ret = dym_array_push_back_KMER_EDGE_PAIR(PairArray, pair);
				}
			}

			if (e->Dest->Type == kmvtRefSeqMiddle) {
				firstReadPoint = i;
				firstSeqPoint = seqIndex;
			}
		}

		if (ret != ERR_SUCCESS)
			break;
	
		if (v->Type != kmvtDummy)
			++seqIndex;
	}

	return ret;
}


typedef struct _DISTANCE_RECORD {
	size_t Distance;
	size_t BackIndex;
} DISTANCE_RECORD, *PDISTANCE_RECORD;

GEN_ARRAY_TYPEDEF(DISTANCE_RECORD);
GEN_ARRAY_IMPLEMENTATION(DISTANCE_RECORD)

static ERR_VALUE _find_best_path(const KMER_GRAPH *Graph, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfVertices, boolean Linear, PKMER_VERTEX **Result, size_t *ResultLength)
{
	DISTANCE_RECORD dr;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	GEN_ARRAY_DISTANCE_RECORD *distances = NULL;
	POINTER_ARRAY_KMER_VERTEX pathArray;

	pointer_array_init_KMER_VERTEX(&pathArray, 140);
	dr.BackIndex = 0xffffffff;
	dr.Distance = 0xffffffff;
	ret = utils_calloc(NumberOfVertices, sizeof(GEN_ARRAY_DISTANCE_RECORD), &distances);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < NumberOfVertices; ++i) {
			const size_t targetSize = pointer_array_size(Vertices + i);
			
			dym_array_init_DISTANCE_RECORD(distances + i, 140);
			ret = dym_array_reserve_DISTANCE_RECORD(distances + i, targetSize);
			if (ret == ERR_SUCCESS) {
				for (size_t j = 0; j < targetSize; ++j)
					dym_array_push_back_no_alloc_DISTANCE_RECORD(distances + i, dr);
			}

			if (ret != ERR_SUCCESS) {
				for (size_t j = 0; j < i; ++j)
					dym_array_finit_DISTANCE_RECORD(distances + j);
			
				break;
			}
		}

		if (ret == ERR_SUCCESS) {
			POINTER_ARRAY_KMER_VERTEX *currentV = Vertices;
			const POINTER_ARRAY_KMER_VERTEX *nextV = Vertices + 1;
			PGEN_ARRAY_DISTANCE_RECORD currentD = distances;
			PGEN_ARRAY_DISTANCE_RECORD nextD = distances + 1;

			for (size_t i = 0; i < gen_array_size(currentV); ++i)
				currentD->Data[i].Distance = 0;

			for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
				assert(gen_array_size(currentD) == pointer_array_size(currentV));
				assert(gen_array_size(nextD) == pointer_array_size(nextV));
				for (size_t j = 0; j < gen_array_size(currentD); ++j) {
					const KMER_VERTEX *u = *pointer_array_const_item_KMER_VERTEX(currentV, j);
					DISTANCE_RECORD *d2 = dym_array_item_DISTANCE_RECORD(nextD, 0);

					dr = *dym_array_item_DISTANCE_RECORD(currentD, j);
					for (size_t k = 0; k < gen_array_size(nextD); ++k) {
						const KMER_VERTEX *v = *pointer_array_const_item_KMER_VERTEX(nextV, k);
						size_t newDistance = dr.Distance;

						newDistance += ((kmer_graph_get_edge(Graph, &u->KMer, &v->KMer) == NULL) ? 1 : 0);
						if (u->Type == kmvtRefSeqMiddle && v->Type == kmvtRefSeqMiddle &&
							u->RefSeqPosition > v->RefSeqPosition)
							newDistance = dr.Distance + 2;

						if (d2->Distance > newDistance) {
							d2->Distance = newDistance;
							d2->BackIndex = j;
						}

						++d2;
					}
				}

				currentD = nextD;
				++nextD;
				currentV = nextV;
				++nextV;
			}

			if (ret == ERR_SUCCESS) {
				size_t minDistance = 0xffffffff;
				size_t minDistanceIndex = 0;
				size_t minDistanceBackIndex = 0xffffffff;

				currentD = distances + NumberOfVertices - 1;
				currentV = Vertices + NumberOfVertices - 1;
				for (size_t i = 0; i < gen_array_size(currentD); ++i) {
					dr = *dym_array_item_DISTANCE_RECORD(currentD, i);
					if (dr.Distance < minDistance) {
						minDistance = dr.Distance;
						minDistanceBackIndex = dr.BackIndex;
						minDistanceIndex = i;
					}
				}


				ret = pointer_array_push_back_KMER_VERTEX(&pathArray, *pointer_array_const_item_KMER_VERTEX(currentV, minDistanceIndex));
				if (ret == ERR_SUCCESS) {
					dr.BackIndex = minDistanceBackIndex;
					dr.Distance = minDistance;
					for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
						--currentV;
						--currentD;
						ret = pointer_array_push_back_KMER_VERTEX(&pathArray, *pointer_array_const_item_KMER_VERTEX(currentV, dr.BackIndex));
						dr = *dym_array_item_DISTANCE_RECORD(currentD, dr.BackIndex);
						if (ret != ERR_SUCCESS)
							break;
					}
				}
			}

			for (size_t i = 0; i < NumberOfVertices; ++i)
				dym_array_finit_DISTANCE_RECORD(distances + i);
		}

		utils_free(distances);
	}

	if (ret == ERR_SUCCESS) {
		PKMER_VERTEX *tmpResult = NULL;
		size_t tmpResultLen = pointer_array_size(&pathArray);
		PKMER_VERTEX *v = pathArray.Data;

		for (size_t i = pointer_array_size(&pathArray) - 1; i > 0; --i) {
			if (v[i - 1]->Type == kmvtRefSeqMiddle && v[i]->Type == kmvtRefSeqMiddle &&
				kmer_graph_get_edge(Graph, &(v[i]->KMer), &(v[i - 1])->KMer) == NULL)
				++tmpResultLen;			
		}

		ret = utils_calloc(tmpResultLen, sizeof(PKMER_VERTEX), &tmpResult);
		if (ret == ERR_SUCCESS) {
			size_t index = 0;
			
			v = pathArray.Data + pointer_array_size(&pathArray);
			for (size_t i = 0; i < pointer_array_size(&pathArray); ++i) {
				const KMER_VERTEX *src = *v;
				const KMER_VERTEX *dest = *(v - 1);
				
				--v;
				if (i > 0 &&
					src->Type == kmvtRefSeqMiddle && dest->Type == kmvtRefSeqMiddle &&
					kmer_graph_get_edge(Graph, &src->KMer, &dest->KMer) == NULL) {
					ret = kmer_graph_add_dummy_vertex(Graph, &src->KMer, &dest->KMer, &(tmpResult[index]));
					if (ret == ERR_ALREADY_EXISTS)
						ret = ERR_SUCCESS;
					
					if (ret != ERR_SUCCESS)
						break;

					++index;
				}

				tmpResult[index] = dest;
				++index;
			}

			*Result = tmpResult;
			*ResultLength = tmpResultLen;
		}
	}

	pointer_array_finit_KMER_VERTEX(&pathArray);

	return ret;
}


static ERR_VALUE _assign_vertice_sets_to_kmers(PKMER_GRAPH Graph, const READ_PART *Part, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfSets, const boolean RegionStartTouch, boolean *Linear)
{
	size_t count = 0;
	PKMER kmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	uint32_t KMerSize = kmer_graph_get_kmer_size(Graph);
	char *beginning = alloca(KMerSize*sizeof(char));
	size_t plusIndex = (RegionStartTouch) ? 0 : KMerSize;

	memset(beginning, 'B', KMerSize*sizeof(char));
	if (!RegionStartTouch)
		memcpy(beginning, Part->ReadSequence, KMerSize*sizeof(char));

	KMER_STACK_ALLOC(kmer, 0, KMerSize, beginning);
	for (size_t i = 0; i < NumberOfSets; ++i) {
		ret = kmer_graph_get_vertices(Graph, kmer, Vertices + i);
		if (ret == ERR_SUCCESS) {
			if (pointer_array_size(Vertices + i) == 0) {
				PKMER_VERTEX v = NULL;

				ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
				if (ret == ERR_SUCCESS)
					pointer_array_push_back_no_alloc_KMER_VERTEX(Vertices + i, v);
			}

			count += pointer_array_size(Vertices + i);
		}

		if (ret != ERR_SUCCESS)
			break;

		kmer_advance(kmer, Part->ReadSequence[i + plusIndex]);
	}

	if (ret == ERR_SUCCESS)
		*Linear = (count == NumberOfSets);

	return ret;
}


static ERR_VALUE _create_short_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t ReadIndex, const size_t ReadPosition)
{
	PKMER_EDGE e = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	ret = kmer_graph_add_edge_ex(Graph, U, V, 1, 1, kmetRead, &e);
	if (ret == ERR_ALREADY_EXISTS) {
		e->Weight++;
		ret = ERR_SUCCESS;
	}

	if (ret == ERR_SUCCESS)
		ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPosition);

	return ret;
}


static ERR_VALUE _create_long_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t StartIndex, const size_t EndIndex, const size_t ExtraOffset, const READ_PART *ReadPart, const size_t ReadIndex)
{
	char *rs = NULL;
	PKMER_EDGE e = NULL;
	size_t rsLen = EndIndex - StartIndex - 1;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(rsLen + 1, sizeof(char), &rs);
	if (ret == ERR_SUCCESS) {
		memcpy(rs, ReadPart->ReadSequence + StartIndex + ExtraOffset, rsLen*sizeof(char));
		rs[rsLen] = '\0';
		ret = kmer_graph_add_edge_ex(Graph, U, V, 1, 0, kmetRead, &e);
		if (ret == ERR_SUCCESS) {
			e->Seq = rs;
			e->SeqLen = rsLen;
			e->SeqType = kmetRead;
			ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + StartIndex + ExtraOffset);
		} else if (ret == ERR_ALREADY_EXISTS) {
			if (rsLen == e->SeqLen && memcmp(e->Seq, rs, rsLen) == 0) {
				e->Weight++;
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + StartIndex + ExtraOffset);
				if (ret == ERR_SUCCESS)
					ret = ERR_ALREADY_EXISTS;
			}
		}

		if (ret != ERR_SUCCESS)
			utils_free(rs);

		if (ret == ERR_ALREADY_EXISTS)
			ret = ERR_SUCCESS;
	}

	return ret;
}

static ERR_VALUE _create_long_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, const size_t NumberOfVertices, const READ_PART *ReadPart, const size_t ReadIndex, const boolean RegionStartTouched, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	PKMER_EDGE e = NULL;
	size_t readGapStart = (size_t)-1;
	size_t readGapSeqStart = (size_t)-1;
	PKMER_VERTEX gapStartV = NULL;
	PKMER_VERTEX gapEndV = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	const size_t posPlus = (RegionStartTouched) ? 0 : kmerSize;
	ERR_VALUE ret = ERR_INTERNAL_ERROR; 
	size_t seqIndex = 0;

	for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
		PKMER_VERTEX v = Vertices[i];
		PKMER_VERTEX w = Vertices[i + 1];

		if ((v->Type == kmvtDummy || v->Type == kmvtRead) && w->Type == kmvtRefSeqMiddle) {
			readGapStart = i;
			readGapSeqStart = seqIndex;
			gapStartV = v;
			ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex + posPlus);
		} else if (
			(v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart) && 
			(w->Type == kmvtRead || w->Type == kmvtDummy)) {
			ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex + posPlus);
			if (ret == ERR_SUCCESS) {
				gapEndV = w;
				if (readGapStart != (size_t)-1)
					ret = _create_long_edge(Graph, gapStartV, gapEndV, readGapSeqStart, seqIndex + 1, posPlus, ReadPart, ReadIndex);
			}

			readGapStart = (size_t)-1;
			readGapSeqStart = (size_t)-1;
		} else if (
			(v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart) &&
			w->Type == kmvtRefSeqMiddle) {
			ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex + posPlus);
		} else if ((v->Type == kmvtRead || v->Type == kmvtDummy) && (w->Type == kmvtRead || w->Type == kmvtDummy)) {
			ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex + posPlus);
		}

		if (w->Type != kmvtDummy)
			++seqIndex;
	}

	if (readGapStart != (size_t)-1) {
		for (size_t i = readGapStart; i < NumberOfVertices - 1; ++i) {
			ret = _create_short_edge(Graph, Vertices[i], Vertices[i + 1], ReadIndex, ReadPart->Offset + readGapSeqStart + posPlus);
			if (ret != ERR_SUCCESS)
				break;

			if (Vertices[i + 1]->Type != kmvtDummy)
				++readGapSeqStart;
		}
	}

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read_v2(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex, const uint64_t RegionStart, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	if (Read->ReadSequenceLen > kmerSize) {
		READ_PART *part = NULL;

		part = &Read->Part;
		if (part->ReadSequenceLength > kmerSize){

			const boolean regionStartTouch = (part->Position - RegionStart == 0);
			const size_t maxNumberOfVertices = part->ReadSequenceLength - (kmerSize - 1) + (regionStartTouch ? kmerSize : 0);
			const size_t maxNumberOfEdges = maxNumberOfVertices - 1;
			POINTER_ARRAY_KMER_VERTEX *vertices = NULL;

			ret = utils_calloc(maxNumberOfVertices, sizeof(POINTER_ARRAY_KMER_VERTEX), &vertices);
			if (ret == ERR_SUCCESS) {
				boolean linear = FALSE;

				for (size_t i = 0; i < maxNumberOfVertices; ++i)
					pointer_array_init_KMER_VERTEX(vertices + i, 140);

				ret = _assign_vertice_sets_to_kmers(Graph, part, vertices, maxNumberOfVertices, regionStartTouch, &linear);
				if (ret == ERR_SUCCESS) {
					size_t pathLength = 0;
					PKMER_VERTEX *path = NULL;

					ret = _find_best_path(Graph, vertices, maxNumberOfVertices, linear, &path, &pathLength);
					if (ret == ERR_SUCCESS) {
						ret = _create_long_edges(Graph, path, pathLength, part, ReadIndex, regionStartTouch, PairArray);
						if (ret == ERR_SUCCESS)
							ret = _capture_connect_candidates(Graph, path, pathLength - 1, PairArray);
					}

					utils_free(path);
				}

				for (size_t i = 0; i < maxNumberOfVertices; ++i)
					pointer_array_finit_KMER_VERTEX(vertices + i);

				utils_free(vertices);
			}
		}
	}

	return ret;
}


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const uint32_t Threshold)
{
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_VERTEX sourceVertex = NULL;
	PKMER_VERTEX destVertex = NULL;
	char *beginningSeq = alloca(kmerSize*sizeof(char));

	memset(beginningSeq, 'B', kmerSize*sizeof(char));
	KMER_STACK_ALLOC(sourceKMer, 0, kmerSize, beginningSeq);
	ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRefSeqStart, &sourceVertex);
	if (ret == ERR_SUCCESS) {
		kmer_graph_set_starting_vertex(Graph, sourceKMer);
		KMER_STACK_ALLOC(destKMer, 0, kmerSize, beginningSeq);
		for (size_t i = 0; i < RefSeqLen; ++i) {
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
				ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, 0, 1, kmetReference, &edge);
				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

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
				ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, 0, 1, kmetReference, &edge);
			}
		}
	}

	if (ret == ERR_SUCCESS)
		kmer_graph_set_ending_vertex(Graph, destKMer);

	return ret;
}


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount, const uint64_t RegionStart, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	for (size_t j = 0; j < ReadCount; ++j) {
		ret = _kmer_graph_parse_read_v2(Graph, Reads, j, RegionStart, PairArray);
		if (ret != ERR_SUCCESS)
			break;

		++Reads;
	}

	return ret;
}
