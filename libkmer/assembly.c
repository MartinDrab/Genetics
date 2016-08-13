
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
} DISTANCE_RECORD, *PDISTANCE_RECORD;

GEN_ARRAY_TYPEDEF(DISTANCE_RECORD);
GEN_ARRAY_IMPLEMENTATION(DISTANCE_RECORD)

static ERR_VALUE _find_best_path(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfVertices, boolean Linear, boolean CreateDummyVertices, PKMER_VERTEX **Result, size_t *ResultLength)
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
			POINTER_ARRAY_KMER_VERTEX *nextV = Vertices + 1;
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
						const KMER_EDGE *existingEdge = kmer_graph_get_edge(Graph, &u->KMer, &v->KMer);

						newDistance += ((existingEdge == NULL || (existingEdge->Type == kmetRead  && read_info_get_count(&existingEdge->ReadInfo) == 0)) ? 1 : 0);
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

		ret = utils_calloc(NumberOfVertices*2, sizeof(PKMER_VERTEX), &tmpResult);
		if (ret == ERR_SUCCESS) {
			size_t index = 0;
			PKMER_VERTEX *v = NULL;

			v = pathArray.Data + pointer_array_size(&pathArray);
			for (size_t i = 0; i < pointer_array_size(&pathArray); ++i) {
				PKMER_VERTEX src = *v;
				PKMER_VERTEX dest = *(v - 1);
				
				--v;
				if (CreateDummyVertices && i > 0 &&
					(src->Type != kmvtRead || dest->Type != kmvtRead) &&
					kmer_graph_get_edge(Graph, &src->KMer, &dest->KMer) == NULL) {
					ret = kmer_graph_add_helper_vertex(Graph, &src->KMer, &dest->KMer, &(tmpResult[index]));
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
			*ResultLength = index;
		}
	}

	pointer_array_finit_KMER_VERTEX(&pathArray);

	return ret;
}


static ERR_VALUE _assign_vertice_sets_to_kmers(PKMER_GRAPH Graph, const char *ReadSequence, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfSets, const boolean RegionStartTouch, boolean *Linear)
{
	size_t count = 0;
	PKMER kmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	uint32_t KMerSize = kmer_graph_get_kmer_size(Graph);
	char *beginning = alloca(KMerSize*sizeof(char));
	size_t plusIndex = (RegionStartTouch) ? 0 : KMerSize;

	memset(beginning, 'B', KMerSize*sizeof(char));
	if (!RegionStartTouch)
		memcpy(beginning, ReadSequence, KMerSize*sizeof(char));

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

		kmer_advance(kmer, ReadSequence[i + plusIndex]);
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

	if (ret == ERR_SUCCESS) {
		e->Seq1Weight++;
		ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPosition, ReadQuality);
	}

	*Edge = e;

	return ret;
}


static ERR_VALUE _create_long_edge(PKMER_GRAPH Graph, PKMER_VERTEX U, PKMER_VERTEX V, const size_t StartIndex, const size_t EndIndex, const size_t ExtraOffset, const READ_PART *ReadPart, const size_t ReadIndex, PKMER_EDGE *NewEdge)
{
	PKMER_EDGE e = NULL;
	size_t rsLen = EndIndex - StartIndex - 1;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
		
	ret = kmer_graph_add_edge_ex(Graph, U, V, kmetRead, &e);
	if (ret == ERR_SUCCESS) {
		char *rs = NULL;

		ret = utils_calloc(rsLen + 1, sizeof(char), &rs);
		if (ret == ERR_SUCCESS) {
			memcpy(rs, ReadPart->ReadSequence + StartIndex + ExtraOffset, rsLen*sizeof(char));
			rs[rsLen] = '\0';
			kmer_edge_add_seq(e, kmetRead, rs, rsLen);
		}
	} else if (ret == ERR_ALREADY_EXISTS) {
		if (rsLen == e->SeqLen && memcmp(ReadPart->ReadSequence + StartIndex + ExtraOffset, e->Seq, rsLen*sizeof(char)) == 0) {
			e->Seq1Weight++;
			ret = ERR_SUCCESS;
		} else {
			printf("ONE EDGE - TWO SEQS\n");
			exit(0);
			ret = ERR_SUCCESS;
		}
	}

	if (ret == ERR_SUCCESS) {
		ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + StartIndex + ExtraOffset, 0);
		*NewEdge = e;
	}

	return ret;
}

static ERR_VALUE _create_read_edges(PKMER_GRAPH Graph, PKMER_VERTEX *Vertices, const size_t NumberOfVertices, const READ_PART *ReadPart, const size_t ReadIndex, const size_t StartOffset, PGEN_ARRAY_KMER_EDGE_PAIR PairArray, const boolean LongEdges)
{
	PKMER_EDGE e = NULL;
	size_t readGapStart = (size_t)-1;
	size_t readGapSeqStart = (size_t)-1;
	PKMER_EDGE gapStartE = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR; 
	size_t seqIndex = 0;

	ret = ERR_SUCCESS;
	if (NumberOfVertices > 0) {
		for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
			PKMER_VERTEX v = Vertices[i];
			PKMER_VERTEX w = Vertices[i + 1];

			if (!LongEdges)
				ret = _create_short_edge(Graph, v, w, ReadIndex, ReadPart->Offset + seqIndex + StartOffset, ReadPart->Quality[seqIndex + StartOffset], &e);
			else e = kmer_graph_get_edge(Graph, &v->KMer, &w->KMer);

			if (v->Type == kmvtRead && w->Type == kmvtRefSeqMiddle) {
				readGapStart = i;
				readGapSeqStart = seqIndex;
				gapStartE = e;
			} else if (
				(v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart) &&
				w->Type == kmvtRead) {
				if (LongEdges && readGapStart != (size_t)-1) {
					PKMER_EDGE connectingEdge = NULL;

					ret = _create_long_edge(Graph, gapStartE->Source, w, readGapSeqStart, seqIndex + 1, StartOffset, ReadPart, ReadIndex, &connectingEdge);
					if (ret == ERR_SUCCESS) {
						KMER_EDGE_PAIR p;

						p.U = gapStartE;
						p.V = e;
						p.ConnectingEdge = connectingEdge;
						p.ReadDistance = seqIndex - readGapSeqStart - 1;
						if (!dym_array_contains_KMER_EDGE_PAIR(PairArray, p))
							ret = dym_array_push_back_KMER_EDGE_PAIR(PairArray, p);
					}

				}

				readGapStart = (size_t)-1;
				readGapSeqStart = (size_t)-1;
			}

			if (ret != ERR_SUCCESS)
				break;

			if (!w->Helper)
				++seqIndex;
		}
	}

	return ret;
}


static ERR_VALUE _produce_single_path(PKMER_GRAPH Graph, const char *ReadSequence, const size_t MaxNumberOfSets, const boolean RegionStartTouch, const boolean CreateDummyVertices, PKMER_VERTEX **Path, size_t *PathLength)
{
	boolean linear = FALSE;
	PPOINTER_ARRAY_KMER_VERTEX vertices = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(MaxNumberOfSets, sizeof(POINTER_ARRAY_KMER_VERTEX), &vertices);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < MaxNumberOfSets; ++i)
			pointer_array_init_KMER_VERTEX(vertices + i, 140);

		ret = _assign_vertice_sets_to_kmers(Graph, ReadSequence, vertices, MaxNumberOfSets, RegionStartTouch, &linear);
		if (ret == ERR_SUCCESS)
			ret = _find_best_path(Graph, vertices, MaxNumberOfSets, linear, CreateDummyVertices, Path, PathLength);
	
		for (size_t i = 0; i < MaxNumberOfSets; ++i)
			pointer_array_finit_KMER_VERTEX(vertices + i);

		utils_free(vertices);
	}

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read_v2(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex, const uint64_t RegionStart, PKMER_VERTEX **Path, size_t *PathLength)
{
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	*Path = NULL;
	*PathLength = 0;
	ret = ERR_SUCCESS;
	if (Read->ReadSequenceLen > kmerSize) {
		const READ_PART *part = NULL;

		part = &Read->Part;
		if (part->ReadSequenceLength > kmerSize) {
			const boolean regionStartTouch = (part->Position - RegionStart == 0);
			const size_t maxNumberOfVertices = part->ReadSequenceLength - (kmerSize - 1) + (regionStartTouch ? kmerSize : 0);
			size_t tmpPathLength = 0;

			PKMER_VERTEX *tmpPath = NULL;
			ret = _produce_single_path(Graph, part->ReadSequence, maxNumberOfVertices, regionStartTouch, TRUE, &tmpPath, &tmpPathLength);
			if (ret == ERR_SUCCESS) {
				ret = _create_read_edges(Graph, tmpPath, tmpPathLength, part, ReadIndex, (regionStartTouch ? 0 : kmerSize), NULL, FALSE);
				if (ret == ERR_SUCCESS) {
					*Path = tmpPath;
					*PathLength = tmpPathLength;
				}
			}
		}
	}

	return ret;
}


static ERR_VALUE _fix_read(PKMER_GRAPH Graph, PONE_READ Read, PKMER_VERTEX **Vertices, size_t *NumberOfVertices, const size_t ReadIndex, const uint32_t Threshold, const uint64_t RegionStart, boolean *Fixed)
{
	boolean fixed = FALSE;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	const boolean regionStartTouch = (Read->Part.Position - RegionStart == 0);
	const size_t startOffset = (regionStartTouch) ? 1 : kmerSize;
	const size_t readSeqLen = Read->Part.ReadSequenceLength - (kmerSize - 1) + (regionStartTouch ? kmerSize : 0);
	size_t tmpVertexCount = *NumberOfVertices;
	PKMER_VERTEX *tmpVertices = *Vertices;
	ERR_VALUE ret = ERR_SUCCESS;

	if (Read->Part.ReadSequenceLength > kmerSize) {
		size_t baseIndex = startOffset;

		do {
			char *pBase = Read->Part.ReadSequence + startOffset;
			PKMER_VERTEX u = tmpVertices[0];
			PKMER_VERTEX v = tmpVertices[1];

			fixed = FALSE;
			for (size_t i = 0; i < tmpVertexCount - 1; ++i) {
				PKMER_EDGE e = kmer_graph_get_edge(Graph, &u->KMer, &v->KMer);

				if (pBase - Read->Part.ReadSequence > baseIndex && read_info_weight(&e->ReadInfo, ReadIndex, pBase - Read->Part.ReadSequence) <= (double)Threshold) {
					double maxValue = 0;
					size_t maxIndex = 0;
					size_t rsIndex = (size_t)-1;

					for (size_t j = 0; j < kmer_vertex_out_degree(u); ++j) {
						const KMER_EDGE *f = kmer_vertex_get_succ_edge(u, j);

						if (f == e)
							continue;

						const double weight = read_info_weight(&f->ReadInfo, -1, -1);

						if (f->Type == kmetReference)
							rsIndex = j;

						if (maxValue <= weight) {
							maxValue = weight;
							maxIndex = j;
						}
					}

					if (maxValue > (double)Threshold || rsIndex != (size_t)-1) {
						PKMER_EDGE newPath = kmer_vertex_get_succ_edge(u, maxIndex);

						if (maxValue <= (double)Threshold)
							maxIndex = rsIndex;

						if (newPath->Dest->Helper)
							newPath = kmer_vertex_get_succ_edge(newPath->Dest, 0);

						assert(!newPath->Dest->Helper);
						char newBase = kmer_get_base(&(newPath->Dest->KMer), kmerSize - 1);
						if (*pBase != newBase) {
							baseIndex = pBase - Read->Part.ReadSequence + 1;
							printf("(%Iu:%Iu): %c -> %c\n", ReadIndex, baseIndex - 1, *pBase, newBase);
							Read->Part.Quality[pBase - Read->Part.ReadSequence] = 255;
							*pBase = newBase;
							fixed = TRUE;
							break;
						}
					}
				}

				if (!v->Helper)
					++pBase;

				u = v;
				v = tmpVertices[i + 2];
			}

			if (fixed) {
				*Fixed = TRUE;
				u = tmpVertices[0];
				v = tmpVertices[1];
				for (size_t i = 0; i < tmpVertexCount - 1; ++i) {
					PKMER_EDGE e = kmer_graph_get_edge(Graph, &u->KMer, &v->KMer);

					read_info_remove_2(&e->ReadInfo, ReadIndex);
					u = v;
					v = tmpVertices[i + 2];
				}

				utils_free(tmpVertices);
				ret = _produce_single_path(Graph, Read->Part.ReadSequence, readSeqLen, regionStartTouch, TRUE, &tmpVertices, &tmpVertexCount);
				if (ret == ERR_SUCCESS)
					ret = _create_read_edges(Graph, tmpVertices, tmpVertexCount, &Read->Part, ReadIndex, regionStartTouch, NULL, FALSE);
			}
		} while (ret == ERR_SUCCESS && fixed);
	}

	if (ret == ERR_SUCCESS) {
		*Vertices = tmpVertices;
		*NumberOfVertices = tmpVertexCount;
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
	ret = kmer_alloc(0, kmerSize, beginningSeq, &sourceKMer);
	if (ret == ERR_SUCCESS) {
		ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRefSeqStart, &sourceVertex);
		if (ret == ERR_SUCCESS) {
			kmer_graph_set_starting_vertex(Graph, sourceKMer);
			ret = kmer_alloc(0, kmerSize, beginningSeq, &destKMer);
			if (ret == ERR_SUCCESS) {
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
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
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
						ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, kmetReference, &edge);
					}

					if (ret == ERR_SUCCESS)
						kmer_graph_set_ending_vertex(Graph, destKMer);
				}

				kmer_free(destKMer);
			}
		}
	
		kmer_free(sourceKMer);
	}

	return ret;
}


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, PONE_READ Reads, const size_t ReadCount, const uint64_t RegionStart, const size_t Threshold, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PONE_READ currentRead = NULL;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	PKMER_VERTEX **paths = NULL;
	size_t *pathLengths = 0;


	ret = utils_calloc(ReadCount, sizeof(PKMER_VERTEX *), &paths);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc(ReadCount, sizeof(size_t), &pathLengths);
		if (ret == ERR_SUCCESS) {
			currentRead = Reads;
			for (size_t i = 0; i < ReadCount; ++i) {
				ret = _kmer_graph_parse_read_v2(Graph, currentRead, i, RegionStart, paths + i, pathLengths + i);
				if (ret != ERR_SUCCESS)
					break;

				++currentRead;
			}
		
			if (ret == ERR_SUCCESS) {
				boolean fixed = FALSE;

				do {
					fixed = FALSE;
					currentRead = Reads;
					for (size_t i = 0; i < ReadCount; ++i) {
						ret = _fix_read(Graph, currentRead, paths + i, pathLengths + i, i, Threshold, RegionStart, &fixed);
						if (ret != ERR_SUCCESS)
							break;

						++currentRead;
					}
				} while (ret == ERR_SUCCESS && fixed);

				if (ret == ERR_SUCCESS) {
					currentRead = Reads;
					for (size_t i = 0; i < ReadCount; ++i) {
						ret = _create_read_edges(Graph, paths[i], pathLengths[i], &currentRead->Part, i, ((currentRead->Part.Position == RegionStart) ? 0 : kmerSize), PairArray, TRUE);
						if (ret != ERR_SUCCESS)
							break;

						++currentRead;
					}
				}
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
