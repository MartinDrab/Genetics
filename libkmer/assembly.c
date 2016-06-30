
#include <malloc.h>
#include <stdlib.h>
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


static boolean _vertex_array_is_homogenous(const POINTER_ARRAY_KMER_VERTEX *Array, const EKMerVertexType Type)
{
	boolean ret = TRUE;
	const size_t count = gen_array_size(Array);

	for (size_t i = 0; i < count; ++i) {
		ret = ((*pointer_array_const_item_KMER_VERTEX(Array, i))->Type == Type);
		if (!ret)
			break;
	}

	return ret;
}


static ERR_VALUE _delete_all_but_ref_edges(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_VERTEX Sources, PPOINTER_ARRAY_KMER_VERTEX Dests, boolean *found)
{
	boolean tmpFound = FALSE;
	const size_t sourceCount = pointer_array_size(Sources);
	const size_t destCount = pointer_array_size(Dests);
	POINTER_ARRAY_KMER_VERTEX survivingSources;
	POINTER_ARRAY_KMER_VERTEX survivingDests;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	pointer_array_init_KMER_VERTEX(&survivingSources, 140);
	pointer_array_init_KMER_VERTEX(&survivingDests, 140);
	for (size_t i = 0; i < sourceCount; ++i) {
		PKMER_VERTEX source = *pointer_array_item_KMER_VERTEX(Sources, i);

		for (size_t j = 0; j < destCount; ++j) {
			PKMER_EDGE e = NULL;
			PKMER_VERTEX dest = *pointer_array_item_KMER_VERTEX(Dests, j);

			e = kmer_graph_get_edge(Graph, source->KMer, dest->KMer);
			if (e != NULL) {
				if (e->Type == kmetReference) {
					tmpFound = TRUE;
					ret = pointer_array_push_back_KMER_VERTEX(&survivingSources, source);
					if (ret == ERR_SUCCESS)
						ret = pointer_array_push_back_KMER_VERTEX(&survivingDests, dest);

					break;
				}
			}
		}

		if (ret != ERR_SUCCESS)
			break;
	}

	if (ret == ERR_SUCCESS && tmpFound) {
		pointer_array_clear_KMER_VERTEX(Sources);
		pointer_array_clear_KMER_VERTEX(Dests);
		ret = pointer_array_push_back_array_KMER_VERTEX(Sources, &survivingSources);
		if (ret == ERR_SUCCESS) {
			ret = pointer_array_push_back_array_KMER_VERTEX(Dests, &survivingDests);
		}
	}

	pointer_array_finit_KMER_VERTEX(&survivingDests);
	pointer_array_finit_KMER_VERTEX(&survivingSources);
	*found = tmpFound;

	return ret;
}


static ERR_VALUE _get_edges(const KMER_GRAPH *Graph, const POINTER_ARRAY_KMER_VERTEX *Sources, const POINTER_ARRAY_KMER_VERTEX *Dests, PPOINTER_ARRAY_KMER_EDGE Edges)
{
	const KMER_EDGE *e = NULL;
	const KMER_VERTEX *s = NULL;
	const KMER_VERTEX *d = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < pointer_array_size(Sources); ++i) {
		s = *pointer_array_const_item_KMER_VERTEX(Sources, i);
		for (size_t j = 0; j < pointer_array_size(Dests); ++j) {
			d = *pointer_array_const_item_KMER_VERTEX(Dests, j);
			e = kmer_graph_get_edge(Graph, s->KMer, d->KMer);
			if (e != NULL) {
				ret = pointer_array_push_back_KMER_EDGE(Edges, e);
				if (ret != ERR_SUCCESS)
					break;
			}
		}
	}
	return ret;
}


static ERR_VALUE _capture_connect_candidates(const KMER_GRAPH *Graph, const POINTER_ARRAY_KMER_VERTEX *Vertices, const size_t NumberOfEdges, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t firstReadPoint = (size_t)-1;
	size_t lastReadPoint = (size_t)-1;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < NumberOfEdges; ++i) {
		PKMER_VERTEX u = *pointer_array_item_KMER_VERTEX(Vertices + i, 0);
		PKMER_VERTEX v = *pointer_array_item_KMER_VERTEX(Vertices + i + 1, 0);
		PKMER_EDGE e = kmer_graph_get_edge(Graph, u->KMer, v->KMer);

		if (e != NULL && e->Type == kmetRead) {
			lastReadPoint = i;
			if (firstReadPoint != (size_t)-1 && e->Source->Type == kmvtRefSeqMiddle) {
				const POINTER_ARRAY_KMER_VERTEX *sourceArray = Vertices + firstReadPoint + 1;
				const POINTER_ARRAY_KMER_VERTEX *destArray = Vertices + lastReadPoint;

				for (size_t j = 0; j < pointer_array_size(sourceArray); ++j) {
					const KMER_VERTEX *u = *pointer_array_const_item_KMER_VERTEX(sourceArray, j);

					assert(u->Type == kmvtRefSeqMiddle);
					for (size_t k = 0; k < pointer_array_size(destArray); ++k) {
						const KMER_VERTEX *v = *pointer_array_const_item_KMER_VERTEX(destArray, k);
						POINTER_ARRAY_KMER_EDGE sourceEdges;
						POINTER_ARRAY_KMER_EDGE destEdges;

						assert(v->Type == kmvtRefSeqMiddle);
						pointer_array_init_KMER_EDGE(&sourceEdges, 140);
						pointer_array_init_KMER_EDGE(&destEdges, 140);
						kmer_vertex_get_certain_edges(u, kmetRead, TRUE, &sourceEdges);
						kmer_vertex_get_certain_edges(v, kmetRead, FALSE, &destEdges);
						for (size_t l = 0; l < pointer_array_size(&sourceEdges); ++l) {
							const KMER_EDGE *se = *pointer_array_const_item_KMER_EDGE(&sourceEdges, l);

							for (size_t m = 0; m < pointer_array_size(&destEdges); ++m) {
								const KMER_EDGE *de = *pointer_array_const_item_KMER_EDGE(&destEdges, m);
								KMER_EDGE_PAIR pair;

								if (se->Dest->Order <= de->Source->Order) {
									pair.U = se;
									pair.V = de;
									pair.ReadDistance = lastReadPoint - firstReadPoint - 1;
									assert(pair.U->Type == kmetRead);
									assert(pair.V->Type == kmetRead);
									if (!dym_array_contains_KMER_EDGE_PAIR(PairArray, pair))
										ret = dym_array_push_back_KMER_EDGE_PAIR(PairArray, pair);
								}

								if (ret != ERR_SUCCESS)
									break;
							}

							if (ret != ERR_SUCCESS)
								break;
						}

						pointer_array_finit_KMER_EDGE(&destEdges);
						pointer_array_finit_KMER_EDGE(&sourceEdges);
						if (ret != ERR_SUCCESS)
							break;
					}

					if (ret != ERR_SUCCESS)
						break;
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			if (e->Dest->Type == kmvtRefSeqMiddle)
				firstReadPoint = lastReadPoint;

			lastReadPoint = (size_t)-1;
		}

		if (ret != ERR_SUCCESS)
			break;
	}

	return ret;
}


typedef struct _DISTANCE_RECORD {
	size_t Distance;
	size_t BackIndex;
} DISTANCE_RECORD, *PDISTANCE_RECORD;

GEN_ARRAY_TYPEDEF(DISTANCE_RECORD);
GEN_ARRAY_IMPLEMENTATION(DISTANCE_RECORD)

static ERR_VALUE _find_best_path(const KMER_GRAPH *Graph, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfVertices)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	DISTANCE_RECORD dr;
	GEN_ARRAY_DISTANCE_RECORD *distances = NULL;

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

						newDistance += ((kmer_graph_get_edge(Graph, u->KMer, v->KMer) == NULL) ? 1 : 0);
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

				{
					const KMER_VERTEX *v = *pointer_array_const_item_KMER_VERTEX(currentV, minDistanceIndex);

					pointer_array_clear_KMER_VERTEX(currentV);
					pointer_array_push_back_no_alloc_KMER_VERTEX(currentV, v);
				}

				dr.BackIndex = minDistanceBackIndex;
				dr.Distance = minDistance;
				for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
					--currentV;
					--currentD;
					const KMER_VERTEX *v = *pointer_array_const_item_KMER_VERTEX(currentV, dr.BackIndex);

					dr = *dym_array_item_DISTANCE_RECORD(currentD, dr.BackIndex);
					pointer_array_clear_KMER_VERTEX(currentV);
					pointer_array_push_back_no_alloc_KMER_VERTEX(currentV, v);
				}
			}

			for (size_t i = 0; i < NumberOfVertices; ++i)
				dym_array_finit_DISTANCE_RECORD(distances + i);
		}

		utils_free(distances);
	}

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

	ret = ERR_SUCCESS;
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

	*Linear = (count == NumberOfSets);

	return ret;
}


static ERR_VALUE _create_long_edges(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfVertices, const READ_PART *ReadPart, const size_t ReadIndex, const boolean RegionStartTouched)
{
	PKMER_EDGE e = NULL;
	size_t readGapStart = (size_t)-1;
	PKMER_VERTEX gapStartV = NULL;
	size_t readGapEnd = (size_t)-1;
	PKMER_VERTEX gapEndV = NULL;
	PPOINTER_ARRAY_KMER_VERTEX current = Vertices;
	PPOINTER_ARRAY_KMER_VERTEX next = Vertices + 1;
	const size_t kmerSize = kmer_graph_get_kmer_size(Graph);
	const size_t posPlus = (RegionStartTouched) ? 0 : kmerSize;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	for (size_t i = 0; i < NumberOfVertices - 1; ++i) {
		PKMER_VERTEX v = *pointer_array_item_KMER_VERTEX(current, 0);
		PKMER_VERTEX w = *pointer_array_item_KMER_VERTEX(next, 0);

		if (v->Type == kmvtRead && w->Type == kmvtRefSeqMiddle) {
			readGapStart = i;
			gapStartV = v;
			ret = kmer_graph_add_edge_ex(Graph, v, w, 1, 1, kmetRead, &e);
			if (ret == ERR_ALREADY_EXISTS) {
				e->Weight++;
				ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS)
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + i + posPlus);
		} else if (
			(v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart) && 
			w->Type == kmvtRead) {
			readGapEnd = i + 1;
			gapEndV = w;
			if (readGapStart != (size_t)-1) {
				char *rs = NULL;
				size_t rsLen = readGapEnd - readGapStart - 1;

				ret = utils_calloc(rsLen + 1, sizeof(char), &rs);
				if (ret == ERR_SUCCESS) {
					memcpy(rs, ReadPart->ReadSequence + readGapStart + posPlus, rsLen*sizeof(char));
					rs[rsLen] = '\0';
					ret = kmer_graph_add_edge_ex(Graph, gapStartV, gapEndV, 1, readGapEnd - readGapStart, kmetRead, &e);
					if (ret == ERR_SUCCESS) {
						e->Seq = rs;
						e->SeqLen = rsLen;
						e->SeqType = kmetRead;
						ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + readGapStart + posPlus);
					} else if (ret == ERR_ALREADY_EXISTS) {
						if (rsLen = e->SeqLen && memcmp(e->Seq, rs, rsLen) == 0) {
							e->Weight++;
							ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + readGapStart + posPlus);
							if (ret == ERR_SUCCESS)
								ret = ERR_ALREADY_EXISTS;
						}
					}

					if (ret != ERR_SUCCESS)
						utils_free(rs);

					if (ret == ERR_ALREADY_EXISTS)
						ret = ERR_SUCCESS;
				}
			}

			readGapStart = (size_t)-1;
			readGapEnd = (size_t)-1;
			ret = kmer_graph_add_edge_ex(Graph, v, w, 1, 1, kmetRead, &e);
			if (ret == ERR_ALREADY_EXISTS) {
				e->Weight++;
				ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS)
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + i + posPlus);
		} else if (
			(v->Type == kmvtRefSeqMiddle || v->Type == kmvtRefSeqStart) &&
			w->Type == kmvtRefSeqMiddle) {
			ret = kmer_graph_add_edge_ex(Graph, v, w, 0, 1, kmetRead, &e);
			if (ret == ERR_ALREADY_EXISTS)
				ret = ERR_SUCCESS;
			
			if (ret == ERR_SUCCESS && (e->Type == kmetRead || readGapStart == (size_t)-1)) {
				ret = ERR_SUCCESS;
				e->Weight++;
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + i + posPlus);
			}
		} else if (v->Type == kmvtRead && w->Type == kmvtRead) {
			ret = kmer_graph_add_edge_ex(Graph, v, w, 1, 1, kmetRead, &e);
			if (ret == ERR_ALREADY_EXISTS) {
				e->Weight++;
				ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS)
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + i + posPlus);
		}

		current = next;
		++next;
	}

	if (readGapStart != (size_t)-1) {
		current = Vertices + readGapStart;
		next = current + 1;
		for (size_t i = readGapStart; i < NumberOfVertices - 1; ++i) {
			PKMER_VERTEX v = *pointer_array_item_KMER_VERTEX(current, 0);
			PKMER_VERTEX w = *pointer_array_item_KMER_VERTEX(next, 0);

			ret = kmer_graph_add_edge_ex(Graph, v, w, 1, 1, kmetRead, &e);
			if (ret == ERR_ALREADY_EXISTS) {
				e->Weight++;
				ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS)
				ret = read_info_add(&e->ReadInfo, ReadIndex, ReadPart->Offset + i + posPlus);

			current = next;
			++next;
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

		part = Read->Parts.Data;
		for (size_t l = 0; l < gen_array_size(&Read->Parts); ++l) {
			if (part->ReadSequenceLength <= kmerSize) {
				++part;
				continue;
			}

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
					if (!linear)
						ret = _find_best_path(Graph, vertices, maxNumberOfVertices);
					
					if (ret == ERR_SUCCESS) {
						ret = _create_long_edges(Graph, vertices, maxNumberOfVertices, part, ReadIndex, regionStartTouch);
						if (ret == ERR_SUCCESS)
							ret = _capture_connect_candidates(Graph, vertices, maxNumberOfEdges, PairArray);
						
					}
				}

				for (size_t i = 0; i < maxNumberOfVertices; ++i)
					pointer_array_finit_KMER_VERTEX(vertices + i);

				utils_free(vertices);
			}

			++part;
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
	char *beginningSeq = alloca(kmerSize*sizeof(char));

	memset(beginningSeq, 'B', kmerSize*sizeof(char));
	KMER_STACK_ALLOC(sourceKMer, 0, kmerSize, beginningSeq);
	ret = kmer_graph_add_vertex(Graph, sourceKMer, kmvtRefSeqStart);
	if (ret == ERR_SUCCESS) {
		kmer_graph_set_starting_vertex(Graph, sourceKMer);
		KMER_STACK_ALLOC(destKMer, 0, kmerSize, beginningSeq);
		for (size_t i = 0; i < RefSeqLen; ++i) {
			kmer_advance(destKMer, RefSeq[i]);
			kmer_set_number(destKMer, 0);
			ret = kmer_graph_add_vertex(Graph, destKMer, kmvtRefSeqMiddle);
			if (ret == ERR_ALREADY_EXISTS) {
				do {
					kmer_set_number(destKMer, kmer_get_number(destKMer) + 1);
					ret = kmer_graph_add_vertex(Graph, destKMer, kmvtRefSeqMiddle);
				} while (ret == ERR_ALREADY_EXISTS);
			}

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				PKMER_VERTEX sourceVertex = NULL;
				PKMER_VERTEX destVertex = NULL;

				sourceVertex = kmer_graph_get_vertex(Graph, sourceKMer);
				destVertex = kmer_graph_get_vertex(Graph, destKMer);
				destVertex->RefSeqPosition = i;
				ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, 0, 1, kmetReference, &edge);
				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

				kmer_advance(sourceKMer, RefSeq[i]);
				kmer_set_number(sourceKMer, kmer_get_number(destKMer));
			}

			if (ret != ERR_SUCCESS)
				break;
		}

		if (ret == ERR_SUCCESS) {
			kmer_advance(destKMer, 'E');
			kmer_set_number(destKMer, 0);
			ret = kmer_graph_add_vertex(Graph, destKMer, kmvtRefSeqEnd);
			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				PKMER_VERTEX sourceVertex = NULL;
				PKMER_VERTEX destVertex = NULL;

				sourceVertex = kmer_graph_get_vertex(Graph, sourceKMer);
				destVertex = kmer_graph_get_vertex(Graph, destKMer);
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


