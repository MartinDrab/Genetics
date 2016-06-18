
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

						newDistance += ((kmer_graph_get_edge(Graph, u->KMer, v->KMer) == NULL) ? 3 : 0);
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


static ERR_VALUE _find_probable_paths(PKMER_GRAPH Graph, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfEdges, boolean *Changed)
{
	boolean changed = FALSE;
	ERR_VALUE ret = ERR_SUCCESS;

	do {
		changed = FALSE;
		for (size_t i = 0; i < NumberOfEdges; ++i) {
			PPOINTER_ARRAY_KMER_VERTEX sources = Vertices + i;
			PPOINTER_ARRAY_KMER_VERTEX dests = Vertices + i + 1;

			if (pointer_array_size(sources) == 0 || pointer_array_size(dests) == 0)
				continue;

			if ((*pointer_array_item_KMER_VERTEX(sources, 0))->Type == kmvtRefSeqMiddle &&
				(*pointer_array_item_KMER_VERTEX(dests, 0))->Type == kmvtRefSeqMiddle) {
				if (!Changed[i]) {
					ret = _delete_all_but_ref_edges(Graph, sources, dests, Changed + i);
					changed |= Changed[i];
				}
			}

			if (ret != ERR_SUCCESS)
				break;
		}
	} while (ret == ERR_SUCCESS && changed);

	return ret;
}


static ERR_VALUE _assign_vertice_sets_to_kmers(PKMER_GRAPH Graph, const READ_PART *Part, PPOINTER_ARRAY_KMER_VERTEX Vertices, const size_t NumberOfSets, const uint32_t KMerSize)
{
	PKMER kmer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	KMER_STACK_ALLOC(kmer, 0, KMerSize, Part->ReadSequence);
	for (size_t i = 0; i < NumberOfSets; ++i) {
		ret = kmer_graph_get_vertices(Graph, kmer, Vertices + i);
		if (ret == ERR_SUCCESS) {
			if (pointer_array_size(Vertices + i) == 0) {
				PKMER_VERTEX v = NULL;

				ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
				if (ret == ERR_SUCCESS)
					pointer_array_push_back_no_alloc_KMER_VERTEX(Vertices + i, v);
			}
		}

		if (ret != ERR_SUCCESS)
			break;

		kmer_advance(kmer, Part->ReadSequence[i + KMerSize]);
	}

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read_v2(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Read->ReadSequenceLen <= kmerSize)
		return ERR_SUCCESS;

			GEN_ARRAY_READ_PART readParts;

			dym_array_init_READ_PART(&readParts, 140);
			ret = read_split(Read, &readParts);
			if (ret == ERR_SUCCESS) {
				READ_PART *part = dym_array_item_READ_PART(&readParts, 0);
				
				for (size_t l = 0; l < gen_array_size(&readParts); ++l) {
					if (part->Offset < Read->StartStrips) {
						if (Read->StartStrips - part->Offset >= part->ReadSequenceLength)
							part->ReadSequenceLength = 0;
						else {
							part->ReadSequenceLength -= (Read->StartStrips - part->Offset);
							part->ReadSequence += (Read->StartStrips - part->Offset);
						}
					}
					
//					if (Read->EndStrips > (Read->ReadSequenceLen - part->Offset - part->ReadSequenceLength))
//						part->ReadSequenceLength -= Read->EndStrips;

					if (part->ReadSequenceLength <= kmerSize)
						continue;
					
					boolean rs[2000];
					const size_t maxNumberOfVertices = part->ReadSequenceLength - kmerSize + 1;
					const size_t maxNumberOfEdges = maxNumberOfVertices - 1;
					POINTER_ARRAY_KMER_VERTEX *vertices = NULL;

					memset(rs, FALSE, sizeof(rs));
					ret = utils_calloc(maxNumberOfVertices, sizeof(POINTER_ARRAY_KMER_VERTEX), &vertices);
					if (ret == ERR_SUCCESS) {
						for (size_t i = 0; i < maxNumberOfVertices; ++i)
							pointer_array_init_KMER_VERTEX(vertices + i, 140);

						ret = _assign_vertice_sets_to_kmers(Graph, part, vertices, maxNumberOfVertices, kmerSize);
						if (ret == ERR_SUCCESS) {
							boolean found[500];

							memset(found, 0, sizeof(found));
//							ret = _find_probable_paths(Graph, vertices, maxNumberOfEdges, found);
							ret = _find_best_path(Graph, vertices, maxNumberOfVertices);
							if (ret == ERR_SUCCESS) {
								PPOINTER_ARRAY_KMER_VERTEX sources = vertices;
								PPOINTER_ARRAY_KMER_VERTEX dests = vertices + 1;

								for (size_t i = 0; i < maxNumberOfEdges; ++i) {
									const size_t sourceCount = pointer_array_size(sources);
									const size_t destCount = pointer_array_size(dests);

									if (sourceCount > 0 && destCount > 0) {
										if (!found[i]) {
											for (size_t j = 0; j < sourceCount; ++j) {
												PKMER_VERTEX source = *pointer_array_item_KMER_VERTEX(sources, j);

												for (size_t k = 0; k < destCount; ++k) {
													PKMER_EDGE e = NULL;
													PKMER_VERTEX dest = *pointer_array_item_KMER_VERTEX(dests, k);

													ret = kmer_graph_add_edge_ex(Graph, source, dest, 1, 1, kmetRead, &e);
													if (ret == ERR_ALREADY_EXISTS) {
														e->Weight++;
														ret = ERR_SUCCESS;
													}

													if (ret == ERR_SUCCESS)
														ret = read_info_add(&e->ReadInfo, ReadIndex, i + kmerSize);

													if (ret != ERR_SUCCESS)
														break;
												}

												if (ret != ERR_SUCCESS)
													break;
											}
										} else {
											for (size_t j = 0; j < sourceCount; ++j) {
												PKMER_VERTEX source = *pointer_array_item_KMER_VERTEX(sources, j);

												for (size_t k = 0; k < destCount; ++k) {
													PKMER_EDGE e = NULL;
													PKMER_VERTEX dest = *pointer_array_item_KMER_VERTEX(dests, k);

													e = kmer_graph_get_edge(Graph, source->KMer, dest->KMer);
													if (e != NULL) {
														e->Weight++;
														ret = read_info_add(&e->ReadInfo, ReadIndex, i + kmerSize);
													}

													if (ret != ERR_SUCCESS)
														break;
												}

												if (ret != ERR_SUCCESS)
													break;
											}
										}
									}

									if (ret != ERR_SUCCESS)
										break;

									sources = dests;
									++dests;
								}

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

			dym_array_finit_READ_PART(&readParts);

	return ret;
}


/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const uint32_t Threshold)
{
	PKMER_EDGE lastEdge = NULL;
	EKMerVertexType vertexType = kmvtRefSeqMiddle;
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	KMER_STACK_ALLOC(sourceKMer, 0, kmerSize, RefSeq);
	kmer_back(sourceKMer, 'B');
	ret = kmer_graph_add_vertex(Graph, sourceKMer, kmvtRefSeqStart);
	if (ret == ERR_SUCCESS) {
		kmer_graph_set_starting_vertex(Graph, sourceKMer);
		KMER_STACK_ALLOC(destKMer, 0, kmerSize, RefSeq);
		kmer_back(destKMer, 'B');
		size_t i = kmerSize - 1;
		while (i <= RefSeqLen) {
			if (i == RefSeqLen) {
				kmer_advance(destKMer, 'E');
				vertexType = kmvtRefSeqEnd;
			} else kmer_advance(destKMer, RefSeq[i]);

			kmer_set_number(destKMer, 0);
			ret = kmer_graph_add_vertex(Graph, destKMer, vertexType);
			if (ret == ERR_ALREADY_EXISTS) {
				do {
					kmer_set_number(destKMer, kmer_get_number(destKMer) + 1);
					ret = kmer_graph_add_vertex(Graph, destKMer, vertexType);
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

				lastEdge = edge;
				kmer_advance(sourceKMer, RefSeq[i]);
				kmer_set_number(sourceKMer, kmer_get_number(destKMer));
			}

			if (ret != ERR_SUCCESS)
				break;

			++i;
		}
	}

	if (ret == ERR_SUCCESS)
		kmer_graph_set_ending_vertex(Graph, destKMer);

	return ret;
}


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount, PGEN_ARRAY_KMER_EDGE_PAIR PairArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	for (size_t j = 0; j < ReadCount; ++j) {
		ret = _kmer_graph_parse_read_v2(Graph, Reads, j, PairArray);
		if (ret != ERR_SUCCESS)
			break;

		++Reads;
	}

	return ret;
}


ERR_VALUE read_split(const ONE_READ *Read, PGEN_ARRAY_READ_PART PartArray)
{
	READ_PART part;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (Read->CIGAR != NULL && *Read->CIGAR != '\0' && *Read->CIGAR != '*') {
		char t;
		unsigned long count = 1;
		const char *c = Read->CIGAR;

		ret = ERR_SUCCESS;
		part.ReadSequence = Read->ReadSequence;
		part.Position = Read->Pos;
		part.Offset = 0;
		part.ReadSequenceLength = 0;
		while (ret == ERR_SUCCESS && *c != '\0') {
			char *tmp;

			count = strtoul(c, &tmp, 10);
			c = tmp;
			t = *c;
			if (t != '\0' && count > 0) {
				++c;
				switch (t) {
				case 'M':
				case 'I':
					part.ReadSequenceLength += count;
					break;
				case 'D':
					break;
				case 'S':
					if (part.ReadSequenceLength > 0) {
						part.Offset = (size_t)(part.Position - Read->Pos);
						ret = dym_array_push_back_READ_PART(PartArray, part);
					}

					part.ReadSequence += count;
					part.Position += count;
					part.ReadSequenceLength = 0;
					break;
				case 'H':
					if (part.ReadSequenceLength > 0) {
						part.Offset = (size_t)(part.Position - Read->Pos);
						ret = dym_array_push_back_READ_PART(PartArray, part);
					}

					part.Position += count;
					part.ReadSequenceLength = 0;
					break;
				default:
					part.Offset = 0;
					part.Position = Read->Pos;
					part.ReadSequence = Read->ReadSequence;
					part.ReadSequenceLength = Read->ReadSequenceLen;
					dym_array_clear_READ_PART(PartArray);
					ret = dym_array_push_back_READ_PART(PartArray, part);
					if (ret == ERR_SUCCESS)
						ret = ERR_NO_MORE_ENTRIES;

					break;
				}
			}
		}

		if (ret == ERR_SUCCESS) {
			if (part.ReadSequenceLength > 0) {
				part.Offset = (size_t)(part.Position - Read->Pos);
				ret = dym_array_push_back_READ_PART(PartArray, part);
			}
		}

		if (ret == ERR_NO_MORE_ENTRIES)
			ret = ERR_SUCCESS;
	} else {
		part.Position = Read->Pos;
		part.ReadSequence = Read->ReadSequence;
		part.ReadSequenceLength = Read->ReadSequenceLen;
		part.Offset = 0;
		ret = dym_array_push_back_READ_PART(PartArray, part);
	}

	return ret;
}
