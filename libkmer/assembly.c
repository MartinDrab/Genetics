
#include <malloc.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
#include "reads.h"
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

//	printf("%Iu, %Iu ---> ", sourceCount, destCount);
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
//	printf("%Iu, %Iu\n", pointer_array_size(Sources), pointer_array_size(Dests));

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read_v2(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex)
{
	const char *readSeq = Read->ReadSequence;
	const size_t readLen = Read->ReadSequenceLen;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE edgeTable = NULL;
	PKMER_TABLE vertexTable = NULL;

	ret = kmer_table_create(kmerSize, 2, utils_next_prime(readLen), NULL, &vertexTable);
	if (ret == ERR_SUCCESS) {
		ret = kmer_edge_table_create(kmerSize, 2, utils_next_prime(readLen), NULL, &edgeTable);
		if (ret == ERR_SUCCESS) {
			boolean rs[2000];
			const size_t maxNumberOfVertices = readLen - kmerSize + 1;
			const size_t maxNumberOfEdges = maxNumberOfVertices - 1;
			POINTER_ARRAY_KMER_VERTEX *vertices = NULL;

			memset(rs, FALSE, sizeof(rs));
			ret = utils_calloc(maxNumberOfVertices, sizeof(POINTER_ARRAY_KMER_VERTEX), &vertices);
			if (ret == ERR_SUCCESS) {
				PKMER kmer = NULL;

				for (size_t i = 0; i < maxNumberOfVertices; ++i)
					pointer_array_init_KMER_VERTEX(vertices + i, 140);

				KMER_STACK_ALLOC(kmer, 0, kmerSize, readSeq);
				for (size_t i = 0; i < maxNumberOfVertices; ++i) {
					ret = kmer_graph_get_vertices(Graph, kmer, vertices + i);
					if (ret == ERR_SUCCESS) {
						if (pointer_array_size(vertices + i) == 0) {
							PKMER_VERTEX v = NULL;

							ret = kmer_graph_add_vertex_ex(Graph, kmer, kmvtRead, &v);
							if (ret == ERR_SUCCESS)
								pointer_array_push_back_no_alloc_KMER_VERTEX(vertices + i, v);
						}
					}

					if (ret != ERR_SUCCESS)
						break;
					
					kmer_advance(kmer, readSeq[i + kmerSize]);
				}


				if (ret == ERR_SUCCESS) {
					/*
					{
						size_t firstReadPoint = (size_t)-1;
						size_t nextReadPoint = (size_t)-1;
						for (size_t i = 0; i < maxNumberOfVertices; ++i) {
							PPOINTER_ARRAY_KMER_VERTEX vs = vertices + i;
							const EKMerVertexType vertexType = (*pointer_array_const_item_KMER_VERTEX(vs, 0))->Type;

							if (vertexType == kmvtRead || (pointer_array_size(vs) == 1 && firstReadPoint != (size_t)-1)) {
								nextReadPoint = i;
								if (firstReadPoint != (size_t)-1) {
									if (nextReadPoint - firstReadPoint > 1) {
										PKMER_EDGE e = NULL;
										size_t length = nextReadPoint - firstReadPoint;
										PKMER_VERTEX s = *pointer_array_item_KMER_VERTEX(vertices + firstReadPoint, 0);
										PKMER_VERTEX d = *pointer_array_item_KMER_VERTEX(vertices + nextReadPoint, 0);

										ret = kmer_graph_add_edge_ex(Graph, s, d, 1, length, kmetRead, &e);
										if (ret == ERR_SUCCESS) {
											e->SeqLen = length - 1;
											ret = utils_calloc(e->SeqLen + 1, sizeof(char), &e->Seq);
											if (ret == ERR_SUCCESS) {
												memcpy(e->Seq, readSeq + kmerSize + i, length - 1);
												e->Seq[e->SeqLen] = '\0';
											}
										}

										if (ret == ERR_ALREADY_EXISTS && e->SeqLen == length - 1 && memcmp(e->Seq, readSeq + kmerSize + i, e->SeqLen) == 0) {
											e->Weight++;
											ret = ERR_SUCCESS;
										}

										if (ret == ERR_SUCCESS) {
											for (size_t j = firstReadPoint + 1; j < nextReadPoint; ++j)
												pointer_array_clear_KMER_VERTEX(vertices + j);
										}

										if (ret == ERR_ALREADY_EXISTS)
											ret = ERR_SUCCESS;
									}
								}

								firstReadPoint = (size_t)-1;
								if (vertexType == kmvtRead)
									firstReadPoint = nextReadPoint;
							}

							if (ret != ERR_SUCCESS)
								break;
						}
					}
					*/
					/*
					{
						if (pointer_array_size(vertices) > 1 && pointer_array_size(vertices + 1) == 1)
							pointer_array_clear_KMER_VERTEX(vertices);

						if (pointer_array_size(vertices + maxNumberOfVertices - 1) > 1 && pointer_array_size(vertices + maxNumberOfVertices - 2) == 1)
							pointer_array_clear_KMER_VERTEX(vertices + maxNumberOfVertices - 1);
					}
					*/
					boolean found[500];
					boolean changed = FALSE;

					memset(found, FALSE, sizeof(found));
					do {
						changed = FALSE;
						for (size_t i = 0; i < maxNumberOfEdges; ++i) {
							PPOINTER_ARRAY_KMER_VERTEX sources = vertices + i;
							PPOINTER_ARRAY_KMER_VERTEX dests = vertices + i + 1;

							if (pointer_array_size(sources) == 0 || pointer_array_size(dests) == 0)
								continue;

							if ((*pointer_array_item_KMER_VERTEX(sources, 0))->Type == kmvtRefSeqMiddle &&
								(*pointer_array_item_KMER_VERTEX(dests, 0))->Type == kmvtRefSeqMiddle) {
								if (!found[i]) {
									ret = _delete_all_but_ref_edges(Graph, sources, dests, found + i);
									changed |= found[i];
								}
							}

							if (ret != ERR_SUCCESS)
								break;
						}
					} while (ret == ERR_SUCCESS && changed);

					for (size_t i = 0; i < maxNumberOfEdges; ++i) {
						PPOINTER_ARRAY_KMER_VERTEX sources = vertices + i;
						PPOINTER_ARRAY_KMER_VERTEX dests = vertices + i + 1;
						const size_t sourceCount = pointer_array_size(sources);
						const size_t destCount = pointer_array_size(dests);

						if (sourceCount == 0 || destCount == 0)
							continue;

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

						if (ret != ERR_SUCCESS)
							break;
					}
				}

				for (size_t i = 0; i < maxNumberOfVertices; ++i)
					pointer_array_finit_KMER_VERTEX(vertices + i);

				utils_free(vertices);
			}

			kmer_edge_table_destroy(edgeTable);
		}

		kmer_table_destroy(vertexTable);
	}

	return ret;
}

static ERR_VALUE _kmer_graph_parse_read(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex)
{
	const char *readSeq = Read->ReadSequence;
	const size_t readLen = Read->ReadSequenceLen;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE edgeTable = NULL;
	PKMER_TABLE vertexTable = NULL;

	ret = kmer_table_create(kmerSize, 2, utils_next_prime(readLen), NULL, &vertexTable);
	if (ret == ERR_SUCCESS) {
		ret = kmer_edge_table_create(kmerSize, 2, utils_next_prime(readLen), NULL, &edgeTable);
		if (ret == ERR_SUCCESS) {
			POINTER_ARRAY_KMER_VERTEX sourceVertices;
			PKMER sourceKMer = NULL;
			PKMER lastRefSeqKmer = NULL;
			char seq[8192];
			size_t seqLen = 0;

			pointer_array_init_KMER_VERTEX(&sourceVertices, 140);
			KMER_STACK_ALLOC(lastRefSeqKmer, 0, kmerSize, NULL);
			KMER_STACK_ALLOC(sourceKMer, 0, kmerSize, readSeq);
			ret = kmer_graph_get_vertices(Graph, sourceKMer, &sourceVertices);
			if (ret == ERR_SUCCESS && pointer_array_size(&sourceVertices) == 0) {
				PKMER_VERTEX v = NULL;

				ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRead, &v);
				if (ret == ERR_SUCCESS)
					ret = pointer_array_push_back_KMER_VERTEX(&sourceVertices, v);
			}

			if (ret == ERR_SUCCESS && pointer_array_size(&sourceVertices) == 1) {
				PKMER destKMer = NULL;
				POINTER_ARRAY_KMER_VERTEX destVertices;

				pointer_array_init_KMER_VERTEX(&destVertices, 140);
				KMER_STACK_ALLOC(destKMer, 0, kmerSize, readSeq);
				size_t i = kmerSize;
				while (i < readLen) {
					if ((*pointer_array_item_KMER_VERTEX(&sourceVertices, 0))->Type == kmvtRefSeqMiddle)
						kmer_init_from_kmer(lastRefSeqKmer, sourceKMer);

					kmer_advance(destKMer, readSeq[i]);
					pointer_array_clear_KMER_VERTEX(&destVertices);
					kmer_set_number(destKMer, 0);
					ret = kmer_graph_get_vertices(Graph, destKMer, &destVertices);
					if (ret == ERR_SUCCESS && pointer_array_size(&destVertices) == 0) {
						PKMER_VERTEX v = NULL;

						ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRead, &v);
						if (ret == ERR_SUCCESS)
							ret = pointer_array_push_back_KMER_VERTEX(&destVertices, v);
					}

					if (ret == ERR_SUCCESS) {
						PKMER_EDGE edge = NULL;
						boolean readEdgeRepeat = FALSE;
						boolean readVertexRepeat = FALSE;

						do {
							ret = kmer_edge_table_insert(edgeTable, sourceKMer, destKMer, (void *)TRUE);
							if (ret == ERR_TABLE_FULL) {
								ret = kmer_edge_table_extend(edgeTable);
								if (ret == ERR_SUCCESS)
									ret = ERR_TABLE_FULL;
							}
						} while (ret == ERR_TABLE_FULL);

						if (ret == ERR_ALREADY_EXISTS) {
							readEdgeRepeat = TRUE;
							ret = ERR_SUCCESS;
						}

						do {
							ret = kmer_table_insert(vertexTable, destKMer, (void *)TRUE);
							if (ret == ERR_TABLE_FULL) {
								ret = kmer_table_extend(vertexTable);
								if (ret == ERR_SUCCESS)
									ret = ERR_TABLE_FULL;
							}
						} while (ret == ERR_TABLE_FULL);

						if (ret == ERR_ALREADY_EXISTS) {
							readVertexRepeat = TRUE;
							ret = ERR_SUCCESS;
						}

						if (!_vertex_array_is_homogenous(&sourceVertices, kmvtRefSeqMiddle) ||
							!_vertex_array_is_homogenous(&destVertices, kmvtRefSeqMiddle)) {
								PKMER_VERTEX sourceVertex = NULL;
								PKMER_VERTEX destVertex = NULL;

								destVertex = *pointer_array_item_KMER_VERTEX(&destVertices, 0);
								if (!kmer_seq_equal(lastRefSeqKmer, destKMer) && (i == readLen - 1 || destVertex->Type == kmvtRefSeqMiddle || !readVertexRepeat)) {
									seq[seqLen] = '\0';
									for (size_t j = 0; j < gen_array_size(&sourceVertices); ++j) {
										sourceVertex = *pointer_array_item_KMER_VERTEX(&sourceVertices, j);
										for (size_t k = 0; k < gen_array_size(&destVertices); ++k) {
											destVertex = *pointer_array_item_KMER_VERTEX(&destVertices, k);
											ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, 1, seqLen + 1, kmetRead, &edge);
											if (ret == ERR_SUCCESS) {
												edge->SeqLen = seqLen;
												ret = utils_copy_string(seq, &edge->Seq);
											}

											if (ret == ERR_ALREADY_EXISTS) {
												if (!readEdgeRepeat)
													edge->Weight++;

												ret = ERR_SUCCESS;
											}

											if (ret == ERR_SUCCESS)
												ret = read_info_add(&edge->ReadInfo, ReadIndex, i);

											if (ret != ERR_SUCCESS)
												break;
										}

										if (ret != ERR_SUCCESS)
											break;
									}

									seqLen = 0;
								} else {
									printf("read repetition detected\n");
									seq[seqLen] = readSeq[i];
									++seqLen;
								}
						} else {
							boolean edgeFound = FALSE;

							for (size_t j = 0; j < gen_array_size(&sourceVertices); ++j) {
								PKMER_VERTEX v = *pointer_array_item_KMER_VERTEX(&sourceVertices, j);

								for (size_t k = 0; k < gen_array_size(&destVertices); ++k) {
									PKMER_VERTEX w = *pointer_array_item_KMER_VERTEX(&destVertices, k);
									PKMER_EDGE e = kmer_graph_get_edge(Graph, v->KMer, w->KMer);

									if (e != NULL) {
										edgeFound = (e->Type == kmetReference);
										if (edgeFound) {
											if (!readEdgeRepeat)
												e->Weight++;

											ret = read_info_add(&e->ReadInfo, ReadIndex, i);
											pointer_array_clear_KMER_VERTEX(&destVertices);
											pointer_array_push_back_no_alloc_KMER_VERTEX(&destVertices, w);
											break;
										}
									}
								}
							}

							if (!edgeFound) {
								for (size_t j = 0; j < gen_array_size(&sourceVertices); ++j) {
									PKMER_VERTEX v = *pointer_array_item_KMER_VERTEX(&sourceVertices, j);

									for (size_t k = 0; k < gen_array_size(&destVertices); ++k) {
										PKMER_VERTEX w = *pointer_array_item_KMER_VERTEX(&destVertices, k);

										if (v == w)
											continue;

										ret = kmer_graph_add_edge_ex(Graph, v, w, 1, 1, kmetRead, &edge);
										if (ret == ERR_ALREADY_EXISTS) {
											if (!readEdgeRepeat)
												++edge->Weight;

											ret = ERR_SUCCESS;
										}

										if (ret == ERR_SUCCESS)
											ret = read_info_add(&edge->ReadInfo, ReadIndex, i);

										if (ret != ERR_SUCCESS)
											break;
									}

									if (ret != ERR_SUCCESS)
										break;
								}
							}
						}

						if (ret == ERR_SUCCESS) {
							kmer_advance(sourceKMer, readSeq[i]);
							pointer_array_exchange_KMER_VERTEX(&sourceVertices, &destVertices);
						}
					}

					if (ret != ERR_SUCCESS)
						break;

					++i;
				}

				pointer_array_finit_KMER_VERTEX(&destVertices);
			}

			pointer_array_finit_KMER_VERTEX(&sourceVertices);
			kmer_edge_table_destroy(edgeTable);
		}
	
		kmer_table_destroy(vertexTable);
	}

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
				if (!kmer_equal(sourceVertex->KMer, sourceKMer))
					__debugbreak();

				if (!kmer_equal(destVertex->KMer, destKMer))
					__debugbreak();

				ret = kmer_graph_add_edge_ex(Graph, sourceVertex, destVertex, Threshold, 1, kmetReference, &edge);
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


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER kmer = NULL;
	GEN_ARRAY_ONE_READ readArray;
	size_t currentIndex = 0;

	dym_array_init_ONE_READ(&readArray, 140);
	ret = dym_array_reserve_ONE_READ(&readArray, ReadCount);
	if (ret == ERR_SUCCESS) {
		ret = kmer_alloc(0, kmer_graph_get_kmer_size(Graph), NULL, &kmer);
		if (ret == ERR_SUCCESS) {
			ret = ERR_SUCCESS;
			for (size_t j = 0; j < ReadCount; ++j) {
				KMER_VERTEX *v = NULL;

				kmer_init(kmer, Reads->ReadSequence);
				v = kmer_graph_get_vertex(Graph, kmer);
				if (v != NULL && v->Type == kmvtRefSeqMiddle) {
					ret = _kmer_graph_parse_read_v2(Graph, Reads, currentIndex);
					++currentIndex;
				} else dym_array_push_back_no_alloc_ONE_READ(&readArray, *Reads);

				if (ret != ERR_SUCCESS)
					break;

				++Reads;
			}

			for (size_t i = 0; i < gen_array_size(&readArray); ++i) {
				ret = _kmer_graph_parse_read_v2(Graph, dym_array_item_ONE_READ(&readArray, i), currentIndex);
				++currentIndex;
				if (ret != ERR_SUCCESS)
					break;
			}

			kmer_free(kmer);
		}
	}

	dym_array_finit_ONE_READ(&readArray);

	return ret;
}
