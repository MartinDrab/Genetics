
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


static boolean _vertex_array_is_homogenous(const DYM_ARRAY *Array, const EKMerVertexType Type)
{
	boolean ret = TRUE;

	for (size_t i = 0; i < dym_array_size(Array); ++i) {
		ret = (((PKMER_VERTEX)dym_array_get(Array, i))->Type == Type);
		if (!ret)
			break;
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

	ret = kmer_edge_table_create(kmerSize, 2, 37, NULL, &edgeTable);
	if (ret == ERR_SUCCESS) {
		DYM_ARRAY sourceVertices;
		PKMER sourceKMer = NULL;

		dym_array_create(&sourceVertices, 140);
		KMER_STACK_ALLOC(sourceKMer, 0, kmerSize, readSeq);
		ret = kmer_graph_get_vertices(Graph, sourceKMer, &sourceVertices);
		if (ret == ERR_SUCCESS && dym_array_size(&sourceVertices) == 0) {
			PKMER_VERTEX v = NULL;
			
			ret = kmer_graph_add_vertex_ex(Graph, sourceKMer, kmvtRead, &v);
			if (ret == ERR_SUCCESS)
				ret = dym_array_push_back(&sourceVertices, v);
		}

		if (ret == ERR_SUCCESS) {
			PKMER destKMer = NULL;
			DYM_ARRAY destVertices;

			dym_array_create(&destVertices, 140);
			KMER_STACK_ALLOC(destKMer, 0, kmerSize, readSeq);
			for (size_t i = kmerSize; i < readLen; ++i) {
				kmer_advance(destKMer, readSeq[i]);
				dym_array_clear(&destVertices);
				kmer_set_number(destKMer, 0);
				ret = kmer_graph_get_vertices(Graph, destKMer, &destVertices);
				if (ret == ERR_SUCCESS && dym_array_size(&destVertices) == 0) {
					PKMER_VERTEX v = NULL;

					ret = kmer_graph_add_vertex_ex(Graph, destKMer, kmvtRead, &v);
					if (ret == ERR_SUCCESS)
						ret = dym_array_push_back(&destVertices, v);
				}

				if (ret == ERR_SUCCESS) {
					PKMER_EDGE edge = NULL;
					boolean readRepeat = FALSE;

					do {
						ret = kmer_edge_table_insert(edgeTable, sourceKMer, destKMer, (void *)TRUE);
						if (ret == ERR_TABLE_FULL) {
							ret = kmer_edge_table_extend(edgeTable);
							if (ret == ERR_SUCCESS)
								ret = ERR_TABLE_FULL;
						}
					} while (ret == ERR_TABLE_FULL);
					
					if (ret == ERR_ALREADY_EXISTS) {
						readRepeat = TRUE;
						ret = ERR_SUCCESS;
					}

					if (!_vertex_array_is_homogenous(&sourceVertices, kmvtRefSeqMiddle) ||
						!_vertex_array_is_homogenous(&destVertices, kmvtRefSeqMiddle)) {						
						for (size_t j = 0; j < dym_array_size(&sourceVertices); ++j) {
							PKMER_VERTEX sourceVertex = (PKMER_VERTEX)dym_array_get(&sourceVertices, j);

							for (size_t k = 0; k < dym_array_size(&destVertices); ++k) {
								PKMER_VERTEX destVertex = (PKMER_VERTEX)dym_array_get(&destVertices, k);

								ret = kmer_graph_add_edge_ex(Graph, sourceVertex->KMer, destVertex->KMer, 1, 1, kmetRead, &edge);
								if (ret == ERR_ALREADY_EXISTS) {
									if (!readRepeat)
										edge->Weight++;

									ret = ERR_SUCCESS;
								}
							}
						}
					} else {
						for (size_t j = 0; j < dym_array_size(&sourceVertices); ++j) {
							PKMER_VERTEX v = (PKMER_VERTEX)dym_array_get(&sourceVertices, j);
						
							for (size_t k = 0; k < v->degreeOut; ++k) {
								edge = kmer_vertex_get_succ_edge(v, k);
								if (!readRepeat)
									edge->Weight++;
							}
						}
					}

					if (ret == ERR_SUCCESS) {
						kmer_advance(sourceKMer, readSeq[i]);
						dym_array_exchange(&sourceVertices, &destVertices);
					}
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			dym_array_destroy(&destVertices);
		}

		dym_array_destroy(&sourceVertices);
		kmer_edge_table_destroy(edgeTable);
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

				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, Threshold, 1, kmetReference, &edge);
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

	ret = ERR_SUCCESS;
	for (size_t j = 0; j < ReadCount; ++j) {
		ret = _kmer_graph_parse_read(Graph, Reads, j);
		if (ret != ERR_SUCCESS)
			break;

		++Reads;
	}

	return ret;
}
