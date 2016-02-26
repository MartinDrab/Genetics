
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


static ERR_VALUE _kmer_graph_parse_read(PKMER_GRAPH Graph, const ONE_READ *Read, const size_t ReadIndex)
{
	char *readSeq = NULL;
	PKMER_EDGE lastEdge = NULL;
	size_t readLen = 0;
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	uint32_t kmerSize = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_EDGE_TABLE edgeTable = NULL;

	kmerSize = kmer_graph_get_kmer_size(Graph);
	ret = kmer_edge_table_create(kmerSize, 2, 37, NULL, &edgeTable);
	if (ret == ERR_SUCCESS) {
		readSeq = Read->ReadSequence;
		readLen = Read->ReadSequenceLen;
		KMER_STACK_ALLOC(sourceKMer, kmerSize, readSeq);
		ret = kmer_graph_add_vertex(Graph, sourceKMer, kmvtRead);
		if (ret == ERR_ALREADY_EXISTS)
			ret = ERR_SUCCESS;

		KMER_STACK_ALLOC(destKMer, kmerSize, readSeq);
		size_t i = kmerSize;
		while (i < readLen) {
			kmer_advance(destKMer, readSeq[i]);
			ret = kmer_graph_add_vertex(Graph, destKMer, kmvtRead);
			if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				boolean readRepeat = FALSE;

				ret = kmer_edge_table_insert(edgeTable, sourceKMer, destKMer, (void *)TRUE);
				if (ret == ERR_ALREADY_EXISTS) {
					readRepeat = TRUE;
					ret = ERR_SUCCESS;
				}
				
				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, 1, 1, kmetRead, &edge);
				if (ret == ERR_ALREADY_EXISTS) {
					if (!readRepeat)
						edge->Weight++;
					
					ret = ERR_SUCCESS;
				}

				if (ret == ERR_SUCCESS) {
					ret = kmer_edge_add_read(edge, ReadIndex, i);
					if (ret == ERR_SUCCESS) {
						if ((lastEdge != NULL && lastEdge->Type == kmetRead) || edge->Type == kmetRead) {
							PKMER_VERTEX v = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, sourceKMer);
							boolean found = FALSE;

							assert(v != NULL);
							for (size_t k = 0; k < v->PassCount; ++k) {
								if (v->Passes[k].Incomming == lastEdge &&
									v->Passes[k].Outgoing == edge) {
									found = TRUE;
									break;
								}

							}

							if (!found)
								ret = kmer_vertex_add_pass(v, lastEdge, edge, vptRead);
						}

						if (ret == ERR_SUCCESS) {
							lastEdge = edge;
							kmer_advance(sourceKMer, readSeq[i]);
						}
					}
				}
			}

			if (ret != ERR_SUCCESS)
				break;

			++i;
		}

		kmer_edge_table_destroy(edgeTable);
	}

	return ret;
}

/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const uint32_t Threshold)
{
	EKMerVertexType vertexType = kmvtRefSeqMiddle;
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	KMER_STACK_ALLOC(sourceKMer, kmerSize, RefSeq);
	kmer_back(sourceKMer, 'B');
	ret = kmer_graph_add_vertex(Graph, sourceKMer, kmvtRefSeqStart);
	if (ret == ERR_SUCCESS) {
		PKMER_EDGE lastEdge = NULL;
		
		kmer_graph_set_starting_vertex(Graph, sourceKMer);
		KMER_STACK_ALLOC(destKMer, kmerSize, RefSeq);
		kmer_back(destKMer, 'B');
		size_t i = kmerSize - 1;
		while (i <= RefSeqLen) {
			if (i == RefSeqLen) {
				kmer_advance(destKMer, 'E');
				vertexType = kmvtRefSeqEnd;
			} else kmer_advance(destKMer, RefSeq[i]);

			boolean destVertexExists = FALSE;
			ret = kmer_graph_add_vertex(Graph, destKMer, vertexType);
			if (ret == ERR_ALREADY_EXISTS) {
				ret = ERR_SUCCESS;
				destVertexExists = TRUE;
			}

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				PKMER_VERTEX sourceVertex = NULL;

				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, Threshold, 1, kmetReference, &edge);
				if (ret == ERR_SUCCESS && destVertexExists)
					++Graph->NumberOfBackwardEdges;
				
				if (ret == ERR_SUCCESS || ret == ERR_ALREADY_EXISTS) {
					++edge->MaxPassCount;
					ret = ERR_SUCCESS;
				}

				sourceVertex = kmer_table_get(Graph->VertexTable, sourceKMer);
				ret = kmer_vertex_add_pass(sourceVertex, lastEdge, edge, vptRefSeq);
				if (ret == ERR_SUCCESS) {
					lastEdge = edge;
					kmer_advance(sourceKMer, RefSeq[i]);
				}
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
