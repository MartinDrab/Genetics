
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

static ERR_VALUE _attempt_skip_vertices(PKMER_GRAPH Graph, PKMER DestKMer, uint32_t *EdgeLength, const uint32_t KMerSize, const char *RefSeq, const size_t Index, const size_t RefSeqLen, const boolean IsRefSequence)
{
	PKMER tmpKMer = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t maxSize = min(Index + KMerSize, RefSeqLen + 1*(IsRefSequence));

	ret = ERR_SUCCESS;
	KMER_STACK_ALLOC(tmpKMer, KMerSize, RefSeq + Index);
	kmer_init_from_kmer(tmpKMer, DestKMer);
	for (size_t j = Index + 1; j < maxSize; ++j) {
		if (IsRefSequence && j == RefSeqLen)
			kmer_advance(tmpKMer, 'E');
		else kmer_advance(tmpKMer, RefSeq[j]);

		ret = kmer_graph_add_vertex(Graph, tmpKMer, (IsRefSequence) ? kmvtRefSeqMiddle : kmvtRead);
		if (ret == ERR_SUCCESS) {
			*EdgeLength = j - Index;
			kmer_init_from_kmer(DestKMer, tmpKMer);
			break;
		} else if (ret != ERR_ALREADY_EXISTS)
			break;
	}

	return ret;
}


static ERR_VALUE _kmer_graph_parse_read(PKMER_GRAPH Graph, const ONE_READ *Read, const boolean SkipVertices)
{
	char *readSeq = NULL;
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
			uint32_t addEdgeLength = 0;

			kmer_advance(destKMer, readSeq[i]);
			ret = kmer_graph_add_vertex(Graph, destKMer, kmvtRead);
			if (ret == ERR_ALREADY_EXISTS) {
				if (SkipVertices)
					ret = _attempt_skip_vertices(Graph, destKMer, &addEdgeLength, kmerSize, readSeq, i, readLen, FALSE);

				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				boolean readRepeat = FALSE;

				ret = kmer_edge_table_insert(edgeTable, sourceKMer, destKMer, (void *)TRUE);
				if (ret == ERR_ALREADY_EXISTS) {
					readRepeat = TRUE;
					ret = ERR_SUCCESS;
				}
				
				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, 1, addEdgeLength + 1, kmetRead, &edge);
				if (ret == ERR_ALREADY_EXISTS) {
					if (!readRepeat)
						edge->Weight++;
					
					ret = ERR_SUCCESS;
				}

				if (addEdgeLength > 0) {
					i += addEdgeLength;
					kmer_init_from_kmer(sourceKMer, destKMer);
				} else kmer_advance(sourceKMer, readSeq[i]);
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

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const boolean SkipVertices, const uint32_t Threshold)
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
			uint32_t addEdgeLength = 0;

			if (i == RefSeqLen) {
				kmer_advance(destKMer, 'E');
				vertexType = kmvtRefSeqEnd;
			} else kmer_advance(destKMer, RefSeq[i]);

			ret = kmer_graph_add_vertex(Graph, destKMer, vertexType);
			if (ret == ERR_ALREADY_EXISTS) {
				if (SkipVertices)
					ret = _attempt_skip_vertices(Graph, destKMer, &addEdgeLength, kmerSize, RefSeq, i, RefSeqLen, TRUE);

				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;
				PKMER_VERTEX sourceVertex = NULL;

				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, Threshold, addEdgeLength + 1, kmetReference, &edge);
				if (ret == ERR_SUCCESS || ret == ERR_ALREADY_EXISTS) {
					++edge->MaxPassCount;
					ret = ERR_SUCCESS;
				}

				sourceVertex = kmer_table_get(Graph->VertexTable, sourceKMer);
				ret = kmer_vertex_add_pass(sourceVertex, lastEdge, edge);
				if (ret == ERR_SUCCESS) {
					lastEdge = edge;
					if (addEdgeLength > 0) {
						i += addEdgeLength;
						kmer_init_from_kmer(sourceKMer, destKMer);
					} else kmer_advance(sourceKMer, RefSeq[i]);
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


ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount, const boolean SkipVertices)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	for (size_t j = 0; j < ReadCount; ++j) {
		ret = _kmer_graph_parse_read(Graph, Reads, SkipVertices);
		if (ret != ERR_SUCCESS)
			break;

		++Reads;
	}

	return ret;
}
