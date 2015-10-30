
#include <malloc.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
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

		ret = kmer_graph_add_vertex(Graph, tmpKMer);
		if (ret == ERR_SUCCESS) {
			*EdgeLength = j - Index;
			kmer_init_from_kmer(DestKMer, tmpKMer);
			break;
		} else if (ret != ERR_ALREADY_EXISTS)
			break;
	}

	return ret;
}

/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/

ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq)
{
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	const uint32_t kmerSize = kmer_graph_get_kmer_size(Graph);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	KMER_STACK_ALLOC(sourceKMer, kmerSize, RefSeq);
	kmer_back(sourceKMer, 'B');
	ret = kmer_graph_add_vertex(Graph, sourceKMer);
	if (ret == ERR_SUCCESS) {
		size_t refSeqLen = strlen(RefSeq);

		KMER_STACK_ALLOC(destKMer, kmerSize, RefSeq);
		kmer_back(destKMer, 'B');
		size_t i = kmerSize - 1;
		while (i <= refSeqLen) {
			uint32_t addEdgeLength = 0;
			
			if (i == refSeqLen)
				kmer_advance(destKMer, 'E');
			else kmer_advance(destKMer, RefSeq[i]);

			ret = kmer_graph_add_vertex(Graph, destKMer);
			if (ret == ERR_ALREADY_EXISTS) {
				ret = _attempt_skip_vertices(Graph, destKMer, &addEdgeLength, kmerSize, RefSeq, i, refSeqLen, TRUE);
				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;
			}

			if (ret == ERR_SUCCESS) {
				PKMER_EDGE edge = NULL;

				ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, 0, addEdgeLength + 1, &edge);
				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

				if (addEdgeLength > 0) {
					i += addEdgeLength;
					kmer_init_from_kmer(sourceKMer, destKMer);
				} else kmer_advance(sourceKMer, RefSeq[i]);
			}

			if (ret != ERR_SUCCESS)
				break;

			++i;
		}
	}

	return ret;
}


ERR_VALUE kmer_graph_parse_read(PKMER_GRAPH Graph, const char *Read)
{
	size_t readLen = 0;
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	uint32_t kmerSize = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	readLen = strlen(Read);
	kmerSize = kmer_graph_get_kmer_size(Graph);
	KMER_STACK_ALLOC(sourceKMer, kmerSize, Read);
	ret = kmer_graph_add_vertex(Graph, sourceKMer);
	if (ret == ERR_ALREADY_EXISTS)
		ret = ERR_SUCCESS;

	KMER_STACK_ALLOC(destKMer, kmerSize, Read);
	size_t i = kmerSize;
	while (i < readLen) {
		uint32_t addEdgeLength = 0;

		kmer_advance(destKMer, Read[i]);
		ret = kmer_graph_add_vertex(Graph, destKMer);
		if (ret == ERR_ALREADY_EXISTS) {
//			ret = _attempt_skip_vertices(Graph, destKMer, &addEdgeLength, kmerSize, Read, i, readLen, FALSE);
			if (ret == ERR_ALREADY_EXISTS)
				ret = ERR_SUCCESS;
		}

		if (ret == ERR_SUCCESS) {
			PKMER_EDGE edge = NULL;

			ret = kmer_graph_add_edge_ex(Graph, sourceKMer, destKMer, 0, addEdgeLength + 1, &edge);
			if (ret == ERR_ALREADY_EXISTS) {
				edge->Weight++;
				ret = ERR_SUCCESS;
			}

			if (addEdgeLength > 0) {
				i += addEdgeLength;
				kmer_init_from_kmer(sourceKMer, destKMer);
			} else kmer_advance(sourceKMer, Read[i]);
		}

		if (ret != ERR_SUCCESS)
			break;

		++i;
	}

	return ret;
}

ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const char **Reads, const size_t ReadCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	for (size_t j = 0; j < ReadCount; ++j) {
		const char *currentRead = Reads[j];

		ret = kmer_graph_parse_read(Graph, currentRead);
		if (ret != ERR_SUCCESS)
			break;
	}

	return ret;
}
