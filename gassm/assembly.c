
#include <malloc.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
#include "assembly.h"











ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq)
{
	PKMER sourceKMer = NULL;
	PKMER destKMer = NULL;
	uint32_t kmerSize = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	kmerSize = kmer_graph_get_kmer_size(Graph);
	KMER_STACK_ALLOC(sourceKMer, kmerSize, RefSeq);
	kmer_back(sourceKMer, 'B');
	ret = kmer_graph_add_vertex(Graph, sourceKMer);
	if (ret == ERR_SUCCESS) {
		size_t refSeqLen = strlen(RefSeq);

		KMER_STACK_ALLOC(destKMer, kmerSize, RefSeq);
		kmer_back(destKMer, 'B');
		for (size_t i = kmerSize - 1; i <= refSeqLen; ++i) {
			if (i == refSeqLen)
				kmer_advance(destKMer, 'E');
			else kmer_advance(destKMer, RefSeq[i]);

			ret = kmer_graph_add_vertex(Graph, destKMer);
			if (ret == ERR_ALREADY_EXISTS)
				ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				ret = kmer_graph_add_edge(Graph, sourceKMer, destKMer, 0);
				if (ret == ERR_ALREADY_EXISTS)
					ret = ERR_SUCCESS;

				kmer_advance(sourceKMer, RefSeq[i]);
			}

			if (ret != ERR_SUCCESS)
				break;
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
	for (size_t i = kmerSize; i < readLen; ++i) {
		kmer_advance(destKMer, Read[i]);
		ret = kmer_graph_add_vertex(Graph, destKMer);
		if (ret == ERR_ALREADY_EXISTS)
			ret = ERR_SUCCESS;

		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_add_edge(Graph, sourceKMer, destKMer, 1);
			if (ret == ERR_ALREADY_EXISTS) {
				PKMER_EDGE edge = NULL;

				edge = kmer_graph_get_edge(Graph, sourceKMer, destKMer);
				edge->Weight++;
				ret = ERR_SUCCESS;
			}

			kmer_advance(sourceKMer, Read[i]);
		}

		if (ret != ERR_SUCCESS)
			break;
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
