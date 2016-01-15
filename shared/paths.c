
#include <math.h>
#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "general-queue.h"
#include "kmer-graph.h"
#include "dym-array.h"
#include "ssw.h"
#include "paths.h"


/************************************************************************/
/*                     HELPER FUNCTIONS                                 */
/************************************************************************/



static ERR_VALUE _place_new_path(PKMER_GRAPH_PATH Paths, size_t *PathCount, const size_t MaxPathCount, const double Weight, const char *Sequence, const size_t Length, const SSW_STATISTICS *SSWStats, const char *OpString)
{
	PKMER_GRAPH_PATH p = Paths;
	PKMER_GRAPH_PATH place = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	size_t tmpPathsCount = 0;

	ret = ERR_SUCCESS;
	tmpPathsCount = *PathCount;
	for (size_t i = 0; i < tmpPathsCount; ++i) {
		if (Weight >= p->Weight) {
			place = p;
			for (size_t j = tmpPathsCount; j > i; --j) {
				if (j == MaxPathCount)
					continue;

				Paths[j] = Paths[j - 1];
			}

			place->Sequence = NULL;
			place->OpString = NULL;
			break;
		}

		++p;
	}

	if (place == NULL && tmpPathsCount < MaxPathCount)
		place = Paths + tmpPathsCount;

	if (place != NULL) {
		place->SSWStatistics = *SSWStats;
		place->Weight = Weight;
		if (place->Sequence != NULL) {
			utils_free(place->Sequence);
			place->Sequence = NULL;
		}

		if (place->OpString != NULL) {
			utils_free(place->OpString);
			place->OpString = NULL;
		}

		place->Length = Length;
		ret = utils_copy_string(Sequence, &place->Sequence);
		if (ret == ERR_SUCCESS) {
			ret = utils_copy_string(OpString, &place->OpString);
			if (ret == ERR_SUCCESS) {
				if (tmpPathsCount < MaxPathCount)
					++tmpPathsCount;
			}

			if (ret != ERR_SUCCESS) {
				utils_free(place->Sequence);
				place->Sequence = NULL;
			}
		}
	}

	*PathCount = tmpPathsCount;

	return ret;
}


static void _stack_node_fill(PPATH_ELEMENT Node, const PKMER_VERTEX Vertex, const PKMER_EDGE Edge, const size_t EdgeIndex, const double Weight)
{
	Node->Edge = Edge;
	Node->EdgeIndex = EdgeIndex;
	Node->Vertex = Vertex;
	Node->Weight = Weight;

	return;
}


#define _path_stack_pop()								\
{														\
	if (shortcut == NULL) {								\
		currentLength--;								\
		if (currentStack->Edge != NULL)					\
			--currentStack->Edge->PassCount;			\
														\
		seq[seqIndex] = '\0';							\
		--seqIndex;										\
	} else {											\
		seqIndex -= shortcut->Length;					\
		seq[seqIndex + 1] = '\0';						\
		currentLength -= shortcut->Length;				\
		shortcut->PassCount--;							\
	}													\
														\
	--currentStack;										\
	v = currentStack->Vertex;							\
}														\

/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/


ERR_VALUE kmer_graph_find_best_paths(PKMER_GRAPH Graph, const char *RegionStart, const size_t RegionLength, const size_t MaxBestPaths, const PATH_SCORING *ScoreWeights, PKMER_GRAPH_PATH *Paths, size_t *Count)
{
	const uint32_t kmerSize = Graph->KMerSize;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PKMER_GRAPH_PATH tmpPaths = NULL;
	size_t tmpPathsCount = 0;
	PKMER_VERTEX v = NULL;
	PPATH_ELEMENT stack = NULL;
	PPATH_ELEMENT currentStack = NULL;
	char *seq = NULL;
	size_t seqIndex = 0;
	size_t currentLength = 0;
	const size_t stackSize = (RegionLength * 2);
	
	ret = utils_calloc(stackSize + 1, sizeof(char), (void **)&seq);
	if (ret == ERR_SUCCESS) {
		v = Graph->StartingVertex;
		memset(seq, 0, (stackSize + 1)*sizeof(char));
		memcpy(seq, v->KMer->Bases + 1, (kmerSize - 1)*sizeof(char));
		seqIndex = kmerSize - 2;
		ret = utils_calloc(MaxBestPaths, sizeof(KMER_GRAPH_PATH), (void **)&tmpPaths);
		if (ret == ERR_SUCCESS) {
			memset(tmpPaths, 0, MaxBestPaths*sizeof(KMER_GRAPH_PATH));
			ret = utils_calloc(stackSize, sizeof(PATH_ELEMENT), (void **)&stack);
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < stackSize; ++i)
					memset(stack + i, 0, sizeof(PATH_ELEMENT));

				currentStack = stack;
				_stack_node_fill(currentStack, Graph->StartingVertex, NULL, 0, 0);
				while (ret == ERR_SUCCESS && (currentStack != stack || stack->EdgeIndex < Graph->StartingVertex->degreeOut)) {
					PKMER_GRAPH_SHORTCUT shortcut = NULL;
					
					if (currentStack->Edge != NULL)
						shortcut = (PKMER_GRAPH_SHORTCUT)currentStack->Edge->Shortcut;

					if (seq[seqIndex] == 'E') {
						size_t opStringLen = 0;
						char *opString = NULL;

						seq[seqIndex] = '\0';
						ret = ssw_clever(seq, seqIndex - 1, RegionStart, RegionLength, ScoreWeights->Match, ScoreWeights->Mismatch, ScoreWeights->Indel, &opString, &opStringLen);
						if (ret == ERR_SUCCESS) {
							SSW_STATISTICS sswStats;
							double indelDiffWeight = 0;

							opstring_statistics(opString, opStringLen, &sswStats);
							indelDiffWeight = (sswStats.TotalInsertions - sswStats.TotalDeletions);
							indelDiffWeight = indelDiffWeight / (1 + sswStats.TotalInsertions + sswStats.TotalDeletions);
							indelDiffWeight = 1 - indelDiffWeight*indelDiffWeight;
							
							currentStack->Weight += ScoreWeights->SSWWeight*log(indelDiffWeight);
							ret = _place_new_path(tmpPaths, &tmpPathsCount, MaxBestPaths, currentStack->Weight, seq, currentLength + kmerSize - 2, &sswStats, opString);
							utils_free(opString);
						}

						_path_stack_pop();
					} else if (currentLength < stackSize) {
						if (currentStack->EdgeIndex < currentStack->Vertex->degreeOut) {
							PKMER_EDGE e = kmer_vertex_get_succ_edge(v, currentStack->EdgeIndex);
							PKMER_VERTEX successor = e->Dest;

							shortcut = (PKMER_GRAPH_SHORTCUT)e->Shortcut;
							++currentStack->EdgeIndex;
							if ((shortcut == NULL && e->PassCount < e->MaxPassCount) ||
								(shortcut != NULL && shortcut->PassCount < shortcut->MaxPassCount)) {
								++currentStack;
								_stack_node_fill(currentStack, ((shortcut != NULL) ? shortcut->EndVertex : successor), e, 0, ((ScoreWeights->EdgeProbabilityWeight*e->Probability) + (currentStack - 1)->Weight));
								if (tmpPathsCount < MaxBestPaths || currentStack->Weight > tmpPaths[MaxBestPaths - 1].Weight) {
									if (shortcut != NULL) {
										shortcut->PassCount++;
										memcpy(seq + seqIndex + 1, shortcut->Sequence, shortcut->Length*sizeof(char));
										seqIndex += shortcut->Length;
										currentLength += shortcut->Length;
									} else {
										e->PassCount++;
										++seqIndex;
										++currentLength;
										seq[seqIndex] = successor->KMer->Bases[kmerSize - 1];
									}

									v = currentStack->Vertex;
								} else {
//									fprintf(stderr, "Path pruned (weight = %lf, length = %u, path count = %u)\n", currentStack->Weight, (uint32_t)seqIndex, tmpPathsCount);
									--currentStack;
								}
							}
						} else {
							_path_stack_pop();
						}
					} else {
						_path_stack_pop();
					}

				}

				utils_free(stack);
			}

			if (ret == ERR_SUCCESS) {
				*Paths = tmpPaths;
				*Count = tmpPathsCount;
			}

			if (ret != ERR_SUCCESS)
				utils_free(tmpPaths);
		}
	
		utils_free(seq);
	}

	return ret;
}


void kmer_graph_path_write(FILE *Stream, const KMER_GRAPH_PATH *Path)
{

	return;
}


void kmer_graph_path_free(PKMER_GRAPH_PATH Path)
{
	if (Path->OpString != NULL)
		utils_free(Path->OpString);

	if (Path->Sequence != NULL)
		utils_free(Path->Sequence);

	return;
}

void kmer_graph_paths_free(PKMER_GRAPH_PATH Paths, const size_t Count)
{
	for (size_t i = 0; i < Count; ++i)
		kmer_graph_path_free(Paths + i);

	utils_free(Paths);

	return;
}
