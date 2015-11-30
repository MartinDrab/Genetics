
#include "err.h"
#include "utils.h"
#include "general-queue.h"
#include "kmer-graph.h"
#include "dym-array.h"
#include "paths.h"

/************************************************************************/
/*                      PUBLIC FUNCTIONS                                */
/************************************************************************/




ERR_VALUE kmer_graph_find_best_paths(PKMER_GRAPH Graph, const size_t BestNumber, const size_t EdgeCount, PKMER_GRAPH_PATH *Paths, size_t *Count)
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

	ret = utils_calloc(EdgeCount + kmerSize + 1, sizeof(char), (void **)&seq);
	if (ret == ERR_SUCCESS) {
		v = Graph->StartingVertex;
		memset(seq, 0, (EdgeCount + kmerSize + 1)*sizeof(char));
		memcpy(seq, v->KMer->Bases + 1, (kmerSize - 1)*sizeof(char));
		seqIndex = kmerSize - 2;
		ret = utils_calloc(BestNumber, sizeof(KMER_GRAPH_PATH), (void **)&tmpPaths);
		if (ret == ERR_SUCCESS) {
			memset(tmpPaths, 0, BestNumber*sizeof(KMER_GRAPH_PATH));
			ret = utils_calloc(EdgeCount + 1, sizeof(PATH_ELEMENT), (void **)&stack);
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < EdgeCount + 1; ++i) {
					memset(stack + i, 0, sizeof(PATH_ELEMENT));
					stack[i].Length = i;
				}

				currentStack = stack;
				currentStack->EdgeIndex = 0;
				currentStack->Vertex = Graph->StartingVertex;
				currentStack->Weight = 1;
				while (ret == ERR_SUCCESS && (currentStack != stack || stack->EdgeIndex < Graph->StartingVertex->degreeOut)) {
					if (currentStack->Length < EdgeCount) {
						if (currentStack->EdgeIndex < currentStack->Vertex->degreeOut) {
							PKMER km = (PKMER)dym_array_data(&v->Successors)[currentStack->EdgeIndex];
							PKMER_EDGE e = kmer_graph_get_edge(Graph, v->KMer, km);
							PKMER_VERTEX successor = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, km);

							++currentStack->EdgeIndex;
							++currentStack;
							currentStack->EdgeIndex = 0;
							currentStack->Vertex = successor;
							currentStack->Weight = (e->Probability) * (currentStack - 1)->Weight;
							if (currentStack->Weight == 0)
								currentStack->Weight = (currentStack - 1)->Weight;

							if (tmpPathsCount < BestNumber || currentStack->Weight > tmpPaths[tmpPathsCount - 1].Weight) {
								++seqIndex;
								seq[seqIndex] = km->Bases[kmerSize - 1];
								v = successor;
							} else --currentStack;
						} else {
							--currentStack;
							seq[seqIndex] = '\0';
							--seqIndex;
							v = currentStack->Vertex;
						}
					} else {
						assert(currentStack->Length == EdgeCount);
						if (seq[seqIndex] == 'E')
							seq[seqIndex] = '\0';

						{
							double w = currentStack->Weight;
							PKMER_GRAPH_PATH place = NULL;
							PKMER_GRAPH_PATH p = tmpPaths;

							for (size_t i = 0; i < tmpPathsCount; ++i) {
								if (w >= p->Weight) {
									place = p;
									for (size_t j = tmpPathsCount; j > i; --j) {
										if (j == BestNumber)
											continue;

										tmpPaths[j] = tmpPaths[j - 1];
									}

									place->Sequence = NULL;
									break;
								}

								++p;
							}

							if (place == NULL && tmpPathsCount < BestNumber)
								place = tmpPaths + tmpPathsCount;

							if (place != NULL) {
								place->Weight = w;
								if (place->Sequence != NULL)
									utils_free(place->Sequence);

								ret = utils_copy_string(seq, &place->Sequence);
								if (tmpPathsCount < BestNumber)
									++tmpPathsCount;
							}
						}

						--currentStack;
						v = currentStack->Vertex;
						seq[seqIndex] = '\0';
						--seqIndex;
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


void kmer_graph_paths_free(PKMER_GRAPH_PATH Paths, const size_t Count)
{
	for (size_t i = 0; i < Count; ++i) {
		if (Paths[i].Sequence != NULL)
			utils_free(Paths[i].Sequence);
	}

	utils_free(Paths);

	return;
}
