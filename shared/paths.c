
#include "err.h"
#include "utils.h"
#include "general-queue.h"
#include "kmer-graph.h"
#include "dym-array.h"
#include "paths.h"


/************************************************************************/
/*                     HELPER FUNCTIONS                                 */
/************************************************************************/



static ERR_VALUE _place_new_path(PKMER_GRAPH_PATH Paths, size_t *PathCount, const size_t MaxPathCount, const double Weight, const char *Sequence, const size_t Length)
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
			break;
		}

		++p;
	}

	if (place == NULL && tmpPathsCount < MaxPathCount)
		place = Paths + tmpPathsCount;

	if (place != NULL) {
		place->Weight = Weight;
		if (place->Sequence != NULL)
			utils_free(place->Sequence);

		place->Length = Length;
		ret = utils_copy_string(Sequence, &place->Sequence);
		if (tmpPathsCount < MaxPathCount)
			++tmpPathsCount;
	}

	*PathCount = tmpPathsCount;

	return ret;
}


static void _stack_node_fill(PPATH_ELEMENT Node, const KMER_VERTEX *Vertex, const KMER_EDGE *Edge, const size_t EdgeIndex, const double Weight)
{
	Node->Edge = Edge;
	Node->EdgeIndex = EdgeIndex;
	Node->Vertex = Vertex;
	Node->Weight = Weight;

	return;
}



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
				_stack_node_fill(currentStack, Graph->StartingVertex, NULL, 0, 0);
				while (ret == ERR_SUCCESS && (currentStack != stack || stack->EdgeIndex < Graph->StartingVertex->degreeOut)) {
					if (currentStack->Length < EdgeCount) {
						if (currentStack->EdgeIndex < currentStack->Vertex->degreeOut) {
							PKMER km = (PKMER)dym_array_data(&v->Successors)[currentStack->EdgeIndex];
							PKMER_EDGE e = kmer_graph_get_edge(Graph, v->KMer, km);
							PKMER_VERTEX successor = (PKMER_VERTEX)kmer_table_get(Graph->VertexTable, km);

							++currentStack->EdgeIndex;
							if (e->MaxPassCount == 0 || e->PassCount < e->MaxPassCount) {
								++currentStack;
								_stack_node_fill(currentStack, successor, e, 0, ((e->Probability) + (currentStack - 1)->Weight));								
								if (tmpPathsCount < BestNumber || currentStack->Weight > tmpPaths[BestNumber - 1].Weight) {
									e->PassCount++;
									++seqIndex;
									seq[seqIndex] = km->Bases[kmerSize - 1];
									v = successor;
								} else --currentStack;
							}
						} else {
							if (currentStack->Edge != NULL)
								--currentStack->Edge->PassCount;
							
							--currentStack;
							seq[seqIndex] = '\0';
							--seqIndex;
							v = currentStack->Vertex;
						}
					} else {
						assert(currentStack->Length == EdgeCount);
						if (seq[seqIndex] == 'E')
							seq[seqIndex] = '\0';

						ret = _place_new_path(tmpPaths, &tmpPathsCount, BestNumber, currentStack->Weight, seq, currentStack->Length + kmerSize - 2);
						if (currentStack->Edge != NULL)
							--currentStack->Edge->PassCount;
						
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
