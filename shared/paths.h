
#ifndef __KMER_GRAPH_PATHS_H__
#define __KMER_GRAPH_PATHS_H__


#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"




typedef struct _PATH_ELEMENT {
	double Weight;
	PKMER_VERTEX Vertex;
	size_t EdgeIndex;
	PKMER_EDGE Edge;
} PATH_ELEMENT, *PPATH_ELEMENT;

typedef struct _KMER_GRAPH_PATH {
	char *Sequence;
	size_t Length;
	double Weight;
} KMER_GRAPH_PATH, *PKMER_GRAPH_PATH;



ERR_VALUE kmer_graph_find_best_paths(PKMER_GRAPH Graph, const size_t BestNumber, const size_t EdgeCount, PKMER_GRAPH_PATH *Paths, size_t *Count);
void kmer_graph_paths_free(PKMER_GRAPH_PATH Paths, const size_t Count);



#endif 
