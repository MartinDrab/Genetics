
#ifndef __KMER_GRAPH_PATHS_H__
#define __KMER_GRAPH_PATHS_H__


#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "kmer-graph.h"
#include "ssw.h"



typedef struct _PATH_ELEMENT {
	double Weight;
	PKMER_VERTEX Vertex;
	size_t EdgeIndex;
	PKMER_EDGE Edge;
} PATH_ELEMENT, *PPATH_ELEMENT;

typedef struct _KMER_GRAPH_PATH {
	char *Sequence;
	char *OpString;
	size_t Length;
	double Weight;
	SSW_STATISTICS SSWStatistics;
} KMER_GRAPH_PATH, *PKMER_GRAPH_PATH;

typedef struct _PATH_SCORING {
	double EdgeProbabilityWeight;
	double SSWWeight;
	int Match;
	int Mismatch;
	int Indel;
} PATH_SCORING, *PPATH_SCORING;



ERR_VALUE kmer_graph_find_best_paths(PKMER_GRAPH Graph, const char *RegionStart, const size_t RegionLength, const size_t MaxBestPaths, const PATH_SCORING *ScoreWeights, PKMER_GRAPH_PATH *Paths, size_t *Count);

void kmer_graph_path_write(FILE *Stream, const KMER_GRAPH_PATH *Path);
void kmer_graph_path_free(PKMER_GRAPH_PATH Path);
void kmer_graph_paths_free(PKMER_GRAPH_PATH Paths, const size_t Count);



#endif 
