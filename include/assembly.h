
#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__



#include "err.h"
#include "reads.h"




ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const boolean SkipVertices);
ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount, const boolean SkipVertices);






#endif 
