
#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__

#include "err.h"
#include "utils.h"
#include "reads.h"


typedef struct _PARSE_OPTIONS {
	boolean ConnectRefSeq;
	boolean ConnectReads;
	boolean MergeBubbles;
	boolean FixReads;
	boolean LinearShrink;
	boolean HelperVertices;
	uint32_t MissingEdgePenalty;
	uint32_t BackwardRefseqPenalty;
	uint32_t ReadThreshold;
} PARSE_OPTIONS, *PPARSE_OPTIONS;



ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const uint32_t Threshold);
ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, PONE_READ Reads, const size_t ReadCount, const size_t Threshold, const PARSE_OPTIONS *Options, PGEN_ARRAY_KMER_EDGE_PAIR PairArray);



#endif 
