
#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__



#include "err.h"
#include "reads.h"


/** Represents a part of a read. */
typedef struct _READ_PART {
	/** Part position within the reference sequence. */
	uint64_t Position;
	/** The read part sequence (it is not null-terminated). */
	const char *ReadSequence;
	/** Length of the part sequence, in bases. */
	size_t ReadSequenceLength;
	/** Position of the part relative to the start of the read. */
	size_t Offset;
} READ_PART, *PREAD_PART;

GEN_ARRAY_TYPEDEF(READ_PART);
GEN_ARRAY_IMPLEMENTATION(READ_PART)


ERR_VALUE kmer_graph_parse_ref_sequence(PKMER_GRAPH Graph, const char *RefSeq, const size_t RefSeqLen, const uint32_t Threshold);
ERR_VALUE kmer_graph_parse_reads(PKMER_GRAPH Graph, const struct _ONE_READ *Reads, const size_t ReadCount, PGEN_ARRAY_KMER_EDGE_PAIR PairArray);
ERR_VALUE read_split(const ONE_READ *Read, PGEN_ARRAY_READ_PART PartArray);



#endif 
