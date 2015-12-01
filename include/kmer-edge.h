
#ifndef __KMER_EDGE_H__
#define __KMER_EDGE_H__


#include <stdint.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"


typedef enum _EKMerEdgeType {
	kmetReference,
	kmetRead,
} EKMerEdgeType, *PEKMerEdgeType;


/** Represents one edge in a kmer graph. */
typedef struct _KMER_EDGE {
	boolean Deleted;
	/** A kmer the edge starts in. */
	PKMER Source;
	/** A kmer the edge leads to. */
	PKMER Dest;
	/** Number of reads going through the edge. */
	long Weight;
	/** Number of kmers the edge skips in order to avoid cycles. Maximum value is
	    the kmer size minus one. */
	uint32_t Length;
	/** Edge creation order. */
	unsigned int Order;
	EKMerEdgeType Type;
	double Probability;
	uint32_t PassCount;
	uint32_t MaxPassCount;
} KMER_EDGE, *PKMER_EDGE;



typedef struct _KMER_EDGE_TABLE {
	size_t Size;
	size_t Inverse;
	size_t X;
	size_t PowX;
	size_t KMerSize;
	PKMER_EDGE Entries;
} KMER_EDGE_TABLE, *PKMER_EDGE_TABLE;




ERR_VALUE kmer_edge_table_create(const size_t KMerSize, const size_t X, const size_t Size, PKMER_EDGE_TABLE *Table);
void kmer_edge_table_destroy(PKMER_EDGE_TABLE Table);
ERR_VALUE kmer_edge_table_extend(PKMER_EDGE_TABLE Table);
ERR_VALUE kmer_edge_table_copy(const PKMER_EDGE_TABLE Source, PKMER_EDGE_TABLE * Copied);

ERR_VALUE kmer_edge_table_insert_ex(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest, PKMER_EDGE *Edge);
ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
ERR_VALUE kmer_edge_table_delete(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
void kmer_edge_table_delete_by_entry(PKMER_EDGE_TABLE Table, PKMER_EDGE Entry);
PKMER_EDGE kmer_edge_table_get(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
ERR_VALUE kmer_edge_table_first(const PKMER_EDGE_TABLE Table, PKMER_EDGE *Slot);
ERR_VALUE kmer_edge_table_next(const PKMER_EDGE_TABLE Table, const PKMER_EDGE Current, PKMER_EDGE *Next);

size_t kmer_edge_hash(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
void kmer_edge_table_print(FILE *Stream, const PKMER_EDGE_TABLE Table);




#endif 
