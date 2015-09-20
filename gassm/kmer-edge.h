
#ifndef __KMER_EDGE_H__
#define __KMER_EDGE_H__


#include <stdint.h>
#include "err.h"
#include "kmer.h"



typedef struct _KMER_EDGE {
	PKMER Source;
	PKMER Dest;
	long Weight;
	unsigned int Order;
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

ERR_VALUE kmer_edge_table_insert(PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
PKMER_EDGE kmer_edge_table_get(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);

size_t kmer_edge_hash(const PKMER_EDGE_TABLE Table, const PKMER Source, const PKMER Dest);
void kmer_edge_table_print(const PKMER_EDGE_TABLE Table);




#endif 
