
#ifndef __PAIRED_READS_H__
#define __PAIRED_READS_H__


#include "khash.h"
#include "reads.h"


ERR_VALUE paired_reads_insert(const ONE_READ *Read);
ERR_VALUE paired_reads_insert_array(const ONE_READ *Reads, const size_t Count);
ERR_VALUE paired_reads_first(khiter_t *Iterator, PPOINTER_ARRAY_ONE_READ *Reads);
ERR_VALUE paired_reads_next(khiter_t Iterator, khiter_t *NewIt, PPOINTER_ARRAY_ONE_READ *Reads);
void paired_reads_connect(void);
void paired_reads_print(FILE *Stream);

ERR_VALUE paired_reads_init(void);
void paired_reads_finit(void);





#endif
