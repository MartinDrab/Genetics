
#ifndef __KMER_READ_INFO_H__
#define __KMER_READ_INFO_H__


#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"


typedef struct _READ_INFO_ENTRY {
	size_t ReadIndex;
	size_t ReadPosition;
	uint8_t Quality;
} READ_INFO_ENTRY, *PREAD_INFO_ENTRY;


GEN_ARRAY_TYPEDEF(READ_INFO_ENTRY);

GEN_ARRAY_IMPLEMENTATION(READ_INFO_ENTRY)

typedef struct _READ_INFO {
	GEN_ARRAY_TYPE(READ_INFO_ENTRY) Array;
} READ_INFO, *PREAD_INFO;

POINTER_ARRAY_TYPEDEF(READ_INFO);
POINTER_ARRAY_IMPLEMENTATION(READ_INFO);



ERR_VALUE read_info_to_indices(const READ_INFO *Info, PGEN_ARRAY_size_t Indices);
void read_info_sort(PREAD_INFO ReadInfo);
size_t read_info_weight(const READ_INFO *Info, const uint8_t *QualityTable);
void read_info_init(PREAD_INFO Info);
ERR_VALUE read_info_assign(PREAD_INFO Info, const GEN_ARRAY_READ_INFO_ENTRY *Array);
void read_info_finit(PREAD_INFO Info);
ERR_VALUE read_info_copy(PREAD_INFO Dest, const READ_INFO *Source);
ERR_VALUE read_info_add_array(PREAD_INFO Info, const GEN_ARRAY_READ_INFO_ENTRY *Array);
ERR_VALUE read_info_add(PREAD_INFO Info, const size_t ReadIndex, const size_t ReadPosition, const uint8_t Quality);
size_t read_info_remove(PREAD_INFO Info, const size_t ReadIndex, const size_t ReadPosition);
void read_info_remove_2(PREAD_INFO Info, const size_t ReadIndex);
ERR_VALUE read_info_intersection(PREAD_INFO Info1, PREAD_INFO Info2, GEN_ARRAY_PTYPE(READ_INFO_ENTRY) Intersection, const size_t ReadDistance);
void read_info_subtract(PREAD_INFO Info, const GEN_ARRAY_TYPE(READ_INFO_ENTRY) *Subtrahend, const size_t Distance);
ERR_VALUE read_info_union(PREAD_INFO Target, const GEN_ARRAY_READ_INFO_ENTRY *Source);
ERR_VALUE read_info_merge(PREAD_INFO Dest, const READ_INFO *Info1, const READ_INFO *Info2);
void read_info_clear(PREAD_INFO Info);

#define	read_info_get_count(aInfo)								gen_array_size(&(aInfo)->Array)
#define read_info_get_entry(aInfo, aIndex)						dym_array_item_READ_INFO_ENTRY((&(aInfo)->Array), (aIndex))


#endif 
