
#ifndef __GASSM_INPUT_FILE_H__
#define __GASSM_INPUT_FILE_H__


#include "gen_dym_array.h"
#include "pointer_array.h"


typedef enum _EActiveRegionType {
	artUnknown,
	artValid,
} EActiveRegionType, *PEActiveRegionType;

typedef struct _ACTIVE_REGION {
	EActiveRegionType Type;
	uint64_t Offset;
	uint64_t Length;
	const char *Sequence;
} ACTIVE_REGION, *PACTIVE_REGION;

GEN_ARRAY_TYPEDEF(ACTIVE_REGION);
GEN_ARRAY_IMPLEMENTATION(ACTIVE_REGION)
POINTER_ARRAY_TYPEDEF(ACTIVE_REGION);
POINTER_ARRAY_IMPLEMENTATION(ACTIVE_REGION)

typedef struct _FASTA_FILE {
	char *FileData;
	size_t DataLength;
	char *CurrentPointer;
} FASTA_FILE, *PFASTA_FILE;

typedef struct _REFSEQ_DATA {
	const char *Sequence;
	size_t Length;
	uint64_t StartPos;
	const char *Name;
} REFSEQ_DATA, *PREFSEQ_DATA;


ERR_VALUE fasta_load(const char *FileName, PFASTA_FILE FastaRecord);
ERR_VALUE fasta_read_seq(PFASTA_FILE FastaRecord, PREFSEQ_DATA Data);
void fasta_free_seq(PREFSEQ_DATA Data);
void fasta_free(PFASTA_FILE FastaRecord);

ERR_VALUE input_get_reads(const char *Filename, const char *InputType, PONE_READ *Reads, size_t *ReadCount);
ERR_VALUE input_filter_reads(const uint32_t KMerSize, const ONE_READ *Source, const size_t SourceCount, const uint64_t RegionStart, const size_t RegionLength, PGEN_ARRAY_ONE_READ NewReads);
void input_free_filtered_reads(PONE_READ Reads, size_t Count);
void input_filter_bad_reads(PONE_READ Reads, size_t *Count, const uint8_t MinQuality, boolean UseCIGAR);
void input_sort_reads(PONE_READ Reads, const size_t Count);
void input_free_reads(PONE_READ Reads, const size_t Count);
ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count);
ERR_VALUE input_get_region_by_offset(const PACTIVE_REGION Regions, const size_t Count, const uint64_t Offset, size_t *Index, uint64_t *RegionOffset);
void input_free_regions(PACTIVE_REGION Regions, const size_t Count);






#endif 
