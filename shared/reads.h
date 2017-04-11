
#ifndef __READS_H__
#define __READS_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"
#include "pointer_array.h"


/** Represents a part of a read. */
typedef struct _READ_PART {
	/** Part position within the reference sequence. */
	uint64_t Position;
	/** The read part sequence (it is not null-terminated). */
	char *ReadSequence;
	/** Length of the part sequence, in bases. */
	size_t ReadSequenceLength;
	uint8_t *Quality;
	/** Position of the part relative to the start of the read. */
	size_t Offset;
} READ_PART, *PREAD_PART;

GEN_ARRAY_TYPEDEF(READ_PART);
GEN_ARRAY_IMPLEMENTATION(READ_PART)

#define READ_TOP_QUALITY						(255 - 33)

struct _ONE_READ;

typedef struct _ONE_READ_PAIRED_LIST {
	struct _ONE_READ *Next;
} ONE_READ_PAIRED_LIST, *PONE_READ_PAIRED_LIST;

typedef struct _ONE_READ_EXTENSION {
	uint16_t Flags;
	char *CIGAR;
	char *RName;
	char *RNext;
	int32_t TLen;
	uint64_t PNext;
	char *TemplateName;
} ONE_READ_EXTENSION, *PONE_READ_EXTENSION;

typedef struct _ONE_READ {
	char *ReadSequence;
	uint32_t ReadSequenceLen;
	uint32_t RealReadSequenceLen;
	uint8_t *Quality;
	uint64_t Pos;
	uint8_t PosQuality;
	uint32_t NumberOfFixes;
	size_t ReadIndex;
	struct _ONE_READ *Parent;
	uint32_t Offset;
	boolean NoStartStrip;
	PONE_READ_EXTENSION Extension;
	boolean NoEndStrip;
} ONE_READ, *PONE_READ;

typedef struct _ASSEMBLY_TASK {
	boolean Allocated;
	const char *Name;
	const char *Reference;
	size_t ReferenceLength;
	const char *Alternate1;
	size_t Alternate1Length;
	const char *Alternate2;
	size_t Alternate2Length;
	const ONE_READ *Reads;
	size_t ReadCount;
	uint64_t RegionStart;
} ASSEMBLY_TASK, *PASSEMBLY_TASK;

GEN_ARRAY_TYPEDEF(ONE_READ);
GEN_ARRAY_IMPLEMENTATION(ONE_READ)
POINTER_ARRAY_TYPEDEF(ONE_READ);
POINTER_ARRAY_IMPLEMENTATION(ONE_READ)



void read_quality_decode(PONE_READ Read);
void read_quality_encode(PONE_READ Read);

void read_write_fastq(FILE *Stream, const ONE_READ *Read);
void read_write_sam(FILE *Stream, const ONE_READ *Read);
ERR_VALUE read_create_from_sam_line(const char *Line, PONE_READ Read);
ERR_VALUE read_create_from_fastq(const char *Block, const char **NewBlock, PONE_READ Read);

void read_destroy(PONE_READ Read);
void _read_destroy_structure(PONE_READ Read);
ERR_VALUE read_concat(PONE_READ Target, const ONE_READ *Source);

void read_set_destroy(PONE_READ ReadSet, const size_t Count);
ERR_VALUE read_set_merge(PONE_READ *Target, const size_t TargetCount, struct _ONE_READ *Source, const size_t SourceCount);
void read_split(PONE_READ Read);
void read_adjust(PONE_READ Read, const uint64_t RegionStart, const size_t RegionLength);
void read_shorten(PONE_READ Read, const size_t Count);
ERR_VALUE read_base_insert(PONE_READ Read, const char Base, size_t Index);
void read_base_delete(PONE_READ Read, size_t Index);

void assembly_task_init(PASSEMBLY_TASK Task, const char *RefSeq, const size_t RefSeqLen, const char *Alternate1, const size_t Alternate1Length, const char *Alternate2, const size_t Alternate2Length, const ONE_READ *ReadSet, const size_t ReadCount);
void assembly_task_set_name(PASSEMBLY_TASK Task, const char *Name);
void assembly_task_finit(PASSEMBLY_TASK Task);



#endif 
