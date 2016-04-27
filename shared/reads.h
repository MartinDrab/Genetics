
#ifndef __READS_H__
#define __READS_H__

#include "err.h"
#include "utils.h"
#include "gen_dym_array.h"


typedef struct _ONE_READ {
	uint16_t Flags;
	char *ReadSequence;
	size_t ReadSequenceLen;
	uint8_t *Quality;
	size_t QualityLen;
	char *CIGAR;
	size_t CIGARLen;
	uint64_t Pos;
	uint8_t PosQuality;
	char *TemplateName;
	size_t TemplateNameLen;
} ONE_READ, *PONE_READ;

typedef struct _ASSEMBLY_TASK {
	boolean Allocated;
	const char *Reference;
	size_t ReferenceLength;
	const char *Alternate1;
	size_t Alternate1Length;
	const char *Alternate2;
	size_t Alternate2Length;
	const ONE_READ *Reads;
	size_t ReadCount;
} ASSEMBLY_TASK, *PASSEMBLY_TASK;

GEN_ARRAY_TYPEDEF(ONE_READ);
GEN_ARRAY_IMPLEMENTATION(ONE_READ)
GEN_ARRAY_TYPEDEF(PONE_READ);
GEN_ARRAY_IMPLEMENTATION(PONE_READ)


ERR_VALUE read_create_from_test_line(const char *Line, const size_t Length, PONE_READ *Read);
ERR_VALUE read_create_from_sam_line(const char *Line, PONE_READ *Read);
ERR_VALUE read_create_from_fasta_seq(const char *Seq, const size_t SeqLen, const char *SeqName, const size_t SeqNameLen, PONE_READ *Read);
ERR_VALUE read_generate_from_sequence(const char *Seq, const size_t SeqLen, const uint32_t ReadLength, PONE_READ *Read);
void read_destroy(PONE_READ Read);

ERR_VALUE read_set_generate_from_sequence(const char *Seq, const size_t SeqLen, const uint32_t ReadLength, const size_t ReadCount, PONE_READ *ReadSet);
void read_set_destroy(PONE_READ ReadSet, const size_t Count);
ERR_VALUE read_set_merge(PONE_READ *Target, const size_t TargetCount, struct _ONE_READ *Source, const size_t SourceCount);

ERR_VALUE read_save(FILE *Stream, const ONE_READ *Read);
ERR_VALUE read_load(FILE *Stream, PONE_READ Read);
ERR_VALUE read_set_save(FILE *Stream, const ONE_READ *ReadSet, const size_t Count);
ERR_VALUE read_set_load(FILE *Stream, PONE_READ *ReadSet, size_t *Count);
ERR_VALUE seq_save(FILE *Stream, const char *RefSeq, const size_t Length);
ERR_VALUE seq_load(FILE *Stream, char **RefSeq, size_t *Length);

void assembly_task_init(PASSEMBLY_TASK Task, const char *RefSeq, const size_t RefSeqLen, const char *Alternate1, const size_t Alternate1Length, const char *Alternate2, const size_t Alternate2Length, const ONE_READ *ReadSet, const size_t ReadCount);
void assembly_task_finit(PASSEMBLY_TASK Task);
ERR_VALUE assembly_task_save(FILE *Stream, const ASSEMBLY_TASK *Task);
ERR_VALUE assembly_task_save_file(const char *FileName, const ASSEMBLY_TASK *Task);
ERR_VALUE assembly_task_load(FILE *Stream, PASSEMBLY_TASK Task);
ERR_VALUE assembly_task_load_file(const char *FileName, PASSEMBLY_TASK Task);



#endif 
