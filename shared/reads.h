
#ifndef __READS_H__
#define __READS_H__

#include "err.h"
#include "utils.h"



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


ERR_VALUE read_create_from_test_line(const char *Line, const size_t Length, PONE_READ *Read);
ERR_VALUE read_create_from_sam_line(const char *Line, PONE_READ *Read);
ERR_VALUE read_create_from_fasta_seq(const char *Seq, const size_t SeqLen, const char *SeqName, const size_t SeqNameLen, PONE_READ *Read);
void read_destroy(PONE_READ Read);






#endif 
