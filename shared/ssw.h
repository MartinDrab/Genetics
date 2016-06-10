
#ifndef __SSW_H__
#define __SSW_H__


#include "err.h"
#include "found-sequence.h"


typedef struct _SSW_STATISTICS {
	int TotalInsertions;
	int TotalDeletions;
	int TotalMismatches;
	int InsertGapCount;
	int DeleteGapCount;
	int MismatchGapCount;
} SSW_STATISTICS, *PSSW_STATISTICS;


typedef enum _EGapType {
	gtMatch,
	gtMismatch,
	gtDeletion,
	gtInsertion,
} EGapType, *PEGapType;


ERR_VALUE ssw_simple(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);
ERR_VALUE ssw_clever(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);
void opstring_statistics(const char *OpString, const size_t OpStringLen, SSW_STATISTICS *Statistics);
ERR_VALUE write_seq_differences(PGEN_ARRAY_VARIANT_CALL VCArray, const char *RefSeq, const size_t RegionStart, const size_t RegionLength, const char *OpString, const char *AltSeq, const FOUND_SEQUENCE_VARIANT *Variant);



#endif 
