
#ifndef __SSW_H__
#define __SSW_H__


#include "err.h"
#include "variant.h"


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



#endif 
