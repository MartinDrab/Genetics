
#ifndef __SSW_H__
#define __SSW_H__


#include "err.h"



typedef enum _EGapType {
	gtMatch,
	gtMismatch,
	gtDeletion,
	gtInsertion,
} EGapType, *PEGapType;


void ssw_clever(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);



#endif 
