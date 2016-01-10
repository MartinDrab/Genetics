
#ifndef __SSW_H__
#define __SSW_H__


#include "err.h"



ERR_VALUE ssw_simple(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);
ERR_VALUE ssw_clever(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);
ERR_VALUE ssw_sse(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);
ERR_VALUE write_differences(const char *RefSeq, const size_t RegionStart, const size_t RegionLength, const char *OpString, const char *AltSeq);



#endif 
