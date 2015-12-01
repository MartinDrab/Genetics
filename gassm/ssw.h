
#ifndef __SSW_H__
#define __SSW_H__


#include "err.h"



ERR_VALUE ssw_simple(const char *A, const size_t ALen, const char *B, const size_t BLen, const int Match, const int Mismatch, const int Indel, char **OperationString, size_t *OperationStringLen);


#endif 
