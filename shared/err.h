
#ifndef __ERR_H__
#define __ERR_H__


typedef int ERR_VALUE;

#define ERR_SUCCESS								0
#define ERR_INTERNAL_ERROR						1
#define ERR_OUT_OF_MEMORY						2
#define ERR_TYPE_MISMATCH						3
#define ERR_NOT_FOUND							4
#define ERR_TABLE_FULL							5
#define ERR_NOT_A_PRIME							6
#define ERR_UNKNOWN_OPTION						7
#define ERR_NOT_AN_INTEGRAL_NUMBER				8
#define ERR_NOT_AN_UNSIGNED_INTEGRAL_NUMBER		9
#define ERR_NOT_A_FLOAT_NUMBER					10
#define ERR_NOT_A_DOUBLE_NUMBER					11
#define ERR_EMPTY_STRING_NOT_ALLOWED			12
#define ERR_STRING_TOO_LONG						13
#define ERR_INVALID_BOOLEAN_VALUE				14
#define ERR_OPTION_VALUE_NOT_FOUND				15
#define ERR_NO_INVERSE							16
#define ERR_MODULUS_TOO_SMALL					17
#define ERR_ALREADY_EXISTS						18
#define ERR_ERRNO_VALUE							19
#define ERR_FERROR								20
#define ERR_UNKNOWN_REFSEQ_INPUT_TYPE			21
#define ERR_UNKNOWN_READS_INPUT_TYPE			22
#define ERR_NO_MORE_ENTRIES						23
#define ERR_OFFSET_TOO_LOW						24
#define ERR_OFFSET_TOO_HIGH						25
#define ERR_NOT_IMPLEMENTED						26

#define ERR_SAM_INVALID_QNAME					27
#define ERR_SAM_INVALID_FLAG					28
#define ERR_SAM_INVALID_RNAME					29
#define ERR_SAM_INVALID_POS 					30
#define ERR_SAM_INVALID_MAPQ					31
#define ERR_SAM_INVALID_CIGAR					32
#define ERR_SAM_INVALID_RNEXT					33
#define ERR_SAM_INVALID_PNEXT					34
#define ERR_SAM_INVALID_TLEN					35
#define ERR_SAM_INVALID_SEQ 					36
#define ERR_SAM_INVALID_QUAL					37
#define ERR_SAM_SEQ_QUAL_LEN_MISMATCH			38


#endif 
