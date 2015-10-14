
#ifndef __TEST_FILE_INPUT_H__
#define __TEST_FILE_INPUT_H__


#include "err.h"
#include "utils.h"



ERR_VALUE parse_test_data(char *Data, size_t DataLength, char **RefSeq, char ***Reads, size_t *ReadsCount);
void free_test_data(char *RefSeq, char **Reads, const size_t ReadsCount);



#endif 
