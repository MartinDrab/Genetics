
#include <stdlib.h>
#include <assert.h>
#include "err.h"
#include "utils.h"



/************************************************************************/
/*                                                                      */
/************************************************************************/

ERR_VALUE parse_test_data(char *Data, size_t DataLength, char **RefSeq, char ***Reads, size_t *ReadsCount)
{
	char *tmpRefSeq = NULL;
	size_t refSeqLen = 0;

	char *tmp = Data;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	// 1) Capture the reference sequence
	while (DataLength && (*tmp != '\r' && *tmp != '\n' && *tmp != 26)) {
		++tmp;
		--DataLength;
	}

	refSeqLen = tmp - Data;
	if (DataLength > 0 && refSeqLen > 0) {
		ret = utils_preallocate_string(refSeqLen, &tmpRefSeq);
		if (ret == ERR_SUCCESS) {
			memcpy(tmpRefSeq, Data, refSeqLen*sizeof(char));
			// 2) Determine number of reads
			size_t tmpReadsCount = 0;
			char **tmpReadArray = NULL;

			while (DataLength > 0 && (*tmp == '\n' || *tmp == '\r')) {
				++tmp;
				--DataLength;
			}

			if (*tmp != 26 && DataLength > 0) {
				char *readStart = tmp;

				while (DataLength > 0) {
					while (DataLength > 0 && (*tmp != '\r' && *tmp != '\n' && *tmp != 26)) {
						++tmp;
						--DataLength;
					}

					++tmpReadsCount;
					while (DataLength > 0 && (*tmp == '\n' || *tmp == '\r')) {
						++tmp;
						--DataLength;
					}
				}

				// 3) Copy the reads
				tmpReadArray = (char **)calloc(tmpReadsCount, sizeof(char *));
				if (tmpReadArray != NULL) {
					size_t index = 0;

					while (ret == ERR_SUCCESS && index < tmpReadsCount) {
						tmp = readStart;
						while ((*tmp != '\r' && *tmp != '\n' && *tmp != 26))
							++tmp;

						char *read = NULL;
						size_t readLen = tmp - readStart;
						ret = utils_preallocate_string(readLen, &read);
						if (ret == ERR_SUCCESS) {
							memcpy(read, readStart, readLen*sizeof(char));
							while ((*tmp == '\n' || *tmp == '\r'))
								++tmp;

							tmpReadArray[index] = read;
							++index;
							readStart = tmp;
						}
					}

					if (ret != ERR_SUCCESS) {
						for (size_t i = 0; i < index; ++i)
							free(tmpReadArray[i]);

						free(tmpReadArray);
					}
				}
				else ret = ERR_OUT_OF_MEMORY;
			}
			else ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				*ReadsCount = tmpReadsCount;
				*Reads = tmpReadArray;
				*RefSeq = tmpRefSeq;
			}

			if (ret != ERR_SUCCESS)
				utils_free_string(tmpRefSeq);
		}
	}

	return ret;
}


void free_test_data(char *RefSeq, char **Reads, const size_t ReadsCount)
{
	for (size_t i = 0; i < ReadsCount; ++i)
		free(Reads[i]);

	free(Reads);
	free(RefSeq);

	return;
}
