
#include <stdio.h>
#include <malloc.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "kmer.h"
#include "kmer-table.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 3);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_TEST_INPUT_FILE, "input.txt");

	return ret;
}

static ERR_VALUE _parse_test_data(char *Data, size_t DataLength, char **RefSeq, char ***Reads, size_t *ReadsCount)
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
				} else ret = ERR_OUT_OF_MEMORY;
			} else ret = ERR_SUCCESS;

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


void _free_test_data(char *RefSeq, char **Reads, const size_t ReadsCount)
{
	for (size_t i = 0; i < ReadsCount; ++i)
		free(Reads[i]);

	free(Reads);
	free(RefSeq);

	return;
}

int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _set_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				uint32_t kmerSize = 0;
				char *testFile = NULL;

				ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
				if (ret == ERR_SUCCESS)
					ret = option_get_String(GASSM_OPTION_TEST_INPUT_FILE, &testFile);
				
				if (ret == ERR_SUCCESS) {
					PKMER_TABLE kmerTable = NULL;

					ret = kmer_table_create(kmerSize, 2, 37, &kmerTable);
					if (ret == ERR_SUCCESS) {
						char *testData = NULL;
						size_t testDataSize = 0;

						ret = utils_file_read(testFile, &testData, &testDataSize);
						if (ret == ERR_SUCCESS) {
							char *refSeq = NULL;
							char **reads = NULL;
							size_t readCount = 0;

							ret = _parse_test_data(testData, testDataSize, &refSeq, &reads, &readCount);
							if (ret == ERR_SUCCESS) {
								printf("Ref. sequence: %s\n", refSeq);
								printf("Number of reads: %u\n", readCount);
								for (size_t i = 0; i < readCount; ++i)
									printf("Read #%u: %s\n", i, reads[i]);

								_free_test_data(refSeq, reads, readCount);
							}
							free(testData);
						}

						kmer_table_destroy(kmerTable);
					}
				}
			}
		}

		options_module_finit();
	}

	return ret;
}
