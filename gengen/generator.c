
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "err.h"
#include "utils.h"
#include "generator.h"


/************************************************************************/
/*                         PUBLIC FUNCTIONS                             */
/************************************************************************/

ERR_VALUE generate_active_region(uint32_t MinLength, uint32_t MaxLength, char *BaseTypes, char **Result, size_t *ResultLength)
{
	size_t tmpLength = 0;
	char *tmpResult = NULL;
	uint32_t baseTypeCount = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	baseTypeCount =(uint32_t)strlen(BaseTypes);
	tmpLength = utils_ranged_rand(MinLength, MaxLength + 1);
	ret = utils_preallocate_string(tmpLength, &tmpResult);
	if (ret == ERR_SUCCESS) {
		char *currentBase = tmpResult;
		
		for (size_t i = 0; i < tmpLength; ++i) {
			*currentBase = BaseTypes[utils_ranged_rand(0, baseTypeCount)];
			++currentBase;
		}

		tmpResult[tmpLength / sizeof(char)] = '\0';
		*Result = tmpResult;
		*ResultLength = tmpLength;
	}

	return ret;
}

void free_active_region(char *Region, const size_t Length)
{
	free(Region);

	return;
}


ERR_VALUE generate_reads(const char *Region, const size_t RegionLength, const char *Nucleotides, const PGENERATOR_READ_OPTIONS ReadOptions, const PGENERATOR_INDEL_OPTIONS IndelOptions, const PGENERATOR_REPLACE_OPTIONS ReplaceOptions, char ***Reads, size_t *ReadCount)
{
	char **tmpArray = NULL;

	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	uint32_t readLength = 0;
	size_t readCount = 0;
	size_t tmpReadStorageSize = 0;
	char *tmpReadStorage = NULL;
	uint32_t nucleotidesCount = 0;

	readCount = utils_ranged_rand(ReadOptions->MinReads, ReadOptions->MaxReads + 1);
	tmpArray = (char **)calloc(readCount, sizeof(char *));
	if (tmpArray != NULL) {
		nucleotidesCount = (uint32_t)strlen(Nucleotides);
		readLength = ReadOptions->ReadLength;
		tmpReadStorageSize = readLength + IndelOptions->MaxIndels;
		ret = utils_preallocate_string(tmpReadStorageSize, &tmpReadStorage);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < readCount; ++i) {
				memset(tmpReadStorage, 0, readLength + IndelOptions->MaxIndels);
				// 1) Randomly pick the read from the active region and copy it to the temporaty storage
				size_t readStart = utils_ranged_rand(0, RegionLength - readLength);
				tmpReadStorageSize = readLength;
				memcpy(tmpReadStorage, Region + readStart, tmpReadStorageSize*sizeof(char));
				// 2) Determine if any indels will be applied
				if (utils_prob_happened(IndelOptions->IndelProbability)) {
					// 2.a) Establish number of the indel operations
					size_t indelCount = utils_ranged_rand(IndelOptions->MinIndels, IndelOptions->MaxIndels + 1);
					// 2.b) Apply them one by one
					for (size_t j = 0; j < indelCount; ++j) {
						size_t place = utils_ranged_rand(0, tmpReadStorageSize);
						if (!IndelOptions->DisableIns && utils_prob_happened(0.5)) {
							// An insertion operation
							memmove(tmpReadStorage + place + 1, tmpReadStorage + place, (tmpReadStorageSize - place)*sizeof(char));
							char newBase = '\0';
							do {
								newBase = Nucleotides[utils_ranged_rand(0, nucleotidesCount)];
							} while (newBase == tmpReadStorage[place]);

							tmpReadStorage[place] = newBase;
							++tmpReadStorageSize;
						}
						else if (!IndelOptions->DisableDels) {
							// A deletion operation
							memmove(tmpReadStorage + place, tmpReadStorage + place + 1, (tmpReadStorageSize - place - 1)*sizeof(char));
							tmpReadStorage[tmpReadStorageSize - 1] = '\0';
							--tmpReadStorageSize;
						}
					}
				}

				// 3) Determine if any replaces will occur.
				if (utils_prob_happened(ReplaceOptions->ReplaceProbability)) {
					// 3.1) Determine the number of replaces
					size_t replaceCount = utils_ranged_rand(ReplaceOptions->MinReplace, ReplaceOptions->MaxReplace + 1);
					for (size_t j = 0; j < replaceCount; ++j) {
						size_t place = utils_ranged_rand(0, tmpReadStorageSize);
						char newBase = '\0';
						do {
							newBase = Nucleotides[utils_ranged_rand(0, nucleotidesCount)];
						} while (newBase == tmpReadStorage[place]);

						tmpReadStorage[place] = newBase;
					}
				}

				// 4) Store the resulting read into the output array
				ret = utils_copy_string(tmpReadStorage, &tmpArray[i]);			
				if (ret != ERR_SUCCESS) {
					for (uint32_t j = 0; j < i; ++j)
						utils_free_string(tmpArray[j]);

					break;
				}
			}

			utils_free_string(tmpReadStorage);
		}
	
		if (ret == ERR_SUCCESS) {
			*Reads = tmpArray;
			*ReadCount = readCount;
		}

		if (ret != ERR_SUCCESS)
			free(tmpArray);
	} else ret = ERR_OUT_OF_MEMORY;

	return ret;
}

void free_reads(char **Reads, const size_t Count)
{
	for (size_t i = 0; i < Count; ++i)
		free(Reads[i]);

	free(Reads);

	return;
}
