
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "gassm.h"
#include "options.h"
#include "input-file.h"


/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/


static UTILS_TYPED_CALLOC_FUNCTION(ACTIVE_REGION)

static boolean _fasta_read_seq_raw(char *Start, size_t Length, char **SeqStart, char **SeqEnd)
{
	boolean ret = FALSE;

	ret = *Start == '>';
	if (ret) {
		while (Length > 0 && *Start != '\n' && *Start != '\r' && *Start != 26) {
			++Start;
			--Length;
		}

		if (Length > 0) {
			while (Length > 0 && (*Start == '\n' || *Start == '\r')) {
				++Start;
				--Length;
			}

			if (Length > 0 && *Start != 26) {
				*SeqStart = Start;
				while (Length > 0) {
					if (*Start == '>' && (*(Start - 1) == '\n' || *(Start - 1) == '\r'))
						break;

					++Start;
					--Length;
				}

				*SeqEnd = Start;
			} else {
				*SeqStart = NULL;
				*SeqEnd = NULL;
			}
		} else { 
			*SeqStart = NULL;
			*SeqEnd = NULL;
		}
	}

	return ret;
}

static ERR_VALUE _fasta_read_seq(char *Start, size_t Length, char **NewStart, char **Seq, size_t *SeqLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *seqStart = NULL;
	char *seqEnd = NULL;
	char *tmpSeq = NULL;
	size_t tmpSeqLen = 0;

	if (_fasta_read_seq_raw(Start, Length, &seqStart, &seqEnd)) {
		*Seq = NULL;
		*SeqLen = 0;
		if (seqStart != NULL) {
			ret = utils_malloc(seqEnd - seqStart + sizeof(char), &tmpSeq);
			if (ret == ERR_SUCCESS) {
				char *tmp = tmpSeq;

				while (seqStart != seqEnd) {
					if (*seqStart != '\n' && *seqStart != '\r' && *seqStart != 26) {
						*tmp = *seqStart;
						++tmp;
						++tmpSeqLen;
					}

					++seqStart;
				}

				tmpSeq[tmpSeqLen] = '\0';
				*NewStart = seqEnd;
				*Seq = tmpSeq;
				*SeqLen = tmpSeqLen;
			}
		} else ret = ERR_SUCCESS;
	}

	return ret;
}

/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE input_get_refseq(char **RefSeq, size_t *RefSeqLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *fileName = NULL;
	char *type = NULL;

	ret = option_get_String(GASSM_OPTION_REFSEQ_INPUT_FILE, &fileName);
	if (ret == ERR_SUCCESS) {
		ret = option_get_String(GASSM_OPTION_REFSEQ_INPUT_TYPE, &type);
		if (ret == ERR_SUCCESS) {
			if (strcasecmp(type, "fasta") == 0) {
				char *data = NULL;
				char *dummy = NULL;
				size_t dataLength = 0;

				ret = utils_file_read(fileName, &data, &dataLength);
				if (ret == ERR_SUCCESS) {
					ret = _fasta_read_seq(data, dataLength, &dummy, RefSeq, RefSeqLen);
					utils_free(data);
				}
			} else ret = ERR_UNKNOWN_REFSEQ_INPUT_TYPE;
		}
	}

	return ret;
}


void input_free_refseq(char *RefSeq, const size_t RefSeqLen)
{
	if (RefSeqLen > 0)
		utils_free(RefSeq);

	return;
}


ERR_VALUE input_get_reads(char ***Reads, size_t *ReadCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	char *fileName = NULL;
	char *type = NULL;

	ret = option_get_String(GASSM_OPTION_READS_INPUT_FILE, &fileName);
	if (ret == ERR_SUCCESS) {
		ret = option_get_String(GASSM_OPTION_READS_INPUT_TYPE, &type);
		if (ret == ERR_SUCCESS) {
			if (strcasecmp(type, "fasta") == 0) {
				char *data = NULL;
				size_t dataLength = 0;

				ret = utils_file_read(fileName, &data, &dataLength);
				if (ret == ERR_SUCCESS) {
					char *start = data;
					char *read = NULL;
					size_t readLen = 0;
					char **tmpReads = NULL;
					size_t tmpReadsCount = 0;

					*Reads = NULL;
					*ReadCount = 0;
					while (ret == ERR_SUCCESS && start != NULL && dataLength > 0) {
						char *oldPos = start;

						ret = _fasta_read_seq(start, dataLength, &start, &read, &readLen);
						if (ret == ERR_SUCCESS) {
							if (read != NULL) {
								char **oldReads = tmpReads;

								dataLength -= (start - oldPos);
								ret = utils_calloc(tmpReadsCount + 1, sizeof(char *), (void **)&tmpReads);
								if (ret == ERR_SUCCESS) {
									memcpy(tmpReads, oldReads, tmpReadsCount*sizeof(char *));
									tmpReads[tmpReadsCount] = read;
									++tmpReadsCount;
									if (oldReads != NULL)
										utils_free(oldReads);
								}
							} else dataLength = 0;

							if (ret != ERR_SUCCESS)
								utils_free(read);
						}

						if (ret != ERR_SUCCESS) {
							for (size_t i = 0; i < tmpReadsCount; ++i)
								utils_free(tmpReads[i]);

							utils_free(tmpReads);
						}
					}

					if (ret == ERR_SUCCESS) {
						*Reads = tmpReads;
						*ReadCount = tmpReadsCount;
					}

					utils_free(data);
				}
			}
			else if (strcasecmp(type, "none") == 0) {
				*Reads = NULL;
				*ReadCount = 0;
				ret = ERR_SUCCESS;
			} else ret = ERR_UNKNOWN_READS_INPUT_TYPE;
		}
	}

	return ret;
}


void input_free_reads(char **Reads, const size_t Count)
{
	if (Count > 0) {
		for (size_t i = 0; i < Count; ++i)
			utils_free(Reads[i]);

		utils_free(Reads);
	}

	return;
}


ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count)
{
	const char *regStart = NULL;
	const char *regEnd = NULL;
	size_t tmpArrayLen = 0;
	PACTIVE_REGION tmpArray = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	*Regions = NULL;
	*Count = 0;
	ret = ERR_SUCCESS;
	regStart = RefSeq;
	regEnd = RefSeq;
	for (size_t i = 0; i < RefSeqLen; ++i) {
		switch (*regEnd) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
				++regEnd;
				break;
			case 'N':
				if (regEnd != regStart) {
					PACTIVE_REGION old = tmpArray;

					ret = utils_calloc_ACTIVE_REGION(tmpArrayLen + 1, &tmpArray);
					if (ret == ERR_SUCCESS) {
						PACTIVE_REGION newReg = tmpArray + tmpArrayLen;

						memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
						newReg->Sequence = regStart;
						newReg->Offset = (regStart - RefSeq);
						newReg->Length = (regEnd - regStart);
						++tmpArrayLen;
						if (old != NULL)
							utils_free(old);
					}
				}

				++regEnd;
				regStart = regEnd;
				break;
			default:
				printf("Unknown character in the reference sequence: %u\n", *regEnd);
				break;
		}

		if (ret != ERR_SUCCESS) {
			if (tmpArray != NULL)
				utils_free(tmpArray);

			break;
		}
	}

	if (ret == ERR_SUCCESS) {
		if (regEnd != regStart) {
			PACTIVE_REGION old = tmpArray;

			ret = utils_calloc_ACTIVE_REGION(tmpArrayLen + 1, &tmpArray);
			if (ret == ERR_SUCCESS) {
				PACTIVE_REGION newReg = tmpArray + tmpArrayLen;
				
				memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
				newReg->Sequence = regStart;
				newReg->Offset = (regStart - RefSeq);
				newReg->Length = (regEnd - regStart);
				tmpArrayLen++;
				if (old != NULL)
					utils_free(old);
			}
		}

		if (ret == ERR_SUCCESS) {
			*Regions = tmpArray;
			*Count = tmpArrayLen;
		}

		if (ret != ERR_SUCCESS) {
			if (tmpArray != NULL)
				utils_free(tmpArray);
		}
	}

	return ret;
}


void input_free_regions(PACTIVE_REGION Regions, const size_t Count)
{
	if (Count > 0)
		utils_free(Regions);

	return;
}
