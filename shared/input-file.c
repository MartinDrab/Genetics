
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "dym-array.h"
#include "gen_dym_array.h"
#include "reads.h"
#include "input-file.h"


/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/


static const char *_read_line(const char *LineStart)
{
	while (*LineStart != '\n' && *LineStart != '\r' && *LineStart != 26 && *LineStart != '\0')
		++LineStart;

	return LineStart;
}


static const char *_advance_to_next_line(const char *LineEnd)
{
	while (*LineEnd == '\n' || *LineEnd == '\r')
		++LineEnd;

	return LineEnd;
}


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
		} else ret = ERR_NO_MORE_ENTRIES;
	} else ret = ERR_NO_MORE_ENTRIES;

	return ret;
}

/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE fasta_load(const char *FileName, PFASTA_FILE FastaRecord)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_file_read(FileName, &FastaRecord->FileData, &FastaRecord->DataLength);
	if (ret == ERR_SUCCESS)
		FastaRecord->CurrentPointer = FastaRecord->FileData;

	return ret;
}



ERR_VALUE fasta_read_seq(PFASTA_FILE FastaRecord, char **Seq, size_t *SeqLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = _fasta_read_seq(FastaRecord->CurrentPointer, FastaRecord->DataLength - (FastaRecord->CurrentPointer - FastaRecord->FileData), &FastaRecord->CurrentPointer, Seq, SeqLen);

	return ret;
}


void fasta_free(PFASTA_FILE FastaRecord)
{
	if (FastaRecord->DataLength > 0)
		utils_free(FastaRecord->FileData);

	return;
}



ERR_VALUE input_get_refseq(const char *FileName, const char *InputType, char **RefSeq, size_t *RefSeqLen)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (strcasecmp(InputType, "fasta") == 0) {
		FASTA_FILE f;

		ret = fasta_load(FileName, &f);
		if (ret == ERR_SUCCESS) {
			ret = fasta_read_seq(&f, RefSeq, RefSeqLen);
			fasta_free(&f);
		}
	} else if (strcasecmp(InputType, "test") == 0) {
		char *data = NULL;
		char *dummy = NULL;
		size_t dataLength = 0;

		ret = utils_file_read(FileName, &data, &dataLength);
		if (ret == ERR_SUCCESS) {
			const char *end = NULL;
			size_t len = 0;

			end = _read_line(data);
			len = (end - data);
			if (len > 0) {
				ret = utils_calloc(len + 1, sizeof(char), &dummy);
				if (ret == ERR_SUCCESS) {
					memcpy(dummy, data, len*sizeof(char));
					dummy[len] = '\0';
					*RefSeq = dummy;
					*RefSeqLen = len;
				}
			} else ret = ERR_INTERNAL_ERROR;

			utils_free(data);
		}
	} else ret = ERR_UNKNOWN_REFSEQ_INPUT_TYPE;

	return ret;
}


void input_free_refseq(char *RefSeq, const size_t RefSeqLen)
{
	if (RefSeqLen > 0)
		utils_free(RefSeq);

	return;
}


ERR_VALUE input_get_reads(const char *Filename, const char *InputType, const uint64_t RegionStart, const uint64_t RegionSize, PONE_READ *Reads, size_t *ReadCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	/*
	if (strcasecmp(InputType, "fasta") == 0) {
		char *data = NULL;
		size_t dataLength = 0;

		ret = utils_file_read(Filename, &data, &dataLength);
		if (ret == ERR_SUCCESS) {
			ret = ERR_NOT_IMPLEMENTED;
			utils_free(data);
		}
	} else */
	if (strcasecmp(InputType, "sam") == 0) {
		char *data = NULL;
		size_t dataLength = 0;

		ret = utils_file_read(Filename, &data, &dataLength);
		if (ret == ERR_SUCCESS) {
			const char *line = data;
			PONE_READ oneRead = NULL;
			const char *lineEnd = _read_line(line);
			DYM_ARRAY readArray;

			dym_array_create(&readArray, 140);
			while (ret == ERR_SUCCESS && line != lineEnd) {
				ret = read_create_from_sam_line(line, &oneRead);
				if (ret == ERR_SUCCESS) {
					if (
						(RegionStart <= oneRead->Pos && oneRead->Pos < RegionStart + RegionSize) ||
						(RegionStart <= oneRead->Pos + oneRead->ReadSequenceLen && oneRead->Pos + oneRead->ReadSequenceLen < RegionStart + RegionSize)) {
						if (oneRead->ReadSequenceLen != (uint64_t)-1) {
							// TODO: Well, let's cut  the read so it does not go outside the active region
						}

						ret = dym_array_push_back(&readArray, oneRead);
					} else ret = ERR_NOT_IN_REGION;

					if (ret != ERR_SUCCESS)
						read_destroy(oneRead);
				
					if (ret == ERR_NOT_IN_REGION)
						ret = ERR_SUCCESS;
				}

				line = _advance_to_next_line(lineEnd);
				lineEnd = _read_line(line);
			}

			if (ret == ERR_SUCCESS) {
				PONE_READ tmpReads = NULL;
				size_t tmpReadCount = dym_array_size(&readArray);

				ret = utils_calloc(tmpReadCount, sizeof(ONE_READ), &tmpReads);
				if (ret == ERR_SUCCESS) {
					PONE_READ tmp = tmpReads;
					
					for (size_t i = 0; i < tmpReadCount; ++i) {
						memcpy(tmp, dym_array_get(&readArray, i), sizeof(ONE_READ));
						utils_free(dym_array_get(&readArray, i));
						++tmp;
					}
					
					*Reads = tmpReads;
					*ReadCount = dym_array_size(&readArray);
				}
			}

			if (ret != ERR_SUCCESS) {
				size_t len = dym_array_size(&readArray);

				for (size_t i = 0; i < len; ++i)
					read_destroy((PONE_READ)dym_array_get(&readArray, i));
			}

			dym_array_destroy(&readArray);
			utils_free(data);
		}
	}
	/*
	else if (strcasecmp(InputType, "none") == 0) {
		*Reads = NULL;
		*ReadCount = 0;
		ret = ERR_SUCCESS;
	} else if (strcasecmp(InputType, "test") == 0) {
		char *data = NULL;
		size_t dataLength = 0;

		ret = utils_file_read(Filename, &data, &dataLength);
		if (ret == ERR_SUCCESS) {
			const char *start = data;
			const char *end = _read_line(start);
			DYM_ARRAY ra;

			dym_array_create(&ra, 140);
			while (ret == ERR_SUCCESS && start != end) {
				PONE_READ r = NULL;

				ret = read_create_from_test_line(start, end - start, &r);
				if (ret == ERR_SUCCESS) {
					ret = dym_array_push_back(&ra, r);
					if (ret != ERR_SUCCESS)
						read_destroy(r);
				}
				
				start = _advance_to_next_line(end);
				end = _read_line(start);
			}

			if (ret == ERR_SUCCESS) {
				ret = dym_array_to_array(&ra, Reads);
				if (ret == ERR_SUCCESS)
					*ReadCount = dym_array_size(&ra);
			}

			if (ret != ERR_SUCCESS) {
				for (size_t i = 0; i < dym_array_size(&ra); ++i)
					read_destroy((PONE_READ)dym_array_data(&ra)[i]);
			}

			dym_array_destroy(&ra);
			utils_free(data);
		}
	}*/ 
	else ret = ERR_UNKNOWN_READS_INPUT_TYPE;
	
	return ret;
}


void input_free_reads(PONE_READ Reads, const size_t Count)
{
	read_set_destroy(Reads, Count);

	return;
}


ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count)
{
	const char *regStart = NULL;
	const char *regEnd = NULL;
	size_t tmpArrayLen = 0;
	uint64_t numberOfNBases = 0;
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
				if (numberOfNBases > 0) {
					PACTIVE_REGION old = tmpArray;

					ret = utils_calloc(tmpArrayLen + 1, sizeof(ACTIVE_REGION), &tmpArray);
					if (ret == ERR_SUCCESS) {
						PACTIVE_REGION newReg = tmpArray + tmpArrayLen;

						memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
						newReg->Sequence = regStart;
						newReg->Offset = (newReg->Sequence - RefSeq);
						newReg->Length = (regEnd - regStart);
						newReg->Type = artUnknown;
						++tmpArrayLen;
						if (old != NULL)
							utils_free(old);
					}

					numberOfNBases = 0;
					regStart = regEnd;
				}

				++regEnd;
				break;
			case 'N':
			case 'M':
			case 'R':
			case 'Y':
			case 'W':
			case 'S':
			case 'K':
			case 'V':
			case 'H':
			case 'D':
			case 'B':
				if (numberOfNBases == 0 && regEnd != regStart) {
					PACTIVE_REGION old = tmpArray;

					ret = utils_calloc(tmpArrayLen + 1, sizeof(ACTIVE_REGION), &tmpArray);
					if (ret == ERR_SUCCESS) {
						PACTIVE_REGION newReg = tmpArray + tmpArrayLen;

						memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
						newReg->Sequence = regStart;
						newReg->Offset = (newReg->Sequence - RefSeq);
						newReg->Length = (regEnd - regStart);
						newReg->Type = artValid;
						++tmpArrayLen;
						if (old != NULL)
							utils_free(old);
					}
				}

				if (numberOfNBases == 0)
					regStart = regEnd;

				++regEnd;
				++numberOfNBases;
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

			ret = utils_calloc(tmpArrayLen + 1, sizeof(ACTIVE_REGION), &tmpArray);
			if (ret == ERR_SUCCESS) {
				PACTIVE_REGION newReg = tmpArray + tmpArrayLen;
				
				memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
				newReg->Sequence = regStart;
				newReg->Offset = (newReg->Sequence - RefSeq);
				newReg->Length = (regEnd - regStart);
				newReg->Type = (numberOfNBases == 0) ? artValid : artUnknown;
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


ERR_VALUE input_get_region_by_offset(const PACTIVE_REGION Regions, const size_t Count, const uint64_t Offset, size_t *Index, uint64_t *RegionOffset)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PACTIVE_REGION cur = Regions;
	uint64_t o = Offset;

	if (Offset >= cur->Offset) {
		ret = ERR_OFFSET_TOO_HIGH;
		for (size_t i = 0; i < Count; ++i) {
			if (cur->Length > o) {
				*RegionOffset = o;
				*Index = i;
				ret = ERR_SUCCESS;
				break;
			} else o -= cur->Length;

			++cur;
		}
	} else ret = ERR_OFFSET_TOO_LOW;

	return ret;
}


void input_free_regions(PACTIVE_REGION Regions, const size_t Count)
{
	if (Count > 0)
		utils_free(Regions);

	return;
}
