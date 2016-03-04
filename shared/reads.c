
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "kmer.h"
#include "reads.h"



/************************************************************************/
/*                          HELPER FUNCTIONS                            */
/************************************************************************/

static const char *_sam_read_field(const char *Start)
{
	while (*Start != '\0' && *Start != '\r' && *Start != '\n' && *Start != '\t' && *Start != 26)
		++Start;

	return Start;
}


static const char *_sam_read_string_field(const char *Start, char **String, size_t *Length)
{
	char *tmpString = NULL;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		if (utils_malloc((len + 1)*sizeof(char), &tmpString) == ERR_SUCCESS) {
			memcpy(tmpString, Start, len*sizeof(char));
			tmpString[len] = '\0';
			if (String != NULL) {
				*String = tmpString;
				*Length = len;
			} else utils_free(tmpString);
		} else end = NULL;
	} else end = NULL;

	return end;
}


static const char *_sam_read_uint_field(const char *Start, uint32_t *Value)
{
	uint32_t tmpValue = 0;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		char *tmpEnd = NULL;

		tmpValue = (uint32_t)strtoul(Start, &tmpEnd, 0);
		if (end == tmpEnd) {
			if (Value != NULL)
				*Value = tmpValue;
		} else end = NULL;
	} else end = NULL;

	return end;
}


static const char *_sam_read_int_field(const char *Start, int32_t *Value)
{
	int32_t tmpValue = 0;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		char *tmpEnd = NULL;

		tmpValue = (int32_t)strtol(Start, &tmpEnd, 0);
		if (end == tmpEnd) {
			if (Value != NULL)
				*Value = tmpValue;
		} else end = NULL;
	}
	else end = NULL;

	return end;
}


static void _read_destroy_structure(PONE_READ Read)
{
	if (Read->Quality != NULL)
		utils_free(Read->Quality);

	if (Read->ReadSequence != NULL)
		utils_free(Read->ReadSequence);

	if (Read->CIGAR != NULL)
		utils_free(Read->CIGAR);

	if (Read->TemplateName != NULL)
		utils_free(Read->TemplateName);

	return;
}

/************************************************************************/
/*                          PUBLIC FUNCTIONS                            */
/************************************************************************/

ERR_VALUE read_create_from_test_line(const char *Line, const size_t Length, PONE_READ *Read)
{
	PONE_READ tmpRead = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(1, sizeof(ONE_READ), &tmpRead);
	if (ret == ERR_SUCCESS) {
		memset(tmpRead, 0, sizeof(ONE_READ));
		tmpRead->Pos = (uint64_t)-1;
		tmpRead->ReadSequenceLen = Length;
		ret = utils_calloc(Length + 1, sizeof(char), &tmpRead->ReadSequence);
		if (ret == ERR_SUCCESS) {
			memcpy(tmpRead->ReadSequence, Line, Length*sizeof(char));
			tmpRead->ReadSequence[Length] = '\0';
			*Read = tmpRead;
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpRead);
	}

	return ret;
}


ERR_VALUE read_generate_from_sequence(const char *Seq, const size_t SeqLen, const uint32_t ReadLength, PONE_READ *Read)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = read_set_generate_from_sequence(Seq, SeqLen, ReadLength, 1, Read);

	return ret;
}


ERR_VALUE read_create_from_sam_line(const char *Line, PONE_READ *Read)
{
	PONE_READ tmpRead = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(1, sizeof(ONE_READ), &tmpRead);
	if (ret == ERR_SUCCESS) {
		uint32_t tmp32;

		memset(tmpRead, 0, sizeof(ONE_READ));
		Line = _sam_read_string_field(Line, &tmpRead->TemplateName, &tmpRead->TemplateNameLen);
		if (Line != NULL && *Line == '\t')
			++Line;
		else ret = ERR_SAM_INVALID_RNAME;

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_uint_field(Line, &tmp32);
			if (Line != NULL && *Line == '\t') {
				++Line;
				tmpRead->Flags = tmp32;
			} else ret = ERR_SAM_INVALID_FLAG;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_string_field(Line, NULL, NULL);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_RNAME;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_uint_field(Line, &tmp32);
			if (Line != NULL && *Line == '\t') {
				++Line;
				tmpRead->Pos = tmp32;
			} else ret = ERR_SAM_INVALID_POS;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_uint_field(Line, &tmp32);
			if (Line != NULL && *Line == '\t') {
				++Line;
				tmpRead->PosQuality = (uint8_t)tmp32;
			} else ret = ERR_SAM_INVALID_MAPQ;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_string_field(Line, &tmpRead->CIGAR, &tmpRead->CIGARLen);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_CIGAR;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_string_field(Line, NULL, NULL);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_RNEXT;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_uint_field(Line, NULL);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_PNEXT;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_int_field(Line, NULL);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_TLEN;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_string_field(Line, &tmpRead->ReadSequence, &tmpRead->ReadSequenceLen);
			if (Line != NULL && *Line == '\t')
				++Line;
			else ret = ERR_SAM_INVALID_SEQ;
		}

		if (ret == ERR_SUCCESS) {
			Line = _sam_read_string_field(Line, (char **)&tmpRead->Quality, &tmpRead->QualityLen);
			if (Line == NULL)
				ret = ERR_SAM_INVALID_QUAL;
		}

		if (ret == ERR_SUCCESS) {
			if (tmpRead->ReadSequenceLen == tmpRead->QualityLen) {
				tmpRead->Pos--;
				for (size_t i = 0; i < tmpRead->QualityLen; ++i)
					tmpRead->Quality[i] -= 33;

				*Read = tmpRead;
			} else ret = ERR_SAM_SEQ_QUAL_LEN_MISMATCH;
		}

		if (ret != ERR_SUCCESS) {
			if (tmpRead->Quality != NULL)
				utils_free(tmpRead->Quality);

			if (tmpRead->ReadSequence != NULL)
				utils_free(tmpRead->ReadSequence);

			if (tmpRead->CIGAR != NULL)
				utils_free(tmpRead->CIGAR);
			
			if (tmpRead->TemplateName != NULL)
				utils_free(tmpRead->TemplateName);

			utils_free(tmpRead);
		}
	}

	return ret;
}

ERR_VALUE read_create_from_fasta_seq(const char *Seq, const size_t SeqLen, const char *SeqName, const size_t SeqNameLen, PONE_READ *Read)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	return ret;
}


void read_destroy(PONE_READ Read)
{
	_read_destroy_structure(Read);
	utils_free(Read);

	return;
}


const char *read_get_kmer_pos(const ONE_READ *Read, const KMER *KMer)
{
	const char *ret = NULL;
	PKMER readKMer = NULL;
	KMER_STACK_ALLOC(readKMer, kmer_get_size(KMer), Read->ReadSequence);
	size_t remainingLength = Read->ReadSequenceLen - kmer_get_size(KMer) + 1;
	size_t pos = 0;

	while (remainingLength > 0) {
		if (kmer_equal(readKMer, KMer)) {
			ret = Read->ReadSequence + pos;
			break;
		}

		++pos;
		--remainingLength;
		kmer_advance(readKMer, pos + kmer_get_size(KMer));
	}

	return ret;
}


ERR_VALUE read_set_generate_from_sequence(const char *Seq, const size_t SeqLen, const uint32_t ReadLength, const size_t ReadCount, PONE_READ *ReadSet)
{
	PONE_READ r = NULL;
	PONE_READ tmpReadSet = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(ReadCount, sizeof(ONE_READ), &tmpReadSet);
	if (ret == ERR_SUCCESS) {
		r = tmpReadSet;
		for (size_t i = 0; i < ReadCount; ++i) {
			memset(r, 0, sizeof(ONE_READ));
			r->Pos = utils_ranged_rand(0, SeqLen - ReadLength + 1);
			r->PosQuality = 254;
			r->ReadSequenceLen = ReadLength;
			ret = utils_calloc(r->ReadSequenceLen + 1, sizeof(char), &r->ReadSequence);
			if (ret == ERR_SUCCESS) {
				memcpy(r->ReadSequence, Seq + r->Pos, r->ReadSequenceLen*sizeof(char));
				r->ReadSequence[r->ReadSequenceLen] = '\0';
				r->QualityLen = r->ReadSequenceLen;
				ret = utils_calloc(r->QualityLen, sizeof(uint8_t), &r->Quality);
				if (ret == ERR_SUCCESS)
					memset(r->Quality, 254, r->QualityLen);

				if (ret != ERR_SUCCESS)
					utils_free(r->ReadSequence);
			}

			if (ret != ERR_SUCCESS) {
				--r;
				for (size_t j = 0; j < i; ++j) {
					_read_destroy_structure(r);
					--r;
				}

				break;
			}

			++r;
		}

		if (ret == ERR_SUCCESS)
			*ReadSet = tmpReadSet;

		if (ret != ERR_SUCCESS)
			utils_free(tmpReadSet);
	}

	return ret;
}


void read_set_destroy(PONE_READ ReadSet, const size_t Count)
{
	PONE_READ tmp = ReadSet;

	for (size_t i = 0; i < Count; ++i) {
		_read_destroy_structure(tmp);
		++tmp;
	}

	return;
}


ERR_VALUE read_set_merge(PONE_READ *Target, const size_t TargetCount, struct _ONE_READ *Source, const size_t SourceCount)
{
	PONE_READ tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(TargetCount + SourceCount, sizeof(ONE_READ), &tmp);
	if (ret == ERR_SUCCESS) {
		memcpy(tmp, *Target, TargetCount*sizeof(ONE_READ));
		memcpy(tmp + TargetCount, Source, SourceCount*sizeof(ONE_READ));
		utils_free(Source);
		utils_free(*Target);
		*Target = tmp;
	}

	return ret;
}
