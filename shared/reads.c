
#include <stdlib.h>
#include "err.h"
#include "utils.h"
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

/************************************************************************/
/*                          PUBLIC FUNCTIONS                            */
/************************************************************************/

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
	if (Read->Quality != NULL)
		utils_free(Read->Quality);

	if (Read->ReadSequence != NULL)
		utils_free(Read->ReadSequence);

	if (Read->CIGAR != NULL)
		utils_free(Read->CIGAR);

	if (Read->TemplateName != NULL)
		utils_free(Read->TemplateName);

	utils_free(Read);

	return;
}
