
#include <inttypes.h>
#include <stdio.h>
#include <ctype.h>
#include <omp.h>
#include "err.h"
#include "utils.h"
#include "reads.h"
#include "fml.h"
#include "librcorrect.h"


UTILS_TYPED_CALLOC_FUNCTION(bseq1_t)



static char *_copy_string(const char *str, const size_t len)
{
	char *ret = NULL;

	ret = malloc((len + 1)*sizeof(char));
	if (ret != NULL) {
		memcpy(ret, str, len*sizeof(char));
		ret[len] = '\0';
	}

	return ret;
}


static ERR_VALUE convert_to_fermilite(PONE_READ Reads, size_t Count, bseq1_t **Result)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	bseq1_t *tmpResult = NULL;

	ret = utils_calloc_bseq1_t(Count, &tmpResult);
	if (ret == ERR_SUCCESS) {
		int i = 0;

#pragma omp parallel for shared(Reads, tmpResult)
		for (i = 0; i < (int)Count; ++i) {
			memset(tmpResult + i, 0, sizeof(tmpResult[i]));
			tmpResult[i].l_seq = Reads[i].ReadSequenceLen;
			read_quality_encode(Reads + i);
			tmpResult[i].seq = _copy_string(Reads[i].ReadSequence, Reads[i].ReadSequenceLen);
			if (Reads[i].Quality != NULL)
				tmpResult[i].qual = _copy_string((char *)Reads[i].Quality, Reads[i].ReadSequenceLen);

			read_quality_decode(Reads + i);
			if (Reads[i].ReadSequenceLen > 0) {
				utils_free(Reads[i].Quality);
				utils_free(Reads[i].ReadSequence);
			}
		}

		if (ret == ERR_SUCCESS)
			*Result = tmpResult;
		
		if (ret != ERR_SUCCESS)
			utils_free(tmpResult);
	}

	return ret;
}


static ERR_VALUE convert_to_gassm2(const bseq1_t *Seqs, size_t Count, PONE_READ Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int i = 0;

#pragma omp parallel for shared(Reads, Seqs)
	for (i = 0; i < (int)Count; ++i) {
		if (Seqs[i].l_seq != Reads[i].ReadSequenceLen) {
			char cigar[20];

			memset(cigar, 0, sizeof(cigar));
			if (Reads[i].Pos > 0)
				snprintf(cigar, sizeof(cigar), "%iM", Seqs[i].l_seq);
			else cigar[0] = '*';

			utils_copy_string(cigar, &Reads[i].Extension->CIGAR);
		}

		Reads[i].ReadSequenceLen = Seqs[i].l_seq;
		ret = utils_copy_string(Seqs[i].seq, &Reads[i].ReadSequence);
		if (ret == ERR_SUCCESS)
			ret = utils_copy_string(Seqs[i].qual, (uint8_t **)&Reads[i].Quality);

		if (Seqs[i].l_seq > 0) {
			free(Seqs[i].seq);
			free(Seqs[i].qual);
		}

		for (uint32_t j = 0; j < Reads[i].ReadSequenceLen; ++j)
			Reads[i].ReadSequence[j] = toupper(Reads[i].ReadSequence[j]);
	}

	return ret;
}


ERR_VALUE libcorrect_state_init(PLIBCORRECT_STATE State, PONE_READ Reads, const size_t Count)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_malloc(sizeof(fml_opt_t), &State->Options);
	if (ret == ERR_SUCCESS) {
		bseq1_t *seqs = NULL;
		fml_opt_t *options = (fml_opt_t *)State->Options;

		State->Count = Count;
		State->Reads = Reads;
		fml_opt_init(options);
		options->n_threads = omp_get_num_procs();
		ret = convert_to_fermilite(Reads, Count, &seqs);
		if (ret == ERR_SUCCESS) {
			State->ConvertedReads = seqs;
			fml_opt_adjust(options, Count, seqs);
		}
			
		if (ret != ERR_SUCCESS)
			utils_free(State->Options);
	}

	return ret;
}


void libcorrect_correct(PLIBCORRECT_STATE State)
{
	bseq1_t *seqs = (bseq1_t *)State->ConvertedReads;
	const fml_opt_t *options = (fml_opt_t *)State->Options;

	fml_correct(options, State->Count, seqs);
	fml_fltuniq(options, State->Count, seqs);

	return;
}


ERR_VALUE libcorrect_correct_stats(const LIBCORRECT_STATE *State, PLIBRCORRECT_STATISTICS Stats)
{
	bseq1_t *seqs = (bseq1_t *)State->ConvertedReads;
	const fml_opt_t *options = (fml_opt_t *)State->Options;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Stats, 0, sizeof(LIBRCORRECT_STATISTICS));
	Stats->K = options->ec_k;
	Stats->TotalReads = State->Count;
	for (size_t i = 0; i < State->Count; ++i) {
		Stats->TotalBases += State->Reads[i].ReadSequenceLen;
		if (Stats->RepairCountDistributionCount < State->Reads[i].ReadSequenceLen)
			Stats->RepairCountDistributionCount = State->Reads[i].ReadSequenceLen;
	}

	ret = utils_calloc_uint64_t(2 * Stats->RepairCountDistributionCount, &Stats->RepairCountDistribution);
	if (ret == ERR_SUCCESS) {
		Stats->RepairBasePositionDistribution = Stats->RepairCountDistribution + Stats->RepairCountDistributionCount;
		for (size_t i = 0; i < State->Count; ++i) {
			const bseq1_t *s = seqs + i;

			if (s->l_seq > 0) {
				uint32_t repairCount = 0;

				if (s->l_seq != State->Reads[i].ReadSequenceLen)
					++Stats->ReadsShortened;

				for (int32_t j = 0; j < s->l_seq; ++j) {
					if (s->seq[j] == 'a' || s->seq[j] == 'c' ||
						s->seq[j] == 'g' || s->seq[j] == 't') {
						++repairCount;
						++Stats->RepairBasePositionDistribution[s->l_seq - j - 1];
					}
				}

				Stats->TotalRepairs += repairCount;
				++Stats->RepairCountDistribution[repairCount];
			} else ++Stats->ReadsRemoved;
		}
	}

	return ret;
}


void librcorrect_stats_free(PLIBRCORRECT_STATISTICS Stats)
{
	utils_free(Stats->RepairCountDistribution);
	memset(Stats, 0, sizeof(LIBRCORRECT_STATISTICS));

	return;
}


ERR_VALUE librcorrect_state_finit(PLIBCORRECT_STATE State)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	bseq1_t *seqs = (bseq1_t *)State->ConvertedReads;

	ret = convert_to_gassm2(seqs, State->Count, State->Reads);
	if (ret == ERR_SUCCESS) {
		for (size_t i = 0; i < State->Count; ++i)
			read_quality_decode(State->Reads + i);
	}

	utils_free(seqs);

	return ret;
}


void libcorrect_stats_print(FILE *Stream, const LIBRCORRECT_STATISTICS *Stats)
{
	fprintf(Stream, "K:                  %u\n", Stats->K);
	fprintf(Stream, "Total reads:        %" PRIu64 "\n", Stats->TotalReads);
	fprintf(Stream, "Removed reads:      %" PRIu64 "\n", Stats->ReadsRemoved);
	fprintf(Stream, "Shortened reads:    %" PRIu64 "\n", Stats->ReadsShortened);
	fprintf(Stream, "Total bases:        %" PRIu64 "\n", Stats->TotalBases);
	fprintf(Stream, "Total repairs:      %" PRIu64 "\n", Stats->TotalRepairs);
	fputs("Repair count distribution:\n", stderr);
	for (uint32_t i = 0; i < Stats->RepairCountDistributionCount; ++i) {
		if (Stats->RepairCountDistribution[i] > 0)
			fprintf(Stream, "%u,\t%" PRIu64 "\t%.3lf %%\n", i, Stats->RepairCountDistribution[i], (double)Stats->RepairCountDistribution[i] * 100 / Stats->TotalReads);
	}

	fputs("Repair base position distribution:\n", Stream);
	for (uint32_t i = 0; i < Stats->RepairCountDistributionCount; ++i) {
		if (Stats->RepairBasePositionDistribution[i] > 0)
			fprintf(Stream, "%u,\t%" PRIu64 "\t%.3lf %%\n", i, Stats->RepairBasePositionDistribution[i], (double)Stats->RepairBasePositionDistribution[i] * 100 / Stats->TotalRepairs);
	}

	return;
}
