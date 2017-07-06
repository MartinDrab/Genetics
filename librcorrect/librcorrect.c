
#include <stdio.h>
#include <ctype.h>
#include <omp.h>
#include "err.h"
#include "utils.h"
#include "reads.h"
#include "fml.h"
#include "librcorrect.h"


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

	ret = utils_calloc(Count, sizeof(bseq1_t), &tmpResult);
	if (ret == ERR_SUCCESS) {
		int i = 0;

#pragma omp parallel for shared(Reads, tmpResult)
		for (i = 0; i < (int)Count; ++i) {
			memset(tmpResult + i, 0, sizeof(tmpResult[i]));
			tmpResult[i].l_seq = Reads[i].ReadSequenceLen;
			read_quality_encode(Reads + i);
			tmpResult[i].seq = _copy_string(Reads[i].ReadSequence, Reads[i].ReadSequenceLen);
			if (Reads[i].Quality != NULL)
				tmpResult[i].qual = _copy_string(Reads[i].Quality, Reads[i].ReadSequenceLen);

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
			ret = utils_copy_string(Seqs[i].qual, &Reads[i].Quality);

		if (Seqs[i].l_seq > 0) {
			free(Seqs[i].seq);
			free(Seqs[i].qual);
		}

		for (size_t j = 0; j < Reads[i].ReadSequenceLen; ++j)
			Reads[i].ReadSequence[j] = toupper(Reads[i].ReadSequence[j]);
	}

	return ret;
}


ERR_VALUE libcorrect_correct(PONE_READ Reads, size_t Count, PLIBRCORRECT_STATISTICS Stats)
{
	fml_opt_t options;
	bseq1_t *seqs = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Stats, 0, sizeof(LIBRCORRECT_STATISTICS));
	fml_opt_init(&options);
	Stats->TotalReads = Count;
	for (size_t i = 0; i < Count; ++i) {
		Stats->TotalBases += Reads[i].ReadSequenceLen;
		if (Stats->RepairCountDistributionCount < Reads[i].ReadSequenceLen)
			Stats->RepairCountDistributionCount = Reads[i].ReadSequenceLen;
	}

	ret = utils_calloc(2*Stats->RepairCountDistributionCount, sizeof(uint64_t), &Stats->RepairCountDistribution);
	if (ret == ERR_SUCCESS) {
		Stats->RepairBasePositionDistribution = Stats->RepairCountDistribution + Stats->RepairCountDistributionCount;
		options.n_threads = omp_get_num_procs();
		ret = convert_to_fermilite(Reads, Count, &seqs);
		if (ret == ERR_SUCCESS) {
			fml_opt_adjust(&options, Count, seqs);
			Stats->K = options.ec_k;
			fml_correct(&options, Count, seqs);
			fml_fltuniq(&options, Count, seqs);
			for (size_t i = 0; i < Count; ++i) {
				const bseq1_t *s = seqs + i;
			
				if (s->l_seq > 0) {
					uint32_t repairCount = 0;

					if (s->l_seq != Reads[i].ReadSequenceLen)
						++Stats->ReadsShortened;

					for (size_t j = 0; j < s->l_seq; ++j) {
						if (s->seq[j] == 'a' || s->seq[j] == 'c' ||
							s->seq[j] == 'g' || s->seq[j] == 't') {
							++repairCount;
							++Stats->RepairBasePositionDistribution[s->l_seq - j];
						}
					}

					Stats->TotalRepairs += repairCount;
					++Stats->RepairCountDistribution[repairCount];
				} else ++Stats->ReadsRemoved;
			}
			
			ret = convert_to_gassm2(seqs, Count, Reads);
			if (ret == ERR_SUCCESS) {
				for (size_t i = 0; i < Count; ++i)
					read_quality_decode(Reads + i);
			}

			utils_free(seqs);
		}

		if (ret != ERR_SUCCESS)
			utils_free(Stats->RepairCountDistribution);
	}

	return ret;
}
