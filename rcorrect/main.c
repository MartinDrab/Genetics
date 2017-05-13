
#include <stdio.h>
#include <ctype.h>
#include <omp.h>
#include "err.h"
#include "utils.h"
#include "reads.h"
#include "input-file.h"
#include "fml.h"


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



ERR_VALUE convert_to_fermilite(PONE_READ Reads, size_t Count, bseq1_t **Result)
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


ERR_VALUE convert_to_gassm2(bseq1_t *Seqs, size_t Count, PONE_READ Reads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	int i = 0;

#pragma omp parallel for shared(Reads, Seqs)
	for (i = 0; i < (int)Count; ++i) {
		char cigar[20];

		memset(cigar, 0, sizeof(cigar));
		if (Reads[i].Pos > 0)
			snprintf(cigar, sizeof(cigar), "%iM", Seqs[i].l_seq);
		else cigar[0] = "*";

		utils_copy_string(cigar, &Reads[i].Extension->CIGAR);
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


int main(int argc, char **argv)
{
	fml_opt_t options;
	bseq1_t *seqs = NULL;
	PONE_READ reads = NULL;
	size_t readCount = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	fml_opt_init(&options);
	options.n_threads = omp_get_num_procs();
	options.ec_k = 31;
	utils_allocator_init(options.n_threads);
	fprintf(stderr, "Loading reads from %s...\n", argv[1]);
	ret = input_get_reads(argv[1], "sam", &reads, &readCount);
	if (ret == ERR_SUCCESS) {
		fprintf(stderr, "Filtering and adjusting reads...\n");
		input_filter_bad_reads(reads, &readCount, 0, 0);
		fprintf(stderr, "Converting to fermi-lite format...\n");
		ret = convert_to_fermilite(reads, readCount, &seqs);
		if (ret == ERR_SUCCESS) {
			fml_opt_adjust(&options, readCount, seqs);
			fprintf(stderr, "Correcting...\n");
			fml_correct(&options, readCount, seqs);
			fprintf(stderr, "Fitting unique k-mers...\n");
			fml_fltuniq(&options, readCount, seqs);
			fprintf(stderr, "Converting back to our format...\n");
			ret = convert_to_gassm2(seqs, readCount, reads);
			fprintf(stderr, "Saving corrected reads...\n");
			for (size_t i = 0; i < readCount; ++i) {
				if (reads[i].ReadSequenceLen > 0 && reads[i].Quality != NULL)
					read_write_sam(stdout, reads + i);

				read_quality_decode(reads + i);
			}

			fprintf(stderr, "Freeing fermi-lite resources...\n");
			utils_free(seqs);
		}
		
		fprintf(stderr, "Freeing our reads...\n");
//		read_set_destroy(reads, readCount);
	}

	return 0;
}
