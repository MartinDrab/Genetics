
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
		ret = utils_calloc(readCount, sizeof(bseq1_t), &seqs);
		if (ret == ERR_SUCCESS) {
			int i = 0;

#pragma omp parallel for shared(reads, seqs)
			for (i = 0; i < (int)readCount; ++i) {
				memset(seqs + i, 0, sizeof(seqs[i]));
				seqs[i].l_seq = reads[i].ReadSequenceLen;
				read_quality_encode(reads + i);
				seqs[i].seq = _copy_string(reads[i].ReadSequence, reads[i].ReadSequenceLen);
				if (reads[i].Quality != NULL)
					seqs[i].qual = _copy_string(reads[i].Quality, reads[i].ReadSequenceLen);

				read_quality_decode(reads + i);
				if (reads[i].ReadSequenceLen > 0) {
					utils_free(reads[i].Quality);
					utils_free(reads[i].ReadSequence);
				}
			}

			fml_opt_adjust(&options, readCount, seqs);
			fprintf(stderr, "Correcting...\n");
			fml_correct(&options, readCount, seqs);
			fprintf(stderr, "Fitting unique k-mers...\n");
			fml_fltuniq(&options, readCount, seqs);
			fprintf(stderr, "Converting back to our format...\n");

#pragma omp parallel for shared(reads, seqs)
			for (i = 0; i < (int)readCount; ++i) {
				char cigar[20];
				
				memset(cigar, 0, sizeof(cigar));
				if (reads[i].Pos > 0)
					snprintf(cigar, sizeof(cigar), "%iM", seqs[i].l_seq);
				else cigar[0] = "*";
				
				utils_copy_string(cigar, &reads[i].Extension->CIGAR);
				reads[i].ReadSequenceLen = seqs[i].l_seq;
				ret = utils_copy_string(seqs[i].seq, &reads[i].ReadSequence);
				if (ret == ERR_SUCCESS)
					ret = utils_copy_string(seqs[i].qual, &reads[i].Quality);

				if (seqs[i].l_seq > 0) {
					free(seqs[i].seq);
					free(seqs[i].qual);
				}

				for (size_t j = 0; j < reads[i].ReadSequenceLen; ++j)
					reads[i].ReadSequence[j] = toupper(reads[i].ReadSequence[j]);
			}
				
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
