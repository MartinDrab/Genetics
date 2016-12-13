
#include <stdio.h>
#include <omp.h>
#include "err.h"
#include "utils.h"
#include "reads.h"
#include "input-file.h"
#include "fml.h"






int main(int argc, char **argv)
{
	float coverage = 0.0;
	fml_opt_t options;
	bseq1_t *seqs = NULL;
	PONE_READ reads = NULL;
	size_t readCount = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	fml_opt_init(&options);
	options.n_threads = omp_get_num_threads();
	printf("Loading reads from %s...\n", argv[1]);
	ret = input_get_reads(argv[1], "sam", &reads, &readCount);
	if (ret == ERR_SUCCESS) {
		printf("Converting to fermi-lite format...\n");
		ret = utils_calloc(readCount, sizeof(bseq1_t), &seqs);
		if (ret == ERR_SUCCESS) {
			for (size_t i = 0; i < readCount; ++i) {
				seqs[i].l_seq = reads[i].ReadSequenceLen;
				read_quality_encode(reads + i);
				ret = utils_copy_string(reads[i].ReadSequence, &seqs[i].seq);
				if (ret == ERR_SUCCESS && reads[i].Quality != NULL)
					ret = utils_copy_string(reads[i].Quality, &seqs[i].qual);

				read_quality_decode(reads + i);
			}

			fml_opt_adjust(&options, readCount, seqs);
			printf("Correcting...\n");
			fml_correct(&options, readCount, seqs);
			printf("Fitting unique k-mers...\n");
			coverage = fml_fltuniq(&options, readCount, seqs);
			printf("Converting back to our format...\n");
			for (size_t i = 0; i < readCount; ++i) {
				reads[i].ReadSequenceLen = seqs[i].l_seq;
				reads[i].QualityLen = seqs[i].l_seq;
				ret = utils_copy_string(seqs[i].seq, &reads[i].ReadSequence);
				if (ret == ERR_SUCCESS)
					ret = utils_copy_string(seqs[i].qual, &reads[i].Quality);

				for (size_t j = 0; j < reads[i].ReadSequenceLen; ++j)
					reads[i].ReadSequence[j] = toupper(reads[i].ReadSequence[j]);

				if (reads[i].ReadSequenceLen > 0 && reads[i].QualityLen > 0)
					read_write_sam(stdout, reads + i);

				read_quality_decode(reads + i);
			}
				
			printf("Freeing fermi-lite resources...\n");
			utils_free(seqs);
		}
		
		printf("Freeing our reads...\n");
		read_set_destroy(reads, readCount);
	}

	return 0;
}