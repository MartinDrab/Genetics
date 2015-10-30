
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "kmer-graph.h"
#include "test-file.h"
#include "assembly.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 3);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_FIL, "refseq.txt");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_TYPE, "test");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_REFSEQ_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_FILE, "reads.txt");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_TYPE, "test");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_READS_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS) {
		ret = option_set_description_const(GASSM_OPTION_KMER_SIZE, GASSM_OPTION_KMER_SIZE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_REFSEQ_INPUT_FIL, GASSM_OPTION_REFSEQ_INPUT_FIL_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_REFSEQ_INPUT_TYPE, GASSM_OPTION_REFSEQ_INPUT_TYPE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_REFSEQ_SKIP_VERT, GASSM_OPTION_REFSEQ_SKIP_VERT_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_READS_INPUT_FILE, GASSM_OPTION_READS_INPUT_FILE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_READS_INPUT_TYPE, GASSM_OPTION_READS_INPUT_TYPE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_READS_SKIP_VERT, GASSM_OPTION_READS_SKIP_VERT_DESC);
		assert(ret == ERR_SUCCESS);
	}

	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _set_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				if (argc == 1 || (argc == 2 && (strcasecmp(argv[1], "--help") == 0 ||
					strcasecmp(argv[1], "-h") == 0 || strcasecmp(argv[1], "/?") == 0 ||
					strcasecmp(argv[1], "-?") == 0))) {
					printf("Usage: gassm [options]\n");
					options_print_help();
				} else {
					uint32_t kmerSize = 0;
					char *testFile = NULL;

					ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_TEST_INPUT_FILE, &testFile);

					if (ret == ERR_SUCCESS) {
						char *testData = NULL;
						size_t testDataSize = 0;

						ret = utils_file_read(testFile, &testData, &testDataSize);
						if (ret == ERR_SUCCESS) {
							char *refSeq = NULL;
							char **reads = NULL;
							size_t readCount = 0;

							ret = parse_test_data(testData, testDataSize, &refSeq, &reads, &readCount);
							if (ret == ERR_SUCCESS) {
								PKMER_GRAPH g = NULL;

								ret = kmer_graph_create(kmerSize, &g);
								if (ret == ERR_SUCCESS) {
									ret = kmer_graph_parse_ref_sequence(g, refSeq);
									if (ret == ERR_SUCCESS) {
										ret = kmer_graph_parse_reads(g, reads, readCount);
										if (ret == ERR_SUCCESS)
											kmer_graph_print(g);
									}

									kmer_graph_destroy(g);
								}

								free_test_data(refSeq, reads, readCount);
							}

							free(testData);
						}
					}
				}
			}
		}

		options_module_finit();
	}

	return ret;
}
