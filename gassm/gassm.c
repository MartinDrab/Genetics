
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
#include "input-file.h"
#include "assembly.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 3);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_FILE, "refseq.txt");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_TYPE, "fasta");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_REFSEQ_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_FILE, "reads.txt");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_TYPE, "none");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_READS_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_OUTPUT_FILE, "-");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_ACTION, "assembly");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_VERBOSE, FALSE);

	if (ret == ERR_SUCCESS) {
		ret = option_set_description_const(GASSM_OPTION_KMER_SIZE, GASSM_OPTION_KMER_SIZE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_REFSEQ_INPUT_FILE, GASSM_OPTION_REFSEQ_INPUT_FIL_DESC);
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
		ret = option_set_description_const(GASSM_OPTION_OUTPUT_FILE, GASSM_OPTION_OUTPUT_FILE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_ACTION, GASSM_OPTION_ACTION_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_VERBOSE, GASSM_OPTION_VERBOSE_DESC);
		assert(ret == ERR_SUCCESS);

		option_set_shortcut(GASSM_OPTION_KMER_SIZE, 'k');
		option_set_shortcut(GASSM_OPTION_REFSEQ_INPUT_FILE, 'i');
		option_set_shortcut(GASSM_OPTION_READS_INPUT_FILE, 'I');
		option_set_shortcut(GASSM_OPTION_REFSEQ_INPUT_TYPE, 't');
		option_set_shortcut(GASSM_OPTION_READS_INPUT_TYPE, 'T');
		option_set_shortcut(GASSM_OPTION_REFSEQ_SKIP_VERT, 's');
		option_set_shortcut(GASSM_OPTION_READS_SKIP_VERT, 'S');
		option_set_shortcut(GASSM_OPTION_OUTPUT_FILE, 'o');
		option_set_shortcut(GASSM_OPTION_ACTION, 'a');
		option_set_shortcut(GASSM_OPTION_VERBOSE, 'v');
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

					ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
					if (ret == ERR_SUCCESS) {
						char *refSeq = NULL;
						size_t refSeqLen = 0;
						char **reads = NULL;
						size_t readCount = 0;

						ret = input_get_refseq(&refSeq, &refSeqLen);
						if (ret == ERR_SUCCESS) {
							ret = input_get_reads(&reads, &readCount);
							if (ret == ERR_SUCCESS) {
								PACTIVE_REGION regions = NULL;
								size_t regionCount = 0;

								ret = input_refseq_to_regions(refSeq, refSeqLen, &regions, &regionCount);
								if (ret == ERR_SUCCESS) {
									char *action = NULL;

									ret = option_get_String(GASSM_OPTION_ACTION, &action);
									if (ret == ERR_SUCCESS) {
										if (strcasecmp(action, "assembly") == 0) {
											PKMER_GRAPH g = NULL;

											ret = kmer_graph_create(kmerSize, &g);
											if (ret == ERR_SUCCESS) {
												ret = kmer_graph_parse_ref_sequence(g, regions[0].Sequence, regions[0].Length);
												if (ret == ERR_SUCCESS)
													ret = kmer_graph_parse_reads(g, reads, readCount);

												if (ret == ERR_SUCCESS) {
//													kmer_graph_delete_1to1_vertices(g);
													kmer_graph_print(stdout, g);
												}

												kmer_graph_destroy(g);
											}
										} else if (strcasecmp(action, "regions") == 0) {
											PACTIVE_REGION ar = regions;

											fprintf(stdout, "%u regions\n", regionCount);
											for (size_t i = 0; i < regionCount; ++i) {
												fprintf(stdout, "  Offset: %llu, length: %llu\n", ar->Offset, ar->Length);
												++ar;
											}
										}
									}

									input_free_regions(regions, regionCount);
								}

								input_free_reads(reads, readCount);
							}

							input_free_refseq(refSeq, refSeqLen);
						}
					}
				}
			}
		}

		options_module_finit();
	}

	return ret;
}
