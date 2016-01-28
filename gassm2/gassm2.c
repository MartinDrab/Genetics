
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "libkmer.h"
#include "gassm2.h"



static ERR_VALUE _init_default_values()
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(PROGRAM_OPTION_KMERSIZE, 5);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_SEQUENCE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_RANDSEQLEN, 100);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_TEST, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_TEST_COUNT, 256);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_HELP, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_PRINT_RESULTS, FALSE);

	option_set_description_const(PROGRAM_OPTION_KMERSIZE, PROGRAM_OPTION_KMERSIZE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQUENCE, PROGRAM_OPTION_SEQUENCE_DESC);
	option_set_description_const(PROGRAM_OPTION_RANDSEQLEN, PROGRAM_OPTION_RANDSEQLEN_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST, PROGRAM_OPTION_TEST_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_COUNT, PROGRAM_OPTION_TEST_COUNT_DESC);
	option_set_description_const(PROGRAM_OPTION_HELP, PROGRAM_OPTION_HELP_DESC);
	option_set_description_const(PROGRAM_OPTION_PRINT_RESULTS, PROGRAM_OPTION_PRINT_RESULTS_DESC);

	option_set_shortcut(PROGRAM_OPTION_KMERSIZE, 'k');
	option_set_shortcut(PROGRAM_OPTION_SEQUENCE, 's');
	option_set_shortcut(PROGRAM_OPTION_RANDSEQLEN, 'r');
	option_set_shortcut(PROGRAM_OPTION_TEST, 't');
	option_set_shortcut(PROGRAM_OPTION_TEST_COUNT, 'c');
	option_set_shortcut(PROGRAM_OPTION_HELP, 'h');
	option_set_shortcut(PROGRAM_OPTION_PRINT_RESULTS, 'p');

	return ret;
}


static ERR_VALUE _capture_program_options(PPROGRAM_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_get_UInt32(PROGRAM_OPTION_KMERSIZE, &Options->KMerSize);
	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_SEQUENCE, &Options->ReferenceSequence);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_RANDSEQLEN, &Options->TestSeqLen);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_HELP, &Options->Help);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_TEST, &Options->Test);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_PRINT_RESULTS, &Options->PrintResults);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_TEST_COUNT, &Options->TestCount);

	return ret;
}


static void _compute_graph(const struct _PROGRAM_OPTIONS *Options)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	if (*Options->ReferenceSequence != '\0') {
		ret = kmer_graph_create(Options->KMerSize, &g);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_parse_ref_sequence(g, Options->ReferenceSequence, strlen(Options->ReferenceSequence), FALSE);
			if (ret == ERR_SUCCESS) {
				size_t l = 0;
				char *s = NULL;

				ret = kmer_graph_get_seq(g, &s, &l);
				if (ret == ERR_SUCCESS) {
					if (strcmp(Options->ReferenceSequence, s)) {
						fprintf(stderr, "REFSEQ   = %s\n", Options->ReferenceSequence);
						fprintf(stderr, "GRAPHSEQ = %s\n", s);
						kmer_graph_print(stderr, g);
					} else {
						if (Options->PrintResults) {
							printf("REFSEQ   = %s\n", Options->ReferenceSequence);
							printf("GRAPHSEQ = %s\n", s);
							kmer_graph_print(stdout, g);
						}
					}

					utils_free(s);
				} else printf("kmer_graph_get_seq(): %u\n", ret);
			} else printf("kmer_graph_parse_ref_sequence(): %u\n", ret);

			kmer_graph_destroy(g);
		} else printf("kmer_graph_create(): %u\n", ret);
	} else printf("Reference sequence not set\n");

	return;
}


static char _rand_nucleotide(void)
{
	char ret;

	switch (rand() % 4) {
		case 0: ret = 'A'; break;
		case 1: ret = 'C'; break;
		case 2: ret = 'G'; break;
		case 3: ret = 'T'; break;
		default: assert(0); break;
	}

	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _init_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				PROGRAM_OPTIONS po;

				ret = _capture_program_options(&po);
				if (ret == ERR_SUCCESS) {
					if (po.Help) {
						options_print_help();
					} else if (po.Test) {
						char *rs = NULL;

						ret = utils_calloc(po.TestSeqLen + 1, sizeof(char), (char **)&rs);
						if (ret == ERR_SUCCESS) {
							for (uint32_t i = 0; i < po.TestCount; ++i) {
								memset(rs, 0, po.TestSeqLen*sizeof(char));
								for (uint32_t j = 0; j < po.TestSeqLen; ++j)
									rs[j] = _rand_nucleotide();

								po.ReferenceSequence = rs;
								_compute_graph(&po);
							}

							utils_free(rs);
						}
					} else _compute_graph(&po);
				}
			}
		}
	
		options_module_finit();
	}

	return ret;
}