
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <stdlib.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "gassm2.h"



static ERR_VALUE _init_default_values()
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(PROGRAM_OPTION_KMERSIZE, 5);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_SEQUENCE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_SEQFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(PROGRAM_OPTION_SEQSTART, (uint64_t)-1);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_SEQLEN, 100);

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
	option_set_description_const(PROGRAM_OPTION_SEQFILE, PROGRAM_OPTION_SEQFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQSTART, PROGRAM_OPTION_SEQSTART_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQLEN, PROGRAM_OPTION_SEQLEN_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST, PROGRAM_OPTION_TEST_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_COUNT, PROGRAM_OPTION_TEST_COUNT_DESC);
	option_set_description_const(PROGRAM_OPTION_HELP, PROGRAM_OPTION_HELP_DESC);
	option_set_description_const(PROGRAM_OPTION_PRINT_RESULTS, PROGRAM_OPTION_PRINT_RESULTS_DESC);

	option_set_shortcut(PROGRAM_OPTION_KMERSIZE, 'k');
	option_set_shortcut(PROGRAM_OPTION_SEQUENCE, 's');
	option_set_shortcut(PROGRAM_OPTION_SEQFILE, 'f');
	option_set_shortcut(PROGRAM_OPTION_SEQSTART, 'S');
	option_set_shortcut(PROGRAM_OPTION_SEQLEN, 'l');
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
	if (ret == ERR_SUCCESS) {
		ret = option_get_String(PROGRAM_OPTION_SEQUENCE, &Options->ReferenceSequence);
		if (ret == ERR_SUCCESS && *Options->ReferenceSequence == '\0') {
			char *seqFile = NULL;
			size_t refSeqLen = 0;

			ret = option_get_String(PROGRAM_OPTION_SEQFILE, &seqFile);
			if (ret == ERR_SUCCESS && *seqFile != '\0')
				ret = input_get_refseq(seqFile, "fasta", &Options->ReferenceSequence, &refSeqLen);
		}
	}

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt64(PROGRAM_OPTION_SEQSTART, &Options->RegionStart);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_SEQLEN, &Options->RegionLength);

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
			ret = kmer_graph_parse_ref_sequence(g, Options->ReferenceSequence, Options->RegionLength, FALSE);
			if (ret == ERR_SUCCESS) {
				size_t l = 0;
				char *s = NULL;

				ret = kmer_graph_get_seq(g, &s, &l);
				if (ret == ERR_SUCCESS) {
					char tmp = Options->ReferenceSequence[Options->RegionLength];
					Options->ReferenceSequence[Options->RegionLength] = '\0';
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

					Options->ReferenceSequence[Options->RegionLength] = tmp;
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


static ERR_VALUE _create_alternatce_sequence(const char *RefSeq, char **Alternate)
{
	return utils_copy_string(RefSeq, Alternate);
}

static ERR_VALUE _generate_reads(const struct _PROGRAM_OPTIONS *Options, const char *Sequence, char ***Reads, size_t *ReadCount)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;



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
						if (*po.ReferenceSequence == '\0') {
							char *rs = NULL;

							ret = utils_calloc(po.RegionLength + 1, sizeof(char), (char **)&rs);
							if (ret == ERR_SUCCESS) {
								printf("Testing with %u random sequences of length %u...\n", po.TestCount, po.RegionLength);
								for (uint32_t i = 0; i < po.TestCount; ++i) {
									memset(rs, 0, po.RegionLength*sizeof(char));
									for (uint32_t j = 0; j < po.RegionLength; ++j)
										rs[j] = _rand_nucleotide();

									po.ReferenceSequence = rs;
									_compute_graph(&po);
								}

								utils_free(rs);
							}
						} else {
							size_t regionCount = 0;
							PACTIVE_REGION regions = NULL;
							size_t refSeqLen = strlen(po.ReferenceSequence);

							ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
							if (ret == ERR_SUCCESS) {
								PACTIVE_REGION pa = regions;
								char *origRefSeq = po.ReferenceSequence;

								printf("Going through a reference sequence with %u regions...\n", regionCount);
								for (size_t i = 0; i < regionCount; ++i) {
									uint64_t remainingLength = 0;

									if (pa->Type == artValid) {
										printf("Region #%u: Offset: %" PRIu64 ", Length %" PRIu64 "\n", i, pa->Offset, pa->Length);
										po.ReferenceSequence = pa->Sequence;
										for (uint64_t j = 0; j < pa->Length; j += po.RegionLength) {
											_compute_graph(&po);
											if (j + po.RegionLength < pa->Length)
												po.ReferenceSequence += po.RegionLength;
											else remainingLength = pa->Length - j;
										}

										if (remainingLength > po.KMerSize) {
											uint32_t tmp = 0;

											tmp = po.RegionLength;
											po.RegionLength = remainingLength;
											_compute_graph(&po);
											po.RegionLength = tmp;
										}
									}

									++pa;
								}

								input_free_regions(regions, regionCount);
							}
						}
					} else {
						size_t regionCount = 0;
						PACTIVE_REGION regions = NULL;
						size_t refSeqLen = strlen(po.ReferenceSequence);

						ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
						if (ret == ERR_SUCCESS) {
							uint32_t index = 0;
							uint64_t regionOffset = 0;

							ret = input_get_region_by_offset(regions, regionCount, po.RegionStart, &index, &regionOffset);
							if (ret == ERR_SUCCESS) {
								PACTIVE_REGION r = regions + index;

								if (r->Type == artValid) {
									po.ReferenceSequence = r->Sequence + regionOffset;
									if (r->Length - regionOffset < po.RegionLength)
										po.RegionLength = r->Length - regionOffset;

									printf("Active region (%" PRIu64 "; %u; %u)...\n", po.RegionStart, po.RegionLength, index);
									_compute_graph(&po);
								} else printf("ERROR: The active region (%" PRIu64 "; %u; %u) does not specify a readable part of the reference sequence\n", po.RegionStart, po.RegionLength, index);
							}

							input_free_regions(regions, regionCount);
						}
					}
				}
			}
		}
	
		options_module_finit();
	}

	return ret;
}
