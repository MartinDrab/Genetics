
#include <omp.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "reads.h"
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

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_READ_COUNT, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_READ_LENGTH, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_TEST_READ_CYCLES, 20);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_TEST_STEP, 1500);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_THRESHOLD, 0);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_READFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_Double(PROGRAM_OPTION_SNP_RATIO, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_Double(PROGRAM_OPTION_INSERT_RATIO, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_Double(PROGRAM_OPTION_DELETE_RATIO, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_ALT1_SEQ, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_ALT2_SEQ, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_CONNECT_READS, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_DISTINCT_PASSES, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_MERGE_UNBRANCHED, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_RESOLVE_BUBBLES, FALSE);


	option_set_description_const(PROGRAM_OPTION_KMERSIZE, PROGRAM_OPTION_KMERSIZE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQUENCE, PROGRAM_OPTION_SEQUENCE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQFILE, PROGRAM_OPTION_SEQFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQSTART, PROGRAM_OPTION_SEQSTART_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQLEN, PROGRAM_OPTION_SEQLEN_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST, PROGRAM_OPTION_TEST_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_COUNT, PROGRAM_OPTION_TEST_COUNT_DESC);
	option_set_description_const(PROGRAM_OPTION_HELP, PROGRAM_OPTION_HELP_DESC);
	option_set_description_const(PROGRAM_OPTION_PRINT_RESULTS, PROGRAM_OPTION_PRINT_RESULTS_DESC);
	option_set_description_const(PROGRAM_OPTION_READ_COUNT, PROGRAM_OPTION_READ_COUNT_DESC);
	option_set_description_const(PROGRAM_OPTION_READ_LENGTH, PROGRAM_OPTION_READ_LENGTH_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_READ_CYCLES, PROGRAM_OPTION_TEST_READ_CYCLES_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_STEP, PROGRAM_OPTION_TEST_STEP_DESC);
	option_set_description_const(PROGRAM_OPTION_THRESHOLD, PROGRAM_OPTION_THRESHOLD_DESC);
	option_set_description_const(PROGRAM_OPTION_READFILE, PROGRAM_OPTION_READFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_SNP_RATIO, PROGRAM_OPTION_SNP_RATIO_DESC);
	option_set_description_const(PROGRAM_OPTION_INSERT_RATIO, PROGRAM_OPTION_INSERT_RATIO_DESC);
	option_set_description_const(PROGRAM_OPTION_DELETE_RATIO, PROGRAM_OPTION_DELETE_RATIO_DESC);
	option_set_description_const(PROGRAM_OPTION_ALT1_SEQ, PROGRAM_OPTION_ALT1_SEQ_DESC);
	option_set_description_const(PROGRAM_OPTION_ALT2_SEQ, PROGRAM_OPTION_ALT2_SEQ_DESC);
	option_set_description_const(PROGRAM_OPTION_DISTINCT_PASSES, PROGRAM_OPTION_DISTINCT_PASSES_DESC);
	option_set_description_const(PROGRAM_OPTION_CONNECT_READS, PROGRAM_OPTION_CONNECT_READS_DESC);
	option_set_description_const(PROGRAM_OPTION_RESOLVE_BUBBLES, PROGRAM_OPTION_RESOLVE_BUBBLES_DESC);
	option_set_description_const(PROGRAM_OPTION_MERGE_UNBRANCHED, PROGRAM_OPTION_MERGE_UNBRANCHED_DESC);

	option_set_shortcut(PROGRAM_OPTION_KMERSIZE, 'k');
	option_set_shortcut(PROGRAM_OPTION_SEQUENCE, 's');
	option_set_shortcut(PROGRAM_OPTION_SEQFILE, 'f');
	option_set_shortcut(PROGRAM_OPTION_SEQSTART, 'S');
	option_set_shortcut(PROGRAM_OPTION_SEQLEN, 'l');
	option_set_shortcut(PROGRAM_OPTION_TEST, 't');
	option_set_shortcut(PROGRAM_OPTION_TEST_COUNT, 'c');
	option_set_shortcut(PROGRAM_OPTION_READ_LENGTH, 'L');
	option_set_shortcut(PROGRAM_OPTION_READ_COUNT, 'C');
	option_set_shortcut(PROGRAM_OPTION_TEST_READ_CYCLES, 'T');
	option_set_shortcut(PROGRAM_OPTION_HELP, 'h');
	option_set_shortcut(PROGRAM_OPTION_PRINT_RESULTS, 'p');
	option_set_shortcut(PROGRAM_OPTION_TEST_STEP, 'e');
	option_set_shortcut(PROGRAM_OPTION_THRESHOLD, 'w');
	option_set_shortcut(PROGRAM_OPTION_READFILE, 'F');
	option_set_shortcut(PROGRAM_OPTION_ALT1_SEQ, 'a');
	option_set_shortcut(PROGRAM_OPTION_ALT2_SEQ, 'A');
	option_set_shortcut(PROGRAM_OPTION_SNP_RATIO, 'r');
	option_set_shortcut(PROGRAM_OPTION_INSERT_RATIO, 'i');
	option_set_shortcut(PROGRAM_OPTION_DELETE_RATIO, 'd');
	option_set_shortcut(PROGRAM_OPTION_DISTINCT_PASSES, '1');
	option_set_shortcut(PROGRAM_OPTION_CONNECT_READS, '2');
	option_set_shortcut(PROGRAM_OPTION_RESOLVE_BUBBLES, '3');
	option_set_shortcut(PROGRAM_OPTION_MERGE_UNBRANCHED, '4');


	return ret;
}


static ERR_VALUE _capture_program_options(PPROGRAM_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Options, 0, sizeof(PROGRAM_OPTIONS));
	ret = option_get_UInt32(PROGRAM_OPTION_KMERSIZE, &Options->KMerSize);
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

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_READ_COUNT, &Options->ReadCount);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_READ_LENGTH, &Options->ReadLength);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_TEST_READ_CYCLES, &Options->TestReadCycles);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_TEST_STEP, &Options->TestStep);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_THRESHOLD, &Options->Threshold);

	if (ret == ERR_SUCCESS)
		ret = option_get_Double(PROGRAM_OPTION_SNP_RATIO, &Options->SNPRatio);

	if (ret == ERR_SUCCESS)
		ret = option_get_Double(PROGRAM_OPTION_INSERT_RATIO, &Options->InsertRatio);

	if (ret == ERR_SUCCESS)
		ret = option_get_Double(PROGRAM_OPTION_DELETE_RATIO, &Options->DeleteRatio);

	if (ret == ERR_SUCCESS) {
		char *rs = NULL;
		ret = option_get_String(PROGRAM_OPTION_SEQUENCE, &rs);
		Options->ReferenceSequence = rs;
		if (ret == ERR_SUCCESS && *Options->ReferenceSequence == '\0')
			ret = option_get_String(PROGRAM_OPTION_SEQFILE, &Options->RefSeqFile);
	}

	if (ret == ERR_SUCCESS) {
		size_t readCount = 0;
		char *readFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_READFILE, &readFile);
		if (ret == ERR_SUCCESS && *readFile != '\0') {
			ret = input_get_reads(readFile, "sam", Options->RegionStart, Options->RegionLength, &Options->Reads, &readCount);
			if (ret == ERR_SUCCESS)
				Options->ReadCount = (uint32_t)readCount;
		}
	}

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_ALT1_SEQ, &Options->AltenrateSequence1);

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_ALT2_SEQ, &Options->AlternateSequence2);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_DISTINCT_PASSES, &Options->MakeDistinctPasses);;

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_CONNECT_READS, &Options->ConnectReads);;

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_RESOLVE_BUBBLES, &Options->ResolveBubbles);;

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(PROGRAM_OPTION_MERGE_UNBRANCHED, &Options->MergeUnbranched);;

	return ret;
}


static void _compare_alternate_sequences(const PROGRAM_OPTIONS *Options, const KMER_GRAPH *Graph, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	boolean notFound = FALSE;
	GEN_ARRAY_PFOUND_SEQUENCE seqArray;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *alternates[2];
	size_t alternateLens[2];

	alternates[0] = Task->Alternate1;
	alternates[1] = Task->Alternate2;
	alternateLens[0] = Task->Alternate1Length;
	alternateLens[1] = Task->Alternate2Length;
	dym_array_init_PFOUND_SEQUENCE(&seqArray, 140);
	ret = kmer_graph_get_seqs(Graph, &seqArray);
	if (ret == ERR_SUCCESS) {
		printf("%u\n", (uint32_t)gen_array_size(&seqArray));
		for (size_t i = 0; i < sizeof(alternateLens) / sizeof(size_t); ++i) {
			boolean found = FALSE;

			for (size_t j = 0; j < gen_array_size(&seqArray); ++j) {
				const FOUND_SEQUENCE *fs = *dym_array_item_PFOUND_SEQUENCE(&seqArray, j);

				found = found_sequence_match(fs, alternates[i], alternateLens[i]);
				if (found)
					break;
			}

			if (!found) {
				notFound = TRUE;
//				fprintf(stderr, "RS:   %s\n", Task->Reference);
//				fprintf(stderr, "ALT1: %s\n", Task->Alternate1);
//				fprintf(stderr, "ALT2: %s\n", Task->Alternate2);
//				kmer_graph_print(stderr, Graph);
//				exit(0);
				break;
			}
		}
	} else printf("ERROR: kmer_graph_get_seqs(): %u\n", ret);
	
	for (size_t j = 0; j < gen_array_size(&seqArray); ++j)
		found_sequence_free(*dym_array_item_PFOUND_SEQUENCE(&seqArray, j));

	dym_array_finit_PFOUND_SEQUENCE(&seqArray);
	if (notFound) {
		++Statistics->FailureCount;
		printf("FAILD\n");
	} else {
		++Statistics->SuccessCount;
		printf("OK\n");
	}

	return;
}


static void _compute_graph(const PROGRAM_OPTIONS *Options, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	Statistics->SuccessCount = 0;
	Statistics->FailureCount = 0;
	Statistics->CannotSucceed = 0;
		ret = kmer_graph_create(Options->KMerSize, &g);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_parse_ref_sequence(g, Task->Reference, Task->ReferenceLength, Options->Threshold);
			if (ret == ERR_SUCCESS) {
				ret = kmer_graph_parse_reads(g, Task->Reads, Task->ReadCount);
				if (ret == ERR_SUCCESS) {
					size_t deletedThings = 0;

					kmer_graph_delete_edges_under_threshold(g, Options->Threshold);
					kmer_graph_delete_trailing_things(g, &deletedThings);
					if (deletedThings == 0) {
						ret = kmer_graph_connect_reads_by_reads(g, Options->Threshold);
						if (ret == ERR_SUCCESS) {
							ret = kmer_graph_connect_reads_by_refseq(g, Options->Threshold);
							if (ret == ERR_SUCCESS) {
								kmer_graph_delete_1to1_vertices(g);
//								ret = kmer_graph_resolve_bubbles(g, Options->Threshold);
								if (ret == ERR_SUCCESS) {
									ret = kmer_graph_detect_uncertainities(g);
									if (ret == ERR_SUCCESS)
										_compare_alternate_sequences(Options, g, Task, Statistics);
									else printf("ERROR: kmer_graph_detect_uncertainities(): %u\n", ret);
								} else printf("ERROR: kmer_graph_resolve_bubbles(): %u\n", ret);
							} else printf("kmer_graph_connect_reads(): %u\n", ret);
						} else printf("kmer_graph_connect_reads_by_reads(): %u\n", ret);
					} else {
						printf("kmer_graph_delete_trailing_things(): deleted something, cannot proceed successfully\n");
						Statistics->CannotSucceed++;
					}
				} else printf("kmer_graph_parse_reads(): %u\n", ret);
			} else printf("kmer_graph_parse_ref_sequence(): %u\n", ret);

			kmer_graph_destroy(g);
		} else printf("kmer_graph_create(): %u\n", ret);

		if (ret != ERR_SUCCESS) {
			++Statistics->FailureCount;
			printf("FAILD\n");
		}

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


static ERR_VALUE _create_alternatce_sequence(const PROGRAM_OPTIONS *Options, const char *RefSeq, const size_t RefSeqLen, char **Alternate, size_t *AlternateLen)
{
	size_t tmpAlternateLen = RefSeqLen*2;
	char *rsCopy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(tmpAlternateLen + 1, sizeof(char), &rsCopy);
	if (ret == ERR_SUCCESS) {
		tmpAlternateLen /= 2;
		memcpy(rsCopy, RefSeq, RefSeqLen*sizeof(char));
		if (Options->SNPRatio > 0 || Options->InsertRatio > 0 || Options->DeleteRatio > 0) {
			size_t spmBorder = ceil(Options->SNPRatio*RefSeqLen);
			size_t insertBorder = ceil(Options->InsertRatio*RefSeqLen);
			size_t deleteBorder = ceil(Options->DeleteRatio*RefSeqLen);
			DYM_ARRAY insertPositions;
			DYM_ARRAY deletePositions;
			
			dym_array_create(&insertPositions, 140);
			dym_array_create(&deletePositions, 140);
			for (size_t i = Options->KMerSize; i < RefSeqLen - Options->KMerSize; ++i) {
				if (utils_ranged_rand(0, RefSeqLen) < spmBorder) {
					char old = rsCopy[i];

					rsCopy[i] = _rand_nucleotide();
					if (rsCopy[i] != old)
						putchar('M');
				}

				if (utils_ranged_rand(0, RefSeqLen) < insertBorder) {
					ret = dym_array_push_back(&insertPositions, (void *)i);
					putchar('I');
				}

				if (ret == ERR_SUCCESS) {
					if (utils_ranged_rand(0, RefSeqLen) < deleteBorder) {
						ret = dym_array_push_back(&deletePositions, (void *)i);
						putchar('D');
					}
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			if (ret == ERR_SUCCESS) {
				tmpAlternateLen = RefSeqLen + dym_array_size(&insertPositions) - dym_array_size(&deletePositions);
				for (size_t i = 0; i < dym_array_size(&insertPositions); ++i) {
					size_t pos = (size_t)dym_array_get(&insertPositions, i) - i;

					memmove(rsCopy + pos + 1, rsCopy + pos, RefSeqLen - pos + 1);
					rsCopy[pos] = _rand_nucleotide();
					for (size_t j = 0; j < dym_array_size(&deletePositions); ++j) {
						size_t *ppos = (size_t *)(dym_array_data(&deletePositions) + j);

						if (*ppos >= pos)
							*ppos++;
					}
				}

				for (size_t i = 0; i < dym_array_size(&deletePositions); ++i) {
					size_t pos = (size_t)dym_array_get(&deletePositions, i) + i;

					memmove(rsCopy + pos, rsCopy + pos + 1, RefSeqLen - pos);
				}
			}

			dym_array_destroy(&deletePositions);
			dym_array_destroy(&insertPositions);
		}
		
		if (Options->SNPRatio > 0)
			printf("\n");

		if (ret == ERR_SUCCESS) {
			rsCopy[tmpAlternateLen] = '\0';
			*Alternate = rsCopy;
			*AlternateLen = tmpAlternateLen;
		}

		if (ret != ERR_SUCCESS)
			utils_free(rsCopy);
	}

	return ret;
}


static ERR_VALUE _test_with_reads(PPROGRAM_OPTIONS Options, const char *RefSeq, PPROGRAM_STATISTICS Statistics)
{
	ASSEMBLY_TASK task;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	memset(Statistics, 0, sizeof(PROGRAM_STATISTICS));
	if (Options->ReadCount > 0 && Options->TestReadCycles > 0) {
		for (size_t j = 0; j < Options->TestReadCycles; ++j) {
			char *alternate = NULL;
			size_t alternateLen = 0;

			ret = ERR_SUCCESS;
			if (*Options->AltenrateSequence1 != '\0') {
				alternate = Options->AltenrateSequence1;
				alternateLen = strlen(Options->AltenrateSequence1);
			} else ret = _create_alternatce_sequence(Options, RefSeq, Options->RegionLength, &alternate, &alternateLen);
			
			if (ret == ERR_SUCCESS) {
				PONE_READ reads1 = NULL;

				ret = read_set_generate_from_sequence(alternate, alternateLen, Options->ReadLength, Options->ReadCount / 2, &reads1);
				if (ret == ERR_SUCCESS) {
					char *alternate2 = NULL;
					size_t alternateLen2 = 0;

					if (*Options->AlternateSequence2 != '\0') {
						alternate2 = Options->AlternateSequence2;
						alternateLen2 = strlen(alternate2);
					} else ret = _create_alternatce_sequence(Options, RefSeq, Options->RegionLength, &alternate2, &alternateLen2);
					
					if (ret == ERR_SUCCESS) {
						PONE_READ reads2 = NULL;

						ret = read_set_generate_from_sequence(alternate2, alternateLen2, Options->ReadLength, Options->ReadCount / 2, &reads2);
						if (ret == ERR_SUCCESS) {
							PONE_READ finalReadSet = reads1;

							ret = read_set_merge(&finalReadSet, Options->ReadCount / 2, reads2, Options->ReadCount / 2);
							if (ret == ERR_SUCCESS) {
								char *rs = NULL;
								PROGRAM_STATISTICS stats;

								reads1 = NULL;
								reads2 = NULL;
								ret = utils_calloc(Options->RegionLength + 1, sizeof(char), &rs);
								if (ret == ERR_SUCCESS) {
									memcpy(rs, RefSeq, Options->RegionLength*sizeof(char));
									rs[Options->RegionLength] = '\0';
									assembly_task_init(&task, rs, Options->RegionLength, alternate, alternateLen, alternate2, alternateLen2, finalReadSet, Options->ReadCount);
									_compute_graph(Options, &task, &stats);
									assembly_task_finit(&task);
									Statistics->FailureCount += stats.FailureCount;
									Statistics->SuccessCount += stats.SuccessCount;
									Statistics->CannotSucceed += stats.CannotSucceed;
									utils_free(rs);
								}

								read_set_destroy(finalReadSet, Options->ReadCount);
							}

							if (reads2 != NULL)
								read_set_destroy(reads2, Options->ReadCount / 2);
						}

						if (*Options->AlternateSequence2 == '\0')
							utils_free(alternate2);
					}

					if (reads1 != NULL)
						read_set_destroy(reads1, Options->ReadCount / 2);
				}

				if (*Options->AltenrateSequence1 == '\0')
					utils_free(alternate);
			}

			if (ret != ERR_SUCCESS)
				break;
		}
	} else {
		assembly_task_init(&task, RefSeq, Options->RegionLength, RefSeq, Options->RegionLength, RefSeq, Options->RegionLength, NULL, 0);
		_compute_graph(Options, &task, Statistics);
		assembly_task_finit(&task);
	}

	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	omp_set_num_threads(1);
	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _init_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				PROGRAM_OPTIONS po;
				PROGRAM_STATISTICS st;

				memset(&st, 0, sizeof(st));
				ret = _capture_program_options(&po);
				if (ret == ERR_SUCCESS) {
					if (po.Help) {
						options_print_help();
					} else if (po.Test) {
						printf("kmer size: %u\n", po.KMerSize);
						if (*po.ReferenceSequence == '\0' && *po.RefSeqFile == '\0') {
							char *rs = NULL;

							ret = utils_calloc(po.RegionLength + 1, sizeof(char), (char **)&rs);
							if (ret == ERR_SUCCESS) {
								printf("kmer size: %u\n", po.KMerSize);
								printf("Testing with %u random sequences of length %u...\n", po.TestCount, po.RegionLength);
								printf("%u test read cycles with %u reads of length %u...\n", po.TestReadCycles, po.ReadCount, po.ReadLength);
								for (uint32_t i = 0; i < po.TestCount; ++i) {									
									PROGRAM_STATISTICS tmpStats;
									
									memset(rs, 0, po.RegionLength*sizeof(char));
									for (uint32_t j = 0; j < po.RegionLength; ++j)
										rs[j] = _rand_nucleotide();

									po.ReferenceSequence = rs;
									ret = _test_with_reads(&po, rs, &tmpStats);
									if (ret == ERR_SUCCESS) {
										st.FailureCount += tmpStats.FailureCount;
										st.SuccessCount += tmpStats.SuccessCount;
										st.CannotSucceed += tmpStats.CannotSucceed;
									} else printf("_test_with_reads(): %u", ret);
								}

								if (ret == ERR_SUCCESS)
									printf("Success: (%" PRIu64 "), Failures: (%" PRIu64 "), Not tried: (%" PRIu64 "), Percentage: (%" PRIu64 ")\n", st.SuccessCount, st.FailureCount, st.CannotSucceed, (uint64_t)(st.SuccessCount * 100 / (st.SuccessCount + st.FailureCount + st.CannotSucceed)));

								utils_free(rs);
							}
						} else {
							size_t refSeqLen = 0;
							boolean explicitSequence = (*po.ReferenceSequence != '\0');
							FASTA_FILE seqFile;

							if (!explicitSequence) {
								ret = fasta_load(po.RefSeqFile, &seqFile);
								if (ret == ERR_SUCCESS) {
									char *rs = NULL;
									ret = fasta_read_seq(&seqFile, &rs, &refSeqLen);
									po.ReferenceSequence = rs;
									if (ret != ERR_SUCCESS)
										fasta_free(&seqFile);
								}
							} else refSeqLen = strlen(po.ReferenceSequence);

							if (ret == ERR_SUCCESS) {
								uint64_t numberOfAttempts = 0;

								do {
									size_t regionCount = 0;
									PACTIVE_REGION regions = NULL;
									
									ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
									if (ret == ERR_SUCCESS) {
										printf("Going through a reference sequence of length %" PRIu64 " MBases with %u regions...\n", (uint64_t)refSeqLen / 1000000, regionCount);
										printf("%u test read cycles with %u reads of length %u...\n", po.TestReadCycles, po.ReadCount, po.ReadLength);
										for (size_t i = 0; i < regionCount; ++i) {
											PACTIVE_REGION pa = regions + i;

											printf("Region #%u: Offset: %" PRIu64 " MBases, Length %" PRIu64 " Mbases\n", i, pa->Offset / 1000000, pa->Length / 1000000);
											if (pa->Type == artValid && pa->Length >= po.RegionLength) {
												int j = 0;
												
												po.ReferenceSequence = pa->Sequence;
#pragma omp parallel for shared(po, st, numberOfAttempts)	
												for (j = 0; j < (int)(pa->Length - po.RegionLength); j += (int)po.RegionLength) {
													const char *refSeq = pa->Sequence + j;
													PROGRAM_STATISTICS tmpstats;

													ret = _test_with_reads(&po, refSeq, &tmpstats);
													if (ret == ERR_SUCCESS) {
														++numberOfAttempts;
														st.FailureCount += tmpstats.FailureCount;
														st.SuccessCount += tmpstats.SuccessCount;
														st.CannotSucceed += tmpstats.CannotSucceed;
													}
												}

												if (pa->Length == po.RegionLength) {
													const char *refSeq = pa->Sequence;
													PROGRAM_STATISTICS tmpstats;

													ret = _test_with_reads(&po, refSeq, &tmpstats);
													if (ret == ERR_SUCCESS) {
														++numberOfAttempts;
														st.FailureCount += tmpstats.FailureCount;
														st.SuccessCount += tmpstats.SuccessCount;
														st.CannotSucceed += tmpstats.CannotSucceed;
													}
												}
											}

											++pa;
										}

										input_free_regions(regions, regionCount);
									}

									if (!explicitSequence) {
										char *rs = NULL;

										utils_free(po.ReferenceSequence);
										ret = fasta_read_seq(&seqFile, &rs, &refSeqLen);
										po.ReferenceSequence = rs;
									}
								} while (ret == ERR_SUCCESS && !explicitSequence);

								if (ret == ERR_NO_MORE_ENTRIES)
									ret = ERR_SUCCESS;

								if (ret == ERR_SUCCESS)
									printf("Success: (%" PRIu64 "), Failures: (%" PRIu64 "), Not tried: (%" PRIu64 "), Percentage: (%" PRIu64 ")\n", st.SuccessCount, st.FailureCount, st.CannotSucceed, (uint64_t)(st.SuccessCount * 100 / (st.SuccessCount + st.FailureCount + st.CannotSucceed)));

								if (!explicitSequence)
									fasta_free(&seqFile);
							}
						}
					} else if (*po.RefSeqFile != '\0') {
						FASTA_FILE seqFile;

						ret = fasta_load(po.RefSeqFile, &seqFile);
						if (ret == ERR_SUCCESS) {
							size_t refSeqLen = 0;
							char *rs = NULL;

							ret = fasta_read_seq(&seqFile, &rs, &refSeqLen);
							po.ReferenceSequence = rs;
							if (ret == ERR_SUCCESS) {
								size_t regionCount = 0;
								PACTIVE_REGION regions = NULL;

								ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
								if (ret == ERR_SUCCESS) {
									size_t index = 0;
									uint64_t regionOffset = 0;

									ret = input_get_region_by_offset(regions, regionCount, po.RegionStart, &index, &regionOffset);
									if (ret == ERR_SUCCESS) {
										PACTIVE_REGION r = regions + index;

										if (r->Type == artValid) {
											ASSEMBLY_TASK task;
											
											po.ReferenceSequence = r->Sequence + regionOffset;
											if (r->Length - regionOffset < po.RegionLength)
												po.RegionLength = r->Length - regionOffset;

											printf("kmer size: %u\n", po.KMerSize);
											printf("Active region (%" PRIu64 "; %u; %u)...\n", po.RegionStart, po.RegionLength, index);
											assembly_task_init(&task, po.ReferenceSequence, po.RegionLength, po.ReferenceSequence, po.RegionLength, po.ReferenceSequence, po.RegionLength, po.Reads, po.ReadCount);
											_compute_graph(&po, &task, &st);
											assembly_task_finit(&task);
										} else printf("ERROR: The active region (%" PRIu64 "; %u; %u) does not specify a readable part of the reference sequence\n", po.RegionStart, po.RegionLength, index);
									}

									input_free_regions(regions, regionCount);
								}
							}

							fasta_free(&seqFile);
						}
					}
				}
			}
		}
	
		options_module_finit();
	}

	return ret;
}
