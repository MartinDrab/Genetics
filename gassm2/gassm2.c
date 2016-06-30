
#include <omp.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "tinydir.h"
#include "err.h"
#include "utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "reads.h"
#include "pointer_array.h"
#include "ssw.h"
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

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_TESTFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_OUTPUT_DIRECTORY, ".");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_VCFFILE, "\0");

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
	option_set_description_const(PROGRAM_OPTION_TESTFILE, PROGRAM_OPTION_TESTFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_OUTPUT_DIRECTORY, PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC);
	option_set_description_const(PROGRAM_OPTION_VCFFILE, PROGRAM_OPTION_VCFFILE_DESC);

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
	option_set_shortcut(PROGRAM_OPTION_TESTFILE, 'g');
	option_set_shortcut(PROGRAM_OPTION_OUTPUT_DIRECTORY, 'o');
	option_set_shortcut(PROGRAM_OPTION_VCFFILE, 'v');

	return ret;
}


static ERR_VALUE _capture_program_options(PPROGRAM_OPTIONS Options)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Options, 0, sizeof(PROGRAM_OPTIONS));
	ret = option_get_UInt32(PROGRAM_OPTION_KMERSIZE, &Options->KMerSize);
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt64(PROGRAM_OPTION_SEQSTART, &Options->RegionStart);

	if (ret == ERR_SUCCESS) {
		char *outputDirectory = NULL;

		ret = option_get_String(PROGRAM_OPTION_OUTPUT_DIRECTORY, &outputDirectory);
		if (ret == ERR_SUCCESS) {
			size_t len = strlen(outputDirectory);

			if (outputDirectory[len] == '\\')
				outputDirectory[len] = '\0';

			Options->OutputDirectoryBase = outputDirectory;
		}
	}

	if (ret == ERR_SUCCESS) {
		char *vcfFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_VCFFILE, &vcfFile);
		if (ret == ERR_SUCCESS)
			Options->VCFFile = vcfFile;
	}

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
		char *readFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_READFILE, &readFile);
		if (ret == ERR_SUCCESS && *readFile != '\0') {
			PONE_READ readSet = NULL;
			size_t readCount = 0;
			
			ret = input_get_reads(readFile, "sam", &readSet, &readCount);
			if (ret == ERR_SUCCESS) {
				ret = input_filter_bad_reads(readSet, readCount, 20, &Options->Reads, &Options->ReadCount);
				read_set_destroy(readSet, readCount);
			}
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

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_TESTFILE, &Options->TestFile);

	return ret;
}


static void _write_differences(FILE *Stream, const char *RefSeq, const char *Alt1, const char *Alt2, const size_t Length)
{
	for (size_t i = 0; i < Length; ++i) {

		if (*RefSeq != *Alt1 || *RefSeq != *Alt2)
			fprintf(Stream, "%Iu:\t%c\t%c\t%c\n", i, *RefSeq, *Alt1, *Alt2);

		++RefSeq;
		++Alt1;
		++Alt2;
	}

	return;
}


typedef enum _EExperimentResult {
	erSuccess,
	erFailure,
	erNotTried,
} EExperimentResult, *PEExperimentResult;


static EExperimentResult _compare_alternate_sequences(const PROGRAM_OPTIONS *Options, const KMER_GRAPH *Graph, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	boolean notFound = FALSE;
	POINTER_ARRAY_FOUND_SEQUENCE seqArray;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *alternates[2];
	EExperimentResult res = erNotTried;
	size_t alternateLens[2];

	alternates[0] = Task->Alternate1;
	alternates[1] = Task->Alternate2;
	alternateLens[0] = Task->Alternate1Length;
	alternateLens[1] = Task->Alternate2Length;
	pointer_array_init_FOUND_SEQUENCE(&seqArray, 140);
	ret = kmer_graph_get_seqs(Graph, &seqArray);
	if (ret == ERR_SUCCESS) {
		kmer_graph_delete_seqs(Graph, &seqArray, Options->Threshold);
		if (Task->Alternate1Length > 0 && Task->Alternate2Length > 0) {
			for (size_t i = 0; i < sizeof(alternateLens) / sizeof(size_t); ++i) {
				boolean found = FALSE;

				for (size_t j = 0; j < gen_array_size(&seqArray); ++j) {
					const FOUND_SEQUENCE *fs = *pointer_array_item_FOUND_SEQUENCE(&seqArray, j);

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
					//				_write_differences(Task->Reference, Task->Alternate1, Task->Alternate2, Task->ReferenceLength);
					//				exit(0);
					break;
				}
			}
		}

		if (ret != ERR_SUCCESS || notFound) {
			++Statistics->FailureCount;
			res = erFailure;
		} else {
			++Statistics->SuccessCount;
			res = erSuccess;
		}
	} else printf("ERROR: kmer_graph_get_seqs(): %u\n", ret);
	
	if (Task->Name != NULL && (Graph->TypedEdgeCount[kmetRead] > 0 ||
		Graph->TypedEdgeCount[kmetVariant] > 0)) {
		FILE *f = NULL;
		char *directory = NULL;
		char graphName[128];
		char diffName[128];
		char vcfName[128];

		switch (res) {
			case erSuccess:
				directory = "succ";
				break;
			case erFailure:
				directory = "fail";
				break;
			case erNotTried:
				directory = "nottried";
				break;
			default:
				assert(FALSE);
				break;
		}

#pragma warning (disable : 4996)											
		sprintf(graphName, "%s" PATH_SEPARATOR "%s" PATH_SEPARATOR "%s-%Iu.graph", Options->OutputDirectoryBase, directory, Task->Name, gen_array_size(&seqArray));
		unlink(graphName);
#pragma warning (disable : 4996)											
		sprintf(diffName, "%s" PATH_SEPARATOR "%s" PATH_SEPARATOR "%s-%Iu.diff", Options->OutputDirectoryBase, directory, Task->Name, gen_array_size(&seqArray));
		unlink(diffName);
#pragma warning (disable : 4996)											
		sprintf(vcfName, "%s" PATH_SEPARATOR "%s" PATH_SEPARATOR "%s-%Iu.vcf", Options->OutputDirectoryBase, directory, Task->Name, gen_array_size(&seqArray));
		unlink(vcfName);
		if (Options->VCFFileHandle != NULL) {
			ret = utils_fopen(graphName, FOPEN_MODE_WRITE, &f);
			if (ret == ERR_SUCCESS) {
				kmer_graph_print(f, Graph);
				utils_fclose(f);
			}

			if (Task->Alternate1Length > 0 && Task->Alternate2Length > 0) {
				ret = utils_fopen(diffName, FOPEN_MODE_WRITE, &f);
				if (ret == ERR_SUCCESS) {
					_write_differences(f, Task->Reference, Task->Alternate1, Task->Alternate2, Task->ReferenceLength);
					utils_fclose(f);
				}
			}
		}

		if (pointer_array_size(&seqArray) == 1 && Options->VCFFileHandle != NULL) {
			PFOUND_SEQUENCE seq = *pointer_array_item_FOUND_SEQUENCE(&seqArray, 0);
			const size_t variantCount = gen_array_size(&seq->Variants);
			const FOUND_SEQUENCE_VARIANT *var = seq->Variants.Data;
			PGEN_ARRAY_VARIANT_CALL vcArray = Options->VCSubArrays + omp_get_thread_num();

			for (size_t i = 0; i < variantCount; ++i) {
				char *opString = NULL;
				size_t opStringLen = 0;

				if (var->Seq1Weight > Options->Threshold && var->Seq1Type == kmetRead) {
					ret = ssw_clever(Task->Reference + var->RefSeqStart, var->RefSeqEnd - var->RefSeqStart, var->Seq1, var->Seq1Len, 2, -1, -1, &opString, &opStringLen);
					if (ret == ERR_SUCCESS) {
						write_seq_differences(vcArray, Task->Reference + var->RefSeqStart, Task->RegionStart + var->RefSeqStart, var->RefSeqEnd - var->RefSeqStart, opString, var->Seq1, var);
						utils_free(opString);
					}
				}

				if (var->Seq2Weight > Options->Threshold && var->Seq2Type == kmetRead) {
					ret = ssw_clever(Task->Reference + var->RefSeqStart, var->RefSeqEnd - var->RefSeqStart, var->Seq2, var->Seq2Len, 2, -1, -1, &opString, &opStringLen);
					if (ret == ERR_SUCCESS) {
						write_seq_differences(vcArray, Task->Reference + var->RefSeqStart, Task->RegionStart + var->RefSeqStart, var->RefSeqEnd - var->RefSeqStart, opString, var->Seq2, var);
						utils_free(opString);
					}
				}

				++var;
			}			
		}
	}
	
	for (size_t j = 0; j < gen_array_size(&seqArray); ++j)
		found_sequence_free(*pointer_array_item_FOUND_SEQUENCE(&seqArray, j));

	pointer_array_finit_FOUND_SEQUENCE(&seqArray);

	return res;
}


static void _on_delete_edge(const KMER_GRAPH *Graph, const KMER_EDGE *Edge, void *Context)
{
	size_t i = 0;
	PGEN_ARRAY_KMER_EDGE_PAIR pairs = (PGEN_ARRAY_KMER_EDGE_PAIR)Context;
	PKMER_EDGE_PAIR p = pairs->Data;

	while (i < gen_array_size(pairs)) {
		if (p->U == Edge || p->V == Edge) {
			dym_array_remove_fastKMER_EDGE_PAIR(pairs, i);
			continue;
		}

		++p;
		++i;
	}

	return;
}

static EExperimentResult _compute_graph(const PROGRAM_OPTIONS *Options, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	EExperimentResult res = erNotTried;
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

		ret = kmer_graph_create(Options->KMerSize, &g);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_parse_ref_sequence(g, Task->Reference, Task->ReferenceLength, Options->Threshold);
			if (ret == ERR_SUCCESS) {
				GEN_ARRAY_KMER_EDGE_PAIR ep;
				
				dym_array_init_KMER_EDGE_PAIR(&ep, 140);
				ret = kmer_graph_parse_reads(g, Task->Reads, Task->ReadCount, Task->RegionStart, &ep);
				if (ret == ERR_SUCCESS) {
					size_t deletedThings = 0;

					g->DeleteEdgeCallback = _on_delete_edge;
					g->DeleteEdgeCallbackContext = &ep;
					kmer_graph_delete_edges_under_threshold(g, Options->Threshold);
					kmer_graph_delete_trailing_things(g, &deletedThings);
					g->DeleteEdgeCallback = NULL;
					if (g->TypedEdgeCount[kmetRead] > 0) {
						size_t changeCount = 0;
						ret = kmer_graph_connect_reads_by_pairs(g, Options->Threshold, &ep, &changeCount);
						if (ret == ERR_SUCCESS) {
							ret = kmer_graph_connect_reads_by_reads(g, Options->Threshold);
							if (ret == ERR_SUCCESS) {
								kmer_graph_delete_1to1_vertices(g);
								kmer_graph_delete_backward_edges(g);
								boolean changed = FALSE;
								do {
									changed = FALSE;
									ret = kmer_graph_detect_uncertainities(g, &changed);
									kmer_graph_delete_trailing_things(g, &deletedThings);
								} while (ret == ERR_SUCCESS && changed);

								if (ret == ERR_SUCCESS) {
									res = _compare_alternate_sequences(Options, g, Task, Statistics);
								} else printf("ERROR: kmer_graph_detect_uncertainities(): %u\n", ret);
							} else printf("kmer_graph_connect_reads(): %u\n", ret);
						} else printf("kmer_graph_connect_reads_by_reads(): %u\n", ret);
					} else {
						++Statistics->CannotSucceed;
						printf("NOTTRIED\n");
					}
				} else printf("kmer_graph_parse_reads(): %u\n", ret);
			
				dym_array_finit_KMER_EDGE_PAIR(&ep);
			} else printf("kmer_graph_parse_ref_sequence(): %u\n", ret);

			kmer_graph_destroy(g);
		} else printf("kmer_graph_create(): %u\n", ret);

		if (ret != ERR_SUCCESS) {
			++Statistics->FailureCount;
			printf("FAILD: %u\n", ret);
		}

	return res;
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


static ERR_VALUE _examine_read_coverage(const ONE_READ *Reads, const size_t ReadCount, const char *RefSeq, const size_t RefSeqLen, const char *Alternate)
{
	uint32_t *c = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(RefSeqLen, sizeof(uint32_t), &c);
	if (ret == ERR_SUCCESS) {
		memset(c, 0, RefSeqLen*sizeof(uint32_t));
		for (size_t i = 0; i < ReadCount; ++i) {
			for (size_t j = 0; j < Reads[i].ReadSequenceLen; ++j)
				c[Reads[i].Pos + j]++;
		}

		printf("Not covered: ");
		for (size_t i = 0; i < RefSeqLen; ++i) {
			if (c[i] == 0) {
				printf("%u ", i);
				if (RefSeq[i] != Alternate[i]) {
					printf("%u: The position has SNPs but is not covered in any read (%c %c)\n", i, RefSeq[i], Alternate[i]);
					ret = ERR_BAD_READ_COVERAGE;
				}
			}
		}

		printf("\n");
		utils_free(c);
	}

	return ret;
}


static unsigned int _taskNumber = 0;

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
					ret = _examine_read_coverage(reads1, Options->ReadCount / 2, RefSeq, Options->RegionLength, alternate);
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
								ret = _examine_read_coverage(reads2, Options->ReadCount / 2, RefSeq, Options->RegionLength, alternate2);
								if (ret == ERR_SUCCESS) {
									PONE_READ finalReadSet = reads1;

									ret = read_set_merge(&finalReadSet, Options->ReadCount / 2, reads2, Options->ReadCount / 2);
									if (ret == ERR_SUCCESS) {
										char *rs = NULL;
										PROGRAM_STATISTICS stats;

										reads1 = NULL;
										reads2 = NULL;
										for (size_t i = 0; i < Options->ReadCount; ++i) {
											boolean tmp = FALSE;

											read_split(finalReadSet + i, 0, 0, &tmp);
										}
										
										ret = utils_calloc(Options->RegionLength + 1, sizeof(char), &rs);
										if (ret == ERR_SUCCESS) {
											char taskFileName[128];
											char succTaskFileName[128];
											char failedTaskFileName[128];
											char notTriedTaskFileName[128];

											memcpy(rs, RefSeq, Options->RegionLength*sizeof(char));
											rs[Options->RegionLength] = '\0';
											assembly_task_init(&task, rs, Options->RegionLength, alternate, alternateLen, alternate2, alternateLen2, finalReadSet, Options->ReadCount);
#pragma warning (disable : 4996)											
											sprintf(taskFileName, "%s" PATH_SEPARATOR "tmp" PATH_SEPARATOR "%09u.task", Options->OutputDirectoryBase, _taskNumber);
#pragma warning (disable : 4996)											
											sprintf(succTaskFileName, "%s" PATH_SEPARATOR "succ" PATH_SEPARATOR "%09u.task", Options->OutputDirectoryBase, _taskNumber);
#pragma warning (disable : 4996)											
											sprintf(notTriedTaskFileName, "%s" PATH_SEPARATOR "nottried" PATH_SEPARATOR "%09u.task", Options->OutputDirectoryBase, _taskNumber);
#pragma warning (disable : 4996)						
											sprintf(failedTaskFileName, "%s" PATH_SEPARATOR "fail" PATH_SEPARATOR "%09u.task", Options->OutputDirectoryBase, _taskNumber);
											++_taskNumber;
											unlink(taskFileName);
											assembly_task_save_file(taskFileName, &task);
											switch (_compute_graph(Options, &task, &stats)) {
												case erSuccess:
													unlink(succTaskFileName);
													rename(taskFileName, succTaskFileName);
													break;
												case erFailure:
													unlink(failedTaskFileName);
													rename(taskFileName, failedTaskFileName);
													break;
												case erNotTried:
													unlink(notTriedTaskFileName);
													rename(taskFileName, notTriedTaskFileName);
													break;
												default:
													assert(FALSE);
													break;
											}

											assembly_task_finit(&task);
											Statistics->FailureCount += stats.FailureCount;
											Statistics->SuccessCount += stats.SuccessCount;
											Statistics->CannotSucceed += stats.CannotSucceed;
											utils_free(rs);
										}

										read_set_destroy(finalReadSet, Options->ReadCount);
									}
								}
								else Statistics->CannotSucceed++;


								if (reads2 != NULL)
									read_set_destroy(reads2, Options->ReadCount / 2);
							}

							if (*Options->AlternateSequence2 == '\0')
								utils_free(alternate2);
						}
					} else Statistics->CannotSucceed++;

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





static ERR_VALUE _obtain_files(PPOINTER_ARRAY_char Array, const size_t MaxCount, tinydir_dir *Dir)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = ERR_SUCCESS;
	pointer_array_clear_char(Array);
	if (Dir->has_next) {
		
		for (size_t i = 0; i < MaxCount; ++i) {
			tinydir_file file;
			char *str = NULL;

			if (ret != ERR_SUCCESS || !Dir->has_next)
				break;

			tinydir_readfile(Dir, &file);
			if (strcmp(file.name, ".") == 0 || strcmp(file.name, "..") == 0) {
				tinydir_next(Dir);
				continue;
			}

			ret = utils_copy_string(file.path, &str);
			if (ret == ERR_SUCCESS)
				pointer_array_push_back_no_alloc_char(Array, str);
			
			if (ret != ERR_SUCCESS || !Dir->has_next)
				break;

			tinydir_next(Dir);
		}
	}

	return ret;
}


ERR_VALUE process_active_region(const PROGRAM_OPTIONS *Options, const uint64_t RegionStart, const char *RefSeq, PGEN_ARRAY_ONE_READ FilteredReads)
{
	boolean indels = FALSE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	ret = input_filter_reads(Options->Reads, Options->ReadCount, RegionStart, Options->RegionLength, &indels, FilteredReads);
	if (ret == ERR_SUCCESS) {
		if (gen_array_size(FilteredReads) > 0 && indels) {
			char taskName[128];
			ASSEMBLY_TASK task;
			PROGRAM_STATISTICS tmpstats;

			sprintf(taskName, "%08" PRIu64 " r%Iu", (uint64_t)RegionStart, gen_array_size(FilteredReads));
			assembly_task_init(&task, RefSeq, Options->RegionLength, NULL, 0, NULL, 0, FilteredReads->Data, gen_array_size(FilteredReads));
			assembly_task_set_name(&task, taskName);
			task.RegionStart = RegionStart;
			ret = _compute_graph(Options, &task, &tmpstats);
			assembly_task_finit(&task);
		}
	}

	dym_array_clear_ONE_READ(FilteredReads);
	
	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	omp_set_num_threads(omp_get_num_procs());
#ifdef _DEBUG
	omp_set_num_threads(1);
#endif
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
							char *rsFasta = NULL;

							if (!explicitSequence) {
								ret = fasta_load(po.RefSeqFile, &seqFile);
								if (ret == ERR_SUCCESS) {
									ret = fasta_read_seq(&seqFile, &rsFasta, &refSeqLen);
									po.ReferenceSequence = rsFasta;
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
										utils_free(rsFasta);
										ret = fasta_read_seq(&seqFile, &rsFasta, &refSeqLen);
										po.ReferenceSequence = rsFasta;
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
					} else if (*po.TestFile != '\0') {
						ASSEMBLY_TASK task;
						ret = assembly_task_load_file(po.TestFile, &task);
						if (ret == ERR_SUCCESS) {							
							PROGRAM_STATISTICS stats;
							char *taskName = strchr(po.TestFile, '/');

							if (taskName == NULL)
								taskName = po.TestFile;
							else taskName++;

							assembly_task_set_name(&task, taskName);
							memset(&stats, 0, sizeof(stats));
							_compute_graph(&po, &task, &stats);
							assembly_task_finit(&task);
						} else {
							tinydir_dir dir;
							PROGRAM_STATISTICS stats[128];

							memset(&stats, 0, sizeof(stats));
							if (tinydir_open(&dir, po.TestFile) == ERR_SUCCESS) {
								printf("The given \"file\" is a directory\n");
								POINTER_ARRAY_char fileNameArray;
								
								pointer_array_init_char(&fileNameArray, 140);
								ret = pointer_array_reserve_char(&fileNameArray, 4096);
								if (ret == ERR_SUCCESS) {
									do {
										_obtain_files(&fileNameArray, 4096, &dir);
										int i = 0;
#pragma omp parallel for shared(po, stats, fileNameArray)	
										for (i = 0; i < (int)pointer_array_size(&fileNameArray); ++i) {
											char *fileName = *pointer_array_item_char(&fileNameArray, i);
											ASSEMBLY_TASK task;

											ret = assembly_task_load_file(fileName, &task);
											if (ret == ERR_SUCCESS) {
												char *taskName = strchr(fileName, '/');

												if (taskName == NULL)
													taskName = fileName;
												else taskName++;

												assembly_task_set_name(&task, taskName);
												task.RegionStart = 0;
												_compute_graph(&po, &task, stats + omp_get_thread_num());
												assembly_task_finit(&task);
											}

											utils_free(fileName);
										}
									} while (dir.has_next);
								}

								pointer_array_finit_char(&fileNameArray);
								tinydir_close(&dir);
							}

							PROGRAM_STATISTICS st;

							memset(&st, 0, sizeof(st));
							for (size_t i = 0; i < sizeof(stats) / sizeof(stats[0]); ++i) {
								st.CannotSucceed += stats[i].CannotSucceed;
								st.FailureCount += stats[i].FailureCount;
								st.SuccessCount += stats[i].SuccessCount;
							}

							printf("Success: (%" PRIu64 "), Failures: (%" PRIu64 "), Not tried: (%" PRIu64 "), Percentage: (%" PRIu64 ")\n", st.SuccessCount, st.FailureCount, st.CannotSucceed, (uint64_t)(st.SuccessCount * 100 / (st.SuccessCount + st.FailureCount + st.CannotSucceed + 1)));
						}
					} else if (*po.RefSeqFile != '\0') {
						ret = paired_reads_init();
						if (ret == ERR_SUCCESS) {
							size_t refSeqLen = 0;
							FASTA_FILE seqFile;
							char *rsFasta = NULL;

							ret = fasta_load(po.RefSeqFile, &seqFile);
							if (ret == ERR_SUCCESS) {
								ret = fasta_read_seq(&seqFile, &rsFasta, &refSeqLen);
								po.ReferenceSequence = rsFasta;
								if (ret != ERR_SUCCESS)
									fasta_free(&seqFile);
							}

							if (ret == ERR_SUCCESS) {
								po.VCFFileHandle = NULL;
								if (*po.VCFFile != '\0') {
									ret = utils_fopen(po.VCFFile, FOPEN_MODE_WRITE, &po.VCFFileHandle);
									if (ret == ERR_SUCCESS)
										dym_array_init_VARIANT_CALL(&po.VCArray, 140);
								}

								if (ret == ERR_SUCCESS) {
									ret = utils_calloc(omp_get_num_procs(), sizeof(GEN_ARRAY_VARIANT_CALL), &po.VCSubArrays);
									if (ret == ERR_SUCCESS) {
										ret = utils_calloc(omp_get_num_procs(), sizeof(GEN_ARRAY_ONE_READ), &po.ReadSubArrays);
										if (ret == ERR_SUCCESS) {
											const size_t numThreads = omp_get_num_procs();
											for (size_t i = 0; i < numThreads; ++i) {
												dym_array_init_VARIANT_CALL(po.VCSubArrays + i, 140);
												dym_array_init_ONE_READ(po.ReadSubArrays + i, 140);
											}

											do {
												size_t regionCount = 0;
												PACTIVE_REGION regions = NULL;

												ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
												if (ret == ERR_SUCCESS) {
													for (size_t i = 0; i < regionCount; ++i) {
														PACTIVE_REGION pa = regions + i;

														if (pa->Type == artValid && pa->Length >= po.RegionLength) {
															int j = 0;
															const int regionLength = po.RegionLength;
															const int intervalLength = pa->Length;
															const int step = po.TestStep;
#pragma omp parallel for shared(po, pa)
															for (j = 0; j < intervalLength - regionLength; j += step)
																process_active_region(&po, pa->Offset + j, pa->Sequence + j, po.ReadSubArrays + omp_get_thread_num());

															process_active_region(&po, pa->Offset + pa->Length - po.RegionLength, pa->Sequence + pa->Length - po.RegionLength, po.ReadSubArrays);
														}

														++pa;
													}

													input_free_regions(regions, regionCount);
												}

												utils_free(rsFasta);
												ret = fasta_read_seq(&seqFile, &rsFasta, &refSeqLen);
												po.ReferenceSequence = rsFasta;
											} while (ret == ERR_SUCCESS);

											if (ret == ERR_NO_MORE_ENTRIES)
												ret = ERR_SUCCESS;

											ret = vc_array_merge(&po.VCArray, po.VCSubArrays, numThreads);
											for (size_t i = 0; i < numThreads; ++i) {
												dym_array_finit_ONE_READ(po.ReadSubArrays + i);
												vc_array_finit(po.VCSubArrays + i);
											}

											utils_free(po.ReadSubArrays);
										}

										utils_free(po.VCSubArrays);
									}

									if (po.VCFFileHandle != NULL) {
										if (ret == ERR_SUCCESS)
											vc_array_print(po.VCFFileHandle, &po.VCArray);

										vc_array_finit(&po.VCArray);
										utils_fclose(po.VCFFileHandle);
									}

								}

								fasta_free(&seqFile);
							}

							paired_reads_finit();
						}
					}
				}
			}
		}
	
		options_module_finit();
	}

	return ret;
}
