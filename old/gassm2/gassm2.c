
#include <omp.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "tinydir.h"
#include "err.h"
#include "utils.h"
#include "utils-lookaside.h"
#include "file-utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "reads.h"
#include "pointer_array.h"
#include "gassm2.h"




static PUTILS_LOOKASIDE *_vertexLAs;
static PUTILS_LOOKASIDE *_edgeLAs;



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
		ret = option_add_String(PROGRAM_OPTION_TESTFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_OUTPUT_DIRECTORY, ".");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_VCFFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_Int32(PROGRAM_OPTION_OMP_THREADS, omp_get_num_procs());

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt8(PROGRAM_OPTION_READ_POS_QUALITY, 20);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_CONNECT_REFSEQ, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_CONNECT_READS, FALSE);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_BUBBLE_MERGING, FALSE);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_READ_FIXING, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_LINEAR_SHRINK, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_HELPER_VERTICES, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_MISSING_EDGE_PENALTY, 3);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_BACKWARD_REFSEQ_PENALTY, 2);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_MAX_PATHS, 10);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_READ_MAX_ERROR_RATE, 20);

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
	option_set_shortcut(PROGRAM_OPTION_TESTFILE, 'g');
	option_set_shortcut(PROGRAM_OPTION_OUTPUT_DIRECTORY, 'o');
	option_set_shortcut(PROGRAM_OPTION_VCFFILE, 'v');

	return ret;
}


static ERR_VALUE _capture_program_options(PPROGRAM_OPTIONS Options)
{
	boolean b = FALSE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Options, 0, sizeof(PROGRAM_OPTIONS));
	ret = option_get_UInt32(PROGRAM_OPTION_KMERSIZE, &Options->KMerSize);
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt64(PROGRAM_OPTION_SEQSTART, &Options->RegionStart);

	if (ret == ERR_SUCCESS)
		ret = option_get_Int32(PROGRAM_OPTION_OMP_THREADS, &Options->OMPThreads);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt8(PROGRAM_OPTION_READ_POS_QUALITY, &Options->ReadPosQuality);

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
			
			fprintf(stderr, "Loading reads from %s...\n", readFile);
			ret = input_get_reads(readFile, "sam", &readSet, &readCount);
			if (ret == ERR_SUCCESS) {
				fprintf(stderr, "Filtering out reads with MAPQ less than %u...\n", Options->ReadPosQuality);
				ret = input_filter_bad_reads(readSet, readCount, Options->ReadPosQuality, &Options->Reads, &Options->ReadCount);
				if (ret == ERR_SUCCESS) {
					fprintf(stderr, "Sorting reads...\n");
					input_sort_reads(Options->Reads, Options->ReadCount);
				}

				read_set_destroy(readSet, readCount);
			}
		}
	}

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_ALT1_SEQ, &Options->AltenrateSequence1);

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_ALT2_SEQ, &Options->AlternateSequence2);

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_TESTFILE, &Options->TestFile);

	option_get_Boolean(PROGRAM_OPTION_NO_CONNECT_REFSEQ, &b);
	Options->ParseOptions.ConnectRefSeq = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_CONNECT_READS, &b);
	Options->ParseOptions.ConnectReads = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_BUBBLE_MERGING, &b);
	Options->ParseOptions.MergeBubbles = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_READ_FIXING, &b);
	Options->ParseOptions.FixReads = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_LINEAR_SHRINK, &b);
	Options->ParseOptions.LinearShrink = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_HELPER_VERTICES, &b);
	Options->ParseOptions.HelperVertices = !b;
	option_get_UInt32(PROGRAM_OPTION_MISSING_EDGE_PENALTY, &Options->ParseOptions.MissingEdgePenalty);
	option_get_UInt32(PROGRAM_OPTION_BACKWARD_REFSEQ_PENALTY, &Options->ParseOptions.BackwardRefseqPenalty);
	option_get_UInt32(PROGRAM_OPTION_MAX_PATHS, &Options->MaxPaths);
	option_get_UInt32(PROGRAM_OPTION_READ_MAX_ERROR_RATE, &Options->ParseOptions.ReadMaxErrorRate);

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



static ERR_VALUE _process_variant_call(const ASSEMBLY_TASK *Task, const size_t RefSeqStart, const size_t RefSeqEnd, const char *AltSeq, const size_t AltSeqLen, const size_t RSWeight, const size_t ReadWeight, PGEN_ARRAY_VARIANT_CALL VCArray)
{
	VARIANT_CALL vc;
	char *altSeqStart = NULL;
	size_t rsPos = Task->RegionStart + RefSeqStart + 1;
	size_t rsLen = RefSeqEnd - RefSeqStart;
	size_t altLen = AltSeqLen;
	char *altSeq = NULL;
	const char *refSeq = Task->Reference + RefSeqStart;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(altLen + 2, sizeof(char), &altSeq);
	if (ret == ERR_SUCCESS) {
		altSeqStart = altSeq;
		memcpy(altSeq + 1, AltSeq, altLen*sizeof(char));
		*altSeq = *(refSeq - 1);
		++altSeq;
		altSeq[altLen] = '\0';

		char *opString = NULL;
		size_t opStringLen = 0;
		ret = ssw_clever(refSeq, rsLen, altSeq, altLen, 2, -1, -1, &opString, &opStringLen);;
		if (ret == ERR_SUCCESS) {
			const char *opIt = opString;
			const char *tmpRS = refSeq;
			const char *tmpAltS = altSeq;
			size_t pos = rsPos;
			boolean nothing = TRUE;

			while (ret == ERR_SUCCESS) {
				switch (*opIt) {
					case '\0':
					case 'M':
						if (!nothing) {
							if (altSeq == tmpAltS || refSeq == tmpRS) {
								--rsPos;
								--refSeq;
								--altSeq;
							}

							ret = variant_call_init("1", rsPos, ".", refSeq, tmpRS - refSeq, altSeq, tmpAltS - altSeq, 60, &vc);
							if (ret == ERR_SUCCESS) {
								vc.RefWeight = RSWeight;
								vc.AltWeight = ReadWeight;
								ret = vc_array_add(VCArray, &vc);
								if (ret == ERR_SUCCESS) {
								}

								if (ret != ERR_SUCCESS) {
									variant_call_finit(&vc);
									if (ret == ERR_ALREADY_EXISTS)
										ret = ERR_SUCCESS;
								}

								rsPos += (tmpRS - refSeq);
								refSeq = tmpRS;
								altSeq = tmpAltS;
							}

							nothing = TRUE;
						} else {
							rsPos++;
							refSeq++;
							altSeq++;
						}

						++tmpRS;
						++tmpAltS;
						break;
					case 'X':
						++tmpRS;
						++tmpAltS;
						nothing = FALSE;
						break;
					case 'I':
						++tmpAltS;
						nothing = FALSE;
						break;
					case 'D':
						++tmpRS;
						nothing = FALSE;
						break;
				}

				if (*opIt == '\0')
					break;

				++opIt;
			}
			
			utils_free(opString);
		}

		utils_free(altSeqStart);
	}

	return ret;
}


static ERR_VALUE _process_variant_calls(PGEN_ARRAY_VARIANT_CALL VCArray, const ASSEMBLY_TASK *Task, const GEN_ARRAY_FOUND_SEQUENCE_VARIANT *VariantArray, const size_t Threshold)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t variantCount = gen_array_size(VariantArray);
	const FOUND_SEQUENCE_VARIANT *var = VariantArray->Data;
	const size_t realThreshold = Threshold * 100;

	ret = ERR_SUCCESS;
	for (size_t i = 0; i < variantCount; ++i) {
		size_t rsPos = var->RefSeqStart + Task->RegionStart + 1;
		size_t rsLen = var->RefSeqEnd - var->RefSeqStart;
		char *altSeq = NULL;
		const char *refSeq = Task->Reference + var->RefSeqStart;

		if (var->RefSeqStart < var->RefSeqEnd) {
			if (var->Seq1Weight > realThreshold && var->Seq1Type == kmetRead)
				ret = _process_variant_call(Task, var->RefSeqStart, var->RefSeqEnd, var->Seq1, var->Seq1Len, var->Seq2Weight, var->Seq1Weight, VCArray);

			if (var->Seq2Weight > realThreshold && var->Seq2Type == kmetRead)
				ret = _process_variant_call(Task, var->RefSeqStart, var->RefSeqEnd, var->Seq2, var->Seq2Len, var->Seq1Weight, var->Seq2Weight, VCArray);
		} else {
			printf("VAR-BACK: %u->%u, %s\n", var->RefSeqStart, var->RefSeqEnd, var->Seq1);
		}

		++var;
	}

	return ret;
}


static EExperimentResult _compare_alternate_sequences(const PROGRAM_OPTIONS *Options, PKMER_GRAPH Graph, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	boolean notFound = FALSE;
	POINTER_ARRAY_FOUND_SEQUENCE seqArray;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *alternates[2];
	EExperimentResult res = erNotTried;
	size_t alternateLens[2];
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT variants;

	dym_array_init_FOUND_SEQUENCE_VARIANT(&variants, 140);
	ret = kmer_graph_get_variants(Graph, &variants);
	if (ret == ERR_SUCCESS) {
		alternates[0] = Task->Alternate1;
		alternates[1] = Task->Alternate2;
		alternateLens[0] = Task->Alternate1Length;
		alternateLens[1] = Task->Alternate2Length;
		pointer_array_init_FOUND_SEQUENCE(&seqArray, 140);
		ret = kmer_graph_get_seqs(Graph, Task->Reference, Options->MaxPaths, &seqArray);
		if (ret == ERR_SUCCESS) {
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
			}
			else {
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
			if (Options->VCFFileHandle != NULL)
			{
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

			if (Options->VCFFileHandle != NULL) {
				PFOUND_SEQUENCE *pseq = seqArray.Data;
				PGEN_ARRAY_VARIANT_CALL vcArray = Options->VCSubArrays + omp_get_thread_num();

				_process_variant_calls(vcArray, Task, &variants, Options->Threshold);
				if (pointer_array_size(&seqArray) <= Options->MaxPaths) {
					for (size_t i = 0; i < pointer_array_size(&seqArray); ++i) {
						_process_variant_calls(vcArray, Task, &(*pseq)->ReadVariants, Options->Threshold);
						++pseq;
					}
				}
			}
		}

		for (size_t j = 0; j < pointer_array_size(&seqArray); ++j)
			found_sequence_free(*pointer_array_item_FOUND_SEQUENCE(&seqArray, j));

		pointer_array_finit_FOUND_SEQUENCE(&seqArray);
	}

	dym_array_finit_FOUND_SEQUENCE_VARIANT(&variants);

	return res;
}


static void _on_delete_edge(const KMER_GRAPH *Graph, const KMER_EDGE *Edge, void *Context)
{
	size_t i = 0;
	PGEN_ARRAY_KMER_EDGE_PAIR pairs = (PGEN_ARRAY_KMER_EDGE_PAIR)Context;
	PKMER_EDGE_PAIR p = pairs->Data;

	while (i < gen_array_size(pairs)) {
		if (p->U == Edge || p->V == Edge) {
			if (p->Edges != NULL) {
				utils_free(p->Edges);
				p->Edges = NULL;
			}

			dym_array_remove_fastKMER_EDGE_PAIR(pairs, i);
			continue;
		}

		++p;
		++i;
	}

	return;
}

static EExperimentResult _compute_graph(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const PARSE_OPTIONS *ParseOptions, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	EExperimentResult res = erNotTried;
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = kmer_graph_create(Options->KMerSize, 2500, 6000, &g);	
	if (ret == ERR_SUCCESS) {
		if (Allocator != NULL)
			g->Allocator = *Allocator;
			ret = kmer_graph_parse_ref_sequence(g, Task->Reference, Task->ReferenceLength);
			if (ret == ERR_SUCCESS) {
				GEN_ARRAY_KMER_EDGE_PAIR ep;
				
				dym_array_init_KMER_EDGE_PAIR(&ep, 140);
				ret = kmer_graph_parse_reads(g, Task->Reads, Task->ReadCount, Options->Threshold, ParseOptions, &ep);
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
							if (ParseOptions->LinearShrink)
								kmer_graph_delete_1to1_vertices(g);
								
							if (ParseOptions->MergeBubbles) {
								boolean changed = FALSE;
								
								do {
									changed = FALSE;
									ret = kmer_graph_detect_uncertainities(g, Task->Reference, &changed);
								} while (ret == ERR_SUCCESS && changed);
							}

							if (ret == ERR_SUCCESS)
								res = _compare_alternate_sequences(Options, g, Task, Statistics);
							else printf("ERROR: kmer_graph_detect_uncertainities(): %u\n", ret);
						} else printf("kmer_graph_connect_reads_by_pairs(): %u\n", ret);
					} else ++Statistics->CannotSucceed;
				} else printf("kmer_graph_parse_reads(): %u\n", ret);
			
				PKMER_EDGE_PAIR p = ep.Data;

				for (size_t i = 0; i < gen_array_size(&ep); ++i) {
					if (p->Edges != NULL)
						utils_free(p->Edges);

					++p;
				}

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
										for (size_t i = 0; i < Options->ReadCount; ++i)
											read_split(finalReadSet + i);
										
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
											switch (_compute_graph(NULL, Options, &Options->ParseOptions, &task, &stats)) {
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
		_compute_graph(NULL, Options, &Options->ParseOptions, &task, Statistics);
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



KHASH_MAP_INIT_STR(kc, size_t)

ERR_VALUE kmer_freq_distribution(const PROGRAM_OPTIONS *Options, const uint32_t KMerSize, const ONE_READ *Reads, const size_t ReadCount, const char *FileName)
{
	int err;
	size_t maxValue = 0;
	khiter_t it;
	size_t kmerCount = 0;
	char *kmerString = NULL;
	khash_t(kc) *table = kh_init(kc);
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(KMerSize + 1, sizeof(char), &kmerString);
	if (ret == ERR_SUCCESS) {
		const ONE_READ *r = Reads;
		
		kmerString[KMerSize] = '\0';
		for (size_t i = 0; i < ReadCount; ++i) {
			const READ_PART *p = &r->Part;
			
			if (r->NumberOfFixes * 100 / r->ReadSequenceLen < Options->ParseOptions.ReadMaxErrorRate) {
				read_split(r);
				if (p->ReadSequenceLength >= KMerSize) {
					for (size_t j = 0; j < p->ReadSequenceLength - KMerSize + 1; ++j) {
						char *s = NULL;

						memcpy(kmerString, p->ReadSequence + j, KMerSize*sizeof(char));
						ret = utils_copy_string(kmerString, &s);
						if (ret == ERR_SUCCESS) {
							it = kh_put(kc, table, s, &err);
							switch (err) {
							case 0:
								kh_value(table, it) += 1;
								if (kh_value(table, it) > maxValue)
									maxValue = kh_value(table, it);

								utils_free(s);
								break;
							case 1:
							case 2:
								kh_value(table, it) = 1;
								break;
							default:
								ret = ERR_OUT_OF_MEMORY;
								break;
							}

							++kmerCount;
							if (ret != ERR_SUCCESS)
								utils_free(s);
						}

						if (ret != ERR_SUCCESS)
							break;
					}
				}
			}

			if (ret != ERR_SUCCESS)
				break;

			++r;
		}

		if (ret == ERR_SUCCESS) {
			size_t *freqArray = NULL;

			++maxValue;
			ret = utils_calloc(maxValue, sizeof(size_t), &freqArray);
			if (ret == ERR_SUCCESS) {
				memset(freqArray, 0, maxValue*sizeof(size_t));
				for (it = kh_begin(table); it != kh_end(table); ++it) {
					if (kh_exist(table, it))
						++freqArray[kh_value(table, it)];
				}

				FILE *f = stdout;
				if (FileName != NULL && *FileName != '\0') {
					unlink(FileName);
					ret = utils_fopen(FileName, FOPEN_MODE_WRITE, &f);
				}

				if (ret == ERR_SUCCESS) {
					fprintf(f, "# Total: %Iu, unique: %u (%lf %%)\n", kmerCount, table->n_occupied, 100*(double)table->n_occupied / (double)kmerCount);
					for (size_t i = 0; i < maxValue; ++i) {
						if (freqArray[i] > 0)
							fprintf(f, "%Iu, %Iu, %lf\n", i, freqArray[i], (double)freqArray[i]*100/ (double)kmerCount);
					}

					if (f != stdout)
						utils_fclose(f);
				}

				utils_free(freqArray);
			}
		}

		utils_free(kmerString);
	}

	int i = 0;
#pragma omp parallel for shared(table)
	for (i = (int)kh_begin(table); i < (int)kh_end(table); ++i) {
		if (kh_exist(table, i))
			utils_free(kh_key(table, i));
	}

	kh_destroy(kc, table);

	return ret;
}


typedef struct _BQ_PAIR {
	char Base;
	uint8_t Quality;
	char *Pointer;
	size_t ReadIndex;
} BQ_PAIR, *PBQ_PAIR;

GEN_ARRAY_TYPEDEF(BQ_PAIR);
GEN_ARRAY_IMPLEMENTATION(BQ_PAIR)

typedef struct _JUNK_BASE {
	uint32_t ReadIndex;
	uint8_t Quality;
	char Base;
	char *Pointer;
	/** Zero-based */
	uint64_t Pos;
	GEN_ARRAY_BQ_PAIR Pairs;
} JUNK_BASE, *PJUNK_BASE;

KHASH_MAP_INIT_INT64(jb, JUNK_BASE)


ERR_VALUE fix_reads(PONE_READ Reads, const uint32_t ReadCount, const uint8_t QualityThreshold, const uint32_t CountMultiplier)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	khiter_t it;
	PONE_READ r = Reads;
	khash_t(jb) *table = kh_init(jb);

	ret = ERR_SUCCESS;
	for (uint32_t i = 0; i < ReadCount; ++i) {
		const READ_PART *part = NULL;
/*
		if (r->PosQuality < 60) {
			++r;
			continue;
		}
*/
			read_split(r);
			part = &r->Part;
			for (size_t k = 0; k < part->ReadSequenceLength; ++k) {
				if (part->Quality[k] <= QualityThreshold) {
					JUNK_BASE junkBase;
					int err;

					junkBase.Base = part->ReadSequence[k];
					junkBase.Pos = part->Position + k;
					junkBase.Quality = part->Quality[k];
					junkBase.ReadIndex = i;
					junkBase.Pointer = part->ReadSequence + k;
					it = kh_put(jb, table, junkBase.Pos, &err);
					switch (err) {
					case 0:
						break;
					case 1:
					case 2:
						kh_value(table, it) = junkBase;
						break;
					case -1:
						ret = ERR_OUT_OF_MEMORY;
						break;
					default:
						printf("Error %i\n", err);
						fflush(stdout);
						break;
					}
				}

				if (ret != ERR_SUCCESS)
					break;
			}

			if (ret != ERR_SUCCESS)
				break;

		++r;
	}

	if (ret == ERR_SUCCESS) {
		for (it = kh_begin(table); it != kh_end(table); ++it) {
			if (kh_exist(table, it)) {
				PJUNK_BASE junkBase = &kh_value(table, it);

				dym_array_init_BQ_PAIR(&junkBase->Pairs, 140);
			}
		}

		r = Reads;
		for (uint32_t i = 0; i < ReadCount; ++i) {
			PREAD_PART part = NULL;

			part = &r->Part;
			uint64_t pos = part->Position;

			for (size_t k = 0; k < part->ReadSequenceLength; ++k) {
				/*
				if (r->PosQuality < 60) {
					++r;
					continue;
				}
				*/
				it = kh_get(jb, table, pos);
				if (kh_exist(table, it)) {
					PJUNK_BASE junkBase = &kh_value(table, it);

					if (junkBase->ReadIndex != i) {
						BQ_PAIR p;

						p.Base = part->ReadSequence[k];
						p.Quality = part->Quality[k];
						p.Pointer = part->ReadSequence + k;
						p.ReadIndex = i;
						ret = dym_array_push_back_BQ_PAIR(&junkBase->Pairs, p);
						if (ret != ERR_SUCCESS)
							break;
					}
				}

				++pos;
			}

			if (ret != ERR_SUCCESS)
				break;

			++r;
		}

		for (it = kh_begin(table); it != kh_end(table); ++it) {
			if (kh_exist(table, it)) {
				PJUNK_BASE junkBase = &kh_value(table, it);
				PBQ_PAIR p = junkBase->Pairs.Data;
				uint32_t totalQualities[] = { 0, 0, 0, 0 };
				char bases[] = { 'A', 'C', 'G', 'T' };
				uint32_t counts[] = { 0, 0, 0, 0 };

				for (size_t k = 0; k < gen_array_size(&junkBase->Pairs); ++k) {
					if (p->Quality > QualityThreshold) {
						switch (p->Base) {
							case 'A': totalQualities[0] += p->Quality; ++counts[0]; break;
							case 'C': totalQualities[1] += p->Quality; ++counts[1]; break;
							case 'G': totalQualities[2] += p->Quality; ++counts[2]; break;
							case 'T': totalQualities[3] += p->Quality; ++counts[3]; break;
						}
					}

					++p;
				}

				uint32_t maxIndex = 0;
				uint32_t maxValue = counts[0];
				for (size_t j = 1; j < sizeof(counts) / sizeof(counts[0]); ++j) {
					if (counts[j] > maxValue) {
						maxValue = counts[j];
						maxIndex = j;
					}
				}

				boolean canRepair = TRUE;
				for (size_t j = 0; j < sizeof(counts) / sizeof(counts[0]); ++j) {
					if (maxIndex == j)
						continue;

					canRepair &= (counts[j] * CountMultiplier <= maxValue);
				}

				if (canRepair) {
					char b = bases[maxIndex];
					if (junkBase->Base != b)
						*junkBase->Pointer = b;

					PBQ_PAIR p = junkBase->Pairs.Data;
					for (size_t k = 0; k < gen_array_size(&junkBase->Pairs); ++k) {
						if (p->Base != b)
							*p->Pointer = b;

						++p;
					}
				}

				dym_array_finit_BQ_PAIR(&junkBase->Pairs);
			}
		}
	}

	kh_destroy(jb, table);

	return ret;
}


static double _readBaseCount = 0;
static size_t _totalRegionLength = 0;
static omp_lock_t _readCoverageLock;


ERR_VALUE process_active_region(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const uint64_t RegionStart, const char *RefSeq, PGEN_ARRAY_ONE_READ FilteredReads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	ret = input_filter_reads(Options->KMerSize, Options->Reads, Options->ReadCount, RegionStart, Options->RegionLength, Options->ParseOptions.ReadMaxErrorRate, FilteredReads);
	if (ret == ERR_SUCCESS) {
		if (gen_array_size(FilteredReads) > 0) {
			char taskName[128];
			ASSEMBLY_TASK task;
			PROGRAM_STATISTICS tmpstats;
			
			size_t coverage = 0;
			{
				const ONE_READ *fr = FilteredReads->Data;
				size_t baseCount = 0;

				omp_set_lock(&_readCoverageLock);
				_totalRegionLength += Options->RegionLength;
				for (size_t i = 0; i < gen_array_size(FilteredReads); ++i) {
					_readBaseCount += fr->Part.ReadSequenceLength;
					baseCount += fr->Part.ReadSequenceLength;
					++fr;
				}

				omp_unset_lock(&_readCoverageLock);
				coverage = baseCount / Options->RegionLength;
			}
			
			sprintf(taskName, "%08" PRIu64 " r%Iu", (uint64_t)RegionStart, gen_array_size(FilteredReads));
			assembly_task_init(&task, RefSeq, Options->RegionLength, NULL, 0, NULL, 0, FilteredReads->Data, gen_array_size(FilteredReads));
			assembly_task_set_name(&task, taskName);
			task.RegionStart = RegionStart;

			PARSE_OPTIONS po = Options->ParseOptions;
//			po.ReadThreshold = (coverage > Options->Threshold * 3) ? coverage / 3 : Options->Threshold;
			po.ReadThreshold = Options->Threshold;
			ret = _compute_graph(Allocator, Options, &po, &task, &tmpstats);
			assembly_task_finit(&task);
		}

		input_back_reads(FilteredReads);
	}

	dym_array_clear_ONE_READ(FilteredReads);

	return ret;
}


static ERR_VALUE process_repair_reads(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const uint64_t RegionStart, const char *RefSeq, PGEN_ARRAY_ONE_READ FilteredReads)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = input_filter_reads(Options->KMerSize, Options->Reads, Options->ReadCount, RegionStart, Options->RegionLength, Options->ParseOptions.ReadMaxErrorRate, FilteredReads);
	if (ret == ERR_SUCCESS) {
		if (gen_array_size(FilteredReads) > 0) {
			size_t coverage = 0;
			const ONE_READ *fr = FilteredReads->Data;
			size_t baseCount = 0;

			for (size_t i = 0; i < gen_array_size(FilteredReads); ++i) {
				baseCount += fr->Part.ReadSequenceLength;
				++fr;
			}

			coverage = baseCount / Options->RegionLength;

			PARSE_OPTIONS po = Options->ParseOptions;

			po.ReadThreshold = (coverage >= Options->Threshold * 4) ? coverage / 4 : Options->Threshold;
			ret = assembly_repair_reads(Allocator, Options->KMerSize, FilteredReads->Data, gen_array_size(FilteredReads), RefSeq, Options->RegionLength, &po);
		}

		input_back_reads(FilteredReads);
	}

	dym_array_clear_ONE_READ(FilteredReads);

	return ret;
}


static long _activeRegionCount = 0;
static volatile long _activeRegionProcessed = 0;


static ERR_VALUE _init_lookasides(const uint32_t KmerSize, PUTILS_LOOKASIDE *VA, PUTILS_LOOKASIDE *EA)
{
	ERR_VALUE ret = ERR_SUCCESS;
	int threadIndex = omp_get_thread_num();

	if (_vertexLAs[threadIndex] == NULL) {
		ret = utils_malloc(sizeof(UTILS_LOOKASIDE), _vertexLAs + threadIndex);
		if (ret == ERR_SUCCESS) {
			ret = utils_lookaside_init(_vertexLAs[threadIndex], sizeof(KMER_VERTEX) + KmerSize, 3000);
			if (ret != ERR_SUCCESS) {
				utils_free(_vertexLAs[threadIndex]);
				_vertexLAs[threadIndex] = NULL;
			}
		}
	}

	if (ret == ERR_SUCCESS) {
		*VA = _vertexLAs[threadIndex];
		if (_edgeLAs[threadIndex] == NULL) {
			ret = utils_malloc(sizeof(UTILS_LOOKASIDE), _edgeLAs + threadIndex);
			if (ret == ERR_SUCCESS) {
				ret = utils_lookaside_init(_edgeLAs[threadIndex], sizeof(KMER_EDGE), 5000);
				if (ret != ERR_SUCCESS) {
					utils_free(_edgeLAs[threadIndex]);
					_edgeLAs[threadIndex] = NULL;
				}
			}
		}

		if (ret == ERR_SUCCESS)
			*EA = _edgeLAs[threadIndex];
	}

	return ret;
}


static PKMER_VERTEX _lookaside_vertex_alloc(struct _KMER_GRAPH *Graph, void *Context)
{
	PKMER_VERTEX ret = NULL;
	PUTILS_LOOKASIDE ll = (PUTILS_LOOKASIDE)Context;

	utils_lookaside_alloc(ll, &ret);

	return ret;
}

static void _lookaside_vertex_free(struct _KMER_GRAPH *Graph, PKMER_VERTEX Vertex, void *Context)
{
	PUTILS_LOOKASIDE ll = (PUTILS_LOOKASIDE)Context;

	utils_lookaside_free(ll, Vertex);

	return;
}


static PKMER_EDGE _lookaside_edge_alloc(struct _KMER_GRAPH *Graph, void *Context)
{
	PKMER_EDGE ret = NULL;
	PUTILS_LOOKASIDE ll = (PUTILS_LOOKASIDE)Context;

	utils_lookaside_alloc(ll, &ret);

	return ret;
}

static void _lookaside_edge_free(struct _KMER_GRAPH *Graph, PKMER_EDGE Edge, void *Context)
{
	PUTILS_LOOKASIDE ll = (PUTILS_LOOKASIDE)Context;

	utils_lookaside_free(ll, Edge);

	return;
}

static void process_active_region_in_parallel(const ACTIVE_REGION *Contig, const PROGRAM_OPTIONS *Options)
{
	int j = 0;
	long done = 0;
	PUTILS_LOOKASIDE el = NULL;
	PUTILS_LOOKASIDE vl = NULL;

	_init_lookasides(Options->KMerSize, &vl, &el);
#pragma omp parallel for shared(Options, Contig, _activeRegionProcessed), private(vl, el)
	for (j = 0; j < Contig->Length - Options->RegionLength; j += (int)Options->TestStep) {
		int threadIndex = omp_get_thread_num();
		KMER_GRAPH_ALLOCATOR ga;

		_init_lookasides(Options->KMerSize, &vl, &el);
		ga.VertexAllocatorContext = vl;
		ga.VertexAllocator = _lookaside_vertex_alloc;
		ga.VertexFreer = _lookaside_vertex_free;
		ga.EdgeAllocatorContext = el;
		ga.EdgeAllocator = _lookaside_edge_alloc;
		ga.EdgeFreer = _lookaside_edge_free;
		process_active_region(&ga, Options, Contig->Offset + j, Contig->Sequence + j, Options->ReadSubArrays + omp_get_thread_num());
		done = utils_atomic_increment(&_activeRegionProcessed);
		if (done % (_activeRegionCount / 100) == 0)
			fprintf(stderr, "%u %%\r", done*100 / _activeRegionCount);
	}

	KMER_GRAPH_ALLOCATOR ga;

	_init_lookasides(Options->KMerSize, &vl, &el);
	ga.VertexAllocatorContext = vl;
	ga.VertexAllocator = _lookaside_vertex_alloc;
	ga.VertexFreer = _lookaside_vertex_free;
	ga.EdgeAllocatorContext = el;
	ga.EdgeAllocator = _lookaside_edge_alloc;
	ga.EdgeFreer = _lookaside_edge_free;
	process_active_region(&ga, Options, Contig->Offset + Contig->Length - Options->RegionLength, Contig->Sequence + Contig->Length - Options->RegionLength, Options->ReadSubArrays);

	return;
}


static void repair_reads_in_parallel(const ACTIVE_REGION *Contig, const PROGRAM_OPTIONS *Options)
{
	const size_t iterations = 2;
	const long count = _activeRegionCount*iterations;
	const uint32_t realStep = ((Options->RegionLength + Options->TestStep - 1) / Options->TestStep)*Options->TestStep;
	long done;
	int threadIndex = omp_get_thread_num();
	ERR_VALUE ret = ERR_SUCCESS;
	PUTILS_LOOKASIDE el = NULL;
	PUTILS_LOOKASIDE vl = NULL;

	_init_lookasides(Options->KMerSize, &vl, &el);
	for (size_t it = 0; it < iterations; ++it) {
		for (uint32_t k = 0; k < realStep; k += Options->TestStep) {
			int j = 0;
#pragma omp parallel for shared(Options, Contig, k, _activeRegionProcessed), private(vl, el)
			for (j = k; j < Contig->Length - Options->RegionLength; j += (int)realStep) {
				KMER_GRAPH_ALLOCATOR ga;
				
				_init_lookasides(Options->KMerSize, &vl, &el);
				ga.VertexAllocatorContext = vl;
				ga.VertexAllocator = _lookaside_vertex_alloc;
				ga.VertexFreer = _lookaside_vertex_free;
				ga.EdgeAllocatorContext = el;
				ga.EdgeAllocator = _lookaside_edge_alloc;
				ga.EdgeFreer = _lookaside_edge_free;
				process_repair_reads(&ga, Options, Contig->Offset + j, Contig->Sequence + j, Options->ReadSubArrays + omp_get_thread_num());
				done = utils_atomic_increment(&_activeRegionProcessed);
				if (done % (count / 100) == 0)
					fprintf(stderr, "%u %%\r", done * 100 / count);
			}
		}

		KMER_GRAPH_ALLOCATOR ga;

		_init_lookasides(Options->KMerSize, &vl, &el);
		ga.VertexAllocatorContext = vl;
		ga.VertexAllocator = _lookaside_vertex_alloc;
		ga.VertexFreer = _lookaside_vertex_free;
		ga.EdgeAllocatorContext = el;
		ga.EdgeAllocator = _lookaside_edge_alloc;
		ga.EdgeFreer = _lookaside_edge_free;
		for (uint32_t k = 0; k < realStep; k += Options->TestStep) {
			process_repair_reads(&ga, Options, Contig->Offset + Contig->Length - Options->RegionLength - k, Contig->Sequence + Contig->Length - Options->RegionLength - k, Options->ReadSubArrays);
			done = utils_atomic_increment(&_activeRegionProcessed);
			if (done % (count / 100) == 0)
				fprintf(stderr, "%u %%\r", done * 100 / count);
		}
	}

	return;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	omp_init_lock(&_readCoverageLock);
#ifdef _MSC_VER
	uint64_t startTime = GetTickCount64();
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
//#ifdef _DEBUG
//					omp_set_num_threads(1);
//#else
					omp_set_num_threads(po.OMPThreads);
//#endif
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
							_compute_graph(NULL, &po, &po.ParseOptions, &task, &stats);
							assembly_task_finit(&task);
						} else {
							tinydir_dir dir;
							PROGRAM_STATISTICS stats[128];

							memset(&stats, 0, sizeof(stats));
							if (tinydir_open(&dir, po.TestFile) == ERR_SUCCESS) {
								POINTER_ARRAY_char fileNameArray;
								
								pointer_array_init_char(&fileNameArray, 140);
								ret = pointer_array_reserve_char(&fileNameArray, 16384);
								if (ret == ERR_SUCCESS) {
									do {
										_obtain_files(&fileNameArray, 16384, &dir);
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
												_compute_graph(NULL, &po, &po.ParseOptions, &task, stats + omp_get_thread_num());
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

							printf("%" PRIu64 " %" PRIu64 " me=%u,br=%u\n", st.SuccessCount, st.FailureCount, po.ParseOptions.MissingEdgePenalty, po.ParseOptions.BackwardRefseqPenalty);
						}
					} else if (*po.RefSeqFile != '\0') {
						printf("K-mer size:                 %u\n", po.KMerSize);
						printf("Active region length:       %u\n", po.RegionLength);
						printf("Reference:                  %s\n", po.RefSeqFile);
						printf("Reads:                      %u\n", po.ReadCount);
						printf("Read coverage threshold:    %u\n", po.Threshold);
						printf("Min. read position quality: %u\n", po.ReadPosQuality);
						printf("OpenMP thread count:        %i\n", po.OMPThreads);
						printf("Output VCF file:            %s\n", po.VCFFile);
						ret = paired_reads_init();
						if (ret == ERR_SUCCESS) {
//							printf("Computing k-mer frequency distribution...\n");
//							kmer_freq_distribution(po.KMerSize, po.Reads, po.ReadCount, "kmer-dist1.csv");
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
										ret = utils_calloc(omp_get_num_procs(), sizeof(PUTILS_LOOKASIDE), &_vertexLAs);
										if (ret == ERR_SUCCESS)
											ret = utils_calloc(omp_get_num_procs(), sizeof(PUTILS_LOOKASIDE), &_edgeLAs);
										
										ret = utils_calloc(omp_get_num_procs(), sizeof(GEN_ARRAY_VARIANT_CALL), &po.VCSubArrays);
										if (ret == ERR_SUCCESS) {
											ret = utils_calloc(omp_get_num_procs(), sizeof(GEN_ARRAY_ONE_READ), &po.ReadSubArrays);
											if (ret == ERR_SUCCESS) {
												const size_t numThreads = omp_get_num_procs();
												for (size_t i = 0; i < numThreads; ++i) {
													dym_array_init_VARIANT_CALL(po.VCSubArrays + i, 140);
													dym_array_init_ONE_READ(po.ReadSubArrays + i, 140);
													_vertexLAs[i] = NULL;
													_edgeLAs[i] = NULL;
												}

												size_t regionCount = 0;
												PACTIVE_REGION regions = NULL;

												ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
												if (ret == ERR_SUCCESS) {
													const ACTIVE_REGION *pa = NULL;

													pa = regions;
													for (size_t i = 0; i < regionCount; ++i) {
														if (pa->Type == artValid && pa->Length >= po.RegionLength)
															_activeRegionCount += (pa->Length / po.TestStep);

														++pa;
													}

													if (po.ParseOptions.FixReads) {
														printf("Repairing reads...\n");
														_activeRegionProcessed = 0;
														pa = regions;
														for (size_t i = 0; i < regionCount; ++i) {
															if (pa->Type == artValid && pa->Length >= po.RegionLength)
																repair_reads_in_parallel(pa, &po);

															++pa;
														}
													}
														
													printf("Calling variants...\n");
													_activeRegionProcessed = 0;
													pa = regions;
													for (size_t i = 0; i < regionCount; ++i) {
														if (pa->Type == artValid && pa->Length >= po.RegionLength)
															process_active_region_in_parallel(pa, &po);
														
														++pa;
													}
														
													input_free_regions(regions, regionCount);
												}

												utils_free(rsFasta);
												ret = vc_array_merge(&po.VCArray, po.VCSubArrays, numThreads);
												for (size_t i = 0; i < numThreads; ++i) {
													dym_array_finit_ONE_READ(po.ReadSubArrays + i);
													vc_array_finit(po.VCSubArrays + i);
												}

												utils_free(po.ReadSubArrays);
											}

											utils_free(po.VCSubArrays);
										}

										utils_free(_edgeLAs);
										utils_free(_vertexLAs);

										if (po.VCFFileHandle != NULL) {
											if (ret == ERR_SUCCESS)
												vc_array_print(po.VCFFileHandle, &po.VCArray);

											vc_array_finit(&po.VCArray);
											utils_fclose(po.VCFFileHandle);
										}

									}

									fasta_free(&seqFile);
								}
							} else printf("fix_reads(): %u\n", ret);

							printf("Computing k-mer frequency distribution...\n");
							printf("Read coverage: %lf\n", _readBaseCount / _totalRegionLength );
							kmer_freq_distribution(&po, po.KMerSize, po.Reads, po.ReadCount, "kmer-dist2.csv");

							paired_reads_finit();
						}
					}
				}
			}
		}
	
		options_module_finit();
	}

#ifdef _MSC_VER
	uint64_t endTime = GetTickCount64();
	printf("Time: %I64u s\n", (endTime - startTime) / 1000);
#endif
	omp_destroy_lock(&_readCoverageLock);

	return ret;
}
