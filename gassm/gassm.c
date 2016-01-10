
#include <omp.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include <windows.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "kmer-graph.h"
#include "input-file.h"
#include "assembly.h"
#include "libkmer.h"
#include "reads.h"
#include "ssw.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 15);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_FILE, "refseq.fa");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_REFSEQ_INPUT_TYPE, "fasta");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_REFSEQ_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_FILE, "reads.sam");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_READS_INPUT_TYPE, "sam");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_READS_SKIP_VERT, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_OUTPUT_FILE, "-");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_ACTION, "assembly");

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(GASSM_OPTION_VERBOSE, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(GASSM_OPTION_REGION_START, 534000); // 534000

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(GASSM_OPTION_REGION_LENGTH, 2000); // 2000

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GASSM_OPTION_WEIGHT_THRESHOLD, 3);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(GASSM_OPTION_NUM_PATHS, 4);

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
		ret = option_set_description_const(GASSM_OPTION_REGION_START, GASSM_OPTION_REGION_START_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_REGION_LENGTH, GASSM_OPTION_REGION_LENGTH_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_WEIGHT_THRESHOLD, GASSM_OPTION_WEIGHT_THRESHOLD_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(GASSM_OPTION_NUM_PATHS, GASSM_OPTION_NUM_PATHS_DESC);
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
		option_set_shortcut(GASSM_OPTION_REGION_START, 'r');
		option_set_shortcut(GASSM_OPTION_REGION_LENGTH, 'l');
		option_set_shortcut(GASSM_OPTION_WEIGHT_THRESHOLD, 'w');
		option_set_shortcut(GASSM_OPTION_NUM_PATHS, 'p');
	}

	return ret;
}

typedef enum _ETimeMeasuredActionType{
	tmatReadRefSeq,
	tmatReadReads,
	tmatRefseqToRegions,
	tmatGraphCreate,
	tmatGraphProcessRefSeq,
	tmatProcessReads,
	tmatThresholds,
	tmatProbabilities,
	tmatPathSearch,
	tmatSSWSimple,
	tmatSSWClever,
	tmatMax,
} EMeasuredActionType, *PEMeasuredActionType;

static int _timeArray[tmatMax];
static char *_timeArrayStr[] = {
	"Reading ref. seq.",
	"Reading reads",
	"Ref. seq. to regions",
	"Graph create",
	"Refseq to graph",
	"Reads to graph",
	"Threshold optimization",
	"Edge probabilities",
	"Path search",
	"SSW simple",
	"SSW clever",
};


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(_timeArray, 0, sizeof(_timeArray));
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
					uint64_t regionStart = 0;
					uint64_t regionLength = 0;
					char *action = NULL;
					boolean refSeqSkipVertices = FALSE;
					boolean readsSkipVertices = FALSE;
					char *refSeqType = NULL;
					char *refSeqFileName = NULL;
					char *readsType = NULL;
					char *readsFileName = NULL;
					uint32_t wightThreshold = 0;
					uint32_t numPaths = 0;
					boolean verbose = FALSE;

					ret = option_get_Boolean(GASSM_OPTION_VERBOSE, &verbose);
					if (ret == ERR_SUCCESS)
						ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
					
					if (ret == ERR_SUCCESS)
						ret = option_get_UInt64(GASSM_OPTION_REGION_START, &regionStart);

					if (ret == ERR_SUCCESS)
						ret = option_get_UInt64(GASSM_OPTION_REGION_LENGTH, &regionLength);
					
					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_ACTION, &action);

					if (ret == ERR_SUCCESS)
						ret = option_get_Boolean(GASSM_OPTION_REFSEQ_SKIP_VERT, &refSeqSkipVertices);

					if (ret == ERR_SUCCESS)
						ret = option_get_Boolean(GASSM_OPTION_READS_SKIP_VERT, &readsSkipVertices);

					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_REFSEQ_INPUT_FILE, &refSeqFileName);

					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_REFSEQ_INPUT_TYPE, &refSeqType);

					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_READS_INPUT_FILE, &readsFileName);

					if (ret == ERR_SUCCESS)
						ret = option_get_String(GASSM_OPTION_READS_INPUT_TYPE, &readsType);

					if (ret == ERR_SUCCESS)
						ret = option_get_UInt32(GASSM_OPTION_WEIGHT_THRESHOLD, &wightThreshold);

					if (ret == ERR_SUCCESS)
						ret = option_get_UInt32(GASSM_OPTION_NUM_PATHS, &numPaths);

					if (verbose)
						options_print();

					if (ret == ERR_SUCCESS) {
						char *refSeq = NULL;
						size_t refSeqLen = 0;
						PONE_READ *reads = NULL;
						size_t readCount = 0;
						int t = GetTickCount();

						ret = input_get_refseq(refSeqFileName, refSeqType, &refSeq, &refSeqLen);
						if (ret == ERR_SUCCESS) {
							_timeArray[tmatReadRefSeq] = GetTickCount() - t;
							t = GetTickCount();
							if (regionLength == 0)
								regionLength = refSeqLen - regionStart;

							ret = input_get_reads(readsFileName, readsType, regionStart, regionLength, &reads, &readCount);
							if (ret == ERR_SUCCESS) {
								PACTIVE_REGION regions = NULL;
								size_t regionCount = 0;

								_timeArray[tmatReadReads] = GetTickCount() - t;
								t = GetTickCount();
								ret = input_refseq_to_regions(refSeq, refSeqLen, &regions, &regionCount);
								if (ret == ERR_SUCCESS) {
									_timeArray[tmatRefseqToRegions] = GetTickCount() - t;
									t = GetTickCount();
									if (strcasecmp(action, "assembly") == 0) {
										PKMER_GRAPH g = NULL;
										
										ret = kmer_graph_create(kmerSize, &g);
										if (ret == ERR_SUCCESS) {
											_timeArray[tmatGraphCreate] = GetTickCount() - t;
											t = GetTickCount();
											ret = kmer_graph_parse_ref_sequence(g, refSeq + regionStart, regionLength, refSeqSkipVertices);
											if (ret == ERR_SUCCESS) {
												_timeArray[tmatGraphProcessRefSeq] = GetTickCount() - t;
												t = GetTickCount();
												ret = kmer_graph_parse_reads(g, reads, readCount, readsSkipVertices);
												if (ret == ERR_SUCCESS) {
													size_t pathCount = 0;
													PKMER_GRAPH_PATH paths = NULL;
													
													_timeArray[tmatProcessReads] = GetTickCount() - t;
													t = GetTickCount();
													if (wightThreshold > 0) {
														kmer_graph_delete_edges_under_threshold(g, wightThreshold);
														_timeArray[tmatThresholds] = GetTickCount() - t;
														t = GetTickCount();
													}

													kmer_graph_compute_edge_probablities(g);
													_timeArray[tmatProbabilities] = GetTickCount() - t;
													t = GetTickCount();
													ret = kmer_graph_find_best_paths(g, numPaths, regionLength - kmerSize + 2, &paths, &pathCount);
													if (ret == ERR_SUCCESS) {
														_timeArray[tmatPathSearch] = GetTickCount() - t;
														t = GetTickCount();
														for (size_t i = 0; i < pathCount; ++i)
															printf("%lE: %u %s\n", paths[i].Weight, paths[i].Length, paths[i].Sequence);

														int64_t totalSSWSimpleTime = 0;
														int64_t totalSSWCleverTime = 0;
														int i;
// #pragma omp parallel for private(i), shared(numPaths, refSeq, regionStart, regionLength, paths, pathCount)
														for (i = 0; i < pathCount; ++i) {
															char *opString = NULL;
															size_t opStringLen = 0;
															
															t = GetTickCount();
															ret = ssw_clever(paths[i].Sequence, paths[i].Length, refSeq + regionStart, regionLength, 2, -1, -1, &opString, &opStringLen);
															if (ret == ERR_SUCCESS) {
																char *opString2 = NULL;
																size_t opStringLen2 = 0;

																totalSSWCleverTime += (GetTickCount() - t);
																t = GetTickCount();
																/*
																ret = ssw_simple(paths[i].Sequence, paths[i].Length, refSeq + regionStart, regionLength, 2, -1, -1, &opString2, &opStringLen2);
																if (ret == ERR_SUCCESS) {
																	totalSSWSimpleTime += (GetTickCount() - t);
																	if (opStringLen != opStringLen2 || strcasecmp(opString, opString2) != 0)
																		fprintf(stderr, "SSW-MISMATCH:\n%s\n%s\n", opString, opString2);						

																	utils_free(opString2);
																}
																*/
																write_differences(refSeq, regionStart, regionLength, opString, paths[i].Sequence);
																utils_free(opString);
																printf("\n");
															}

															if (ret != ERR_SUCCESS)
																break;
														}

														if (ret == ERR_SUCCESS) {
															_timeArray[tmatSSWSimple] = totalSSWSimpleTime / pathCount;
															_timeArray[tmatSSWClever] = totalSSWCleverTime / pathCount;
														}

														kmer_graph_paths_free(paths, pathCount);
													}

													if (ret == ERR_SUCCESS) {
														for (size_t i = 0; i < tmatMax; ++i) {
															fprintf(stderr, "%s: %u ms\n", _timeArrayStr[i], _timeArray[i]);
														}

														fprintf(stderr, "\n");
														kmer_graph_delete_1to1_vertices(g);
														kmer_graph_print(stdout, g);
													}
												}
											}

											kmer_graph_destroy(g);
										}
									} else if (strcasecmp(action, "regions") == 0) {
										PACTIVE_REGION ar = regions;

										fprintf(stdout, "%u regions\n", regionCount);
										for (size_t i = 0; i < regionCount; ++i) {
											fprintf(stdout, "  Type: %s, Offset: %llu, length: %llu\n", (ar->Type == artUnknown ? "Unknown" : "Valid" ), ar->Offset, ar->Length);
											++ar;
										}
									} else if (strcasecmp(action, "sswtest") == 0) {
										char *opString = NULL;
										size_t opStringLen = 0;

										ret = ssw_simple("TGAAAAAAA", 9, "AAAAAATA", 8, 2, -1, -1, &opString, &opStringLen);
										if (ret == ERR_SUCCESS) {
											printf("%s", opString);
											utils_free(opString);
										}

										getchar();
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
