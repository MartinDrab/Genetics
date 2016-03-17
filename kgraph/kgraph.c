
#include <stdio.h>
#include <omp.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "kgraph.h"



static ERR_VALUE _options_set_default(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_String(KGRAPH_OPTION_INPUT_FILE, "input.fa");
	if (ret == ERR_SUCCESS)
		ret = option_add_String(KGRAPH_OPTION_INPUT_TYPE, "fasta");

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(KGRAPH_OPTION_REGION_START, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(KGRAPH_OPTION_REGION_LENGTH, 0);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(KGRAPH_OPTION_UNIT_SIZE, 2000);

	if (ret == ERR_SUCCESS)
		ret = option_add_Float(KGRAPH_OPTION_UNIT_OVERLAP, 0.25);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(KGRAPH_OPTION_MIN_K, 5);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(KGRAPH_OPTION_MAX_K, 50);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(KGRAPH_OPTION_STEP_K, 1);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(KGRAPH_OPTION_SKIP_VERTICES, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(KGRAPH_OPTION_HELP, FALSE);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(KGRAPH_OPTION_VERBOSE, FALSE);

	if (ret == ERR_SUCCESS) {
		ret = option_set_description_const(KGRAPH_OPTION_INPUT_FILE, KGRAPH_OPTION_INPUT_FILE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_INPUT_TYPE, KGRAPH_OPTION_INPUT_TYPE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_REGION_START, KGRAPH_OPTION_REGION_START_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_REGION_LENGTH, KGRAPH_OPTION_REGION_LENGTH_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_UNIT_SIZE, KGRAPH_OPTION_UNIT_SIZE_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_UNIT_OVERLAP, KGRAPH_OPTION_UNIT_OVERLAP_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_MIN_K, KGRAPH_OPTION_MIN_K_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_MAX_K, KGRAPH_OPTION_MAX_K_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_STEP_K, KGRAPH_OPTION_STEP_K_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_SKIP_VERTICES, KGRAPH_OPTION_SKIP_VERTICES_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_HELP, KGRAPH_OPTION_HELP_DESC);
		assert(ret == ERR_SUCCESS);
		ret = option_set_description_const(KGRAPH_OPTION_VERBOSE, KGRAPH_OPTION_VERBOSE_DESC);
		assert(ret == ERR_SUCCESS);
	}


	return ret;
}


ERR_VALUE _options_to_record(PKGRAPH_OPTIONS_RECORD Record)
{
	ERR_VALUE ret = ERR_SUCCESS;

	ret = option_get_String(KGRAPH_OPTION_INPUT_FILE, &Record->InputFile);
	if (ret == ERR_SUCCESS)
		ret = option_get_String(KGRAPH_OPTION_INPUT_TYPE, &Record->Inputtype);
	
	if (ret == ERR_SUCCESS)
		ret = option_get_UInt64(KGRAPH_OPTION_REGION_START, &Record->Regionstart);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt64(KGRAPH_OPTION_REGION_LENGTH, &Record->regionLength);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(KGRAPH_OPTION_MIN_K, &Record->MinK);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(KGRAPH_OPTION_MAX_K, &Record->MaxK);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(KGRAPH_OPTION_STEP_K, &Record->StepK);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(KGRAPH_OPTION_UNIT_SIZE, &Record->UnitSize);

	if (ret == ERR_SUCCESS)
		ret = option_get_Float(KGRAPH_OPTION_UNIT_OVERLAP_DESC, &Record->UnitOverlap);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(KGRAPH_OPTION_SKIP_VERTICES, &Record->SkipVertices);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(KGRAPH_OPTION_VERBOSE, &Record->Verbose);

	if (ret == ERR_SUCCESS)
		ret = option_get_Boolean(KGRAPH_OPTION_HELP, &Record->Help);

	return ret;
}


ERR_VALUE _compute_backwards_edge_count(const char *Region, const size_t RegionLength, const uint32_t MinK, const uint32_t MaxK, const uint32_t StepK, const boolean SkipVertices, uint32_t *ResultArray)
{
	int i = MinK;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

#pragma omp parallel for shared(ret, Region, RegionLength, MinK, MaxK, StepK, SkipVertices, ResultArray)
	for (i = MinK; i <= MaxK; i += StepK) {
		PKMER_GRAPH g = NULL;
		const uint32_t index = (i - MinK) / StepK;

		ret = kmer_graph_create(i, &g);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_parse_ref_sequence(g, Region, RegionLength, 0);
			if (ret == ERR_SUCCESS)
				ResultArray[index] = g->NumberOfBackwardEdges;

			kmer_graph_destroy(g);
		}
	}

	return ret;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = options_module_init(29);
	if (ret == ERR_SUCCESS) {
		ret = _options_set_default();
		if (ret == ERR_SUCCESS) {
			options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS) {
				KGRAPH_OPTIONS_RECORD or;

				ret = _options_to_record(&or);
				if (ret == ERR_SUCCESS) {
					char *seq = NULL;
					size_t len = 0;

					if (or.Verbose)
						options_print();

					if (!or.Help) {
						ret = input_get_refseq(or.InputFile, or.Inputtype, &seq, &len);
						if (ret == ERR_SUCCESS) {
							size_t regionCount = 0;
							PACTIVE_REGION regions = NULL;

							ret = input_refseq_to_regions(seq, len, &regions, &regionCount);
							if (ret == ERR_SUCCESS) {
								size_t regionIndex = 0;
								uint64_t regionOffset = 0;

								if (or.regionLength == 0)
									or.regionLength = len - or.Regionstart;

								ret = input_get_region_by_offset(regions, regionCount, or.Regionstart, &regionIndex, &regionOffset);
								if (ret == ERR_SUCCESS) {
									PACTIVE_REGION current = regions + regionIndex;
									const uint64_t stepSize = (uint64_t)(or.UnitSize * (1 - or.UnitOverlap));
									uint32_t *resultArray = NULL;
									const size_t kCount = (or.MaxK - or.MinK + or.StepK) / or.StepK;

									ret = utils_calloc(kCount, sizeof(uint32_t), &resultArray);
									if (ret == ERR_SUCCESS) {
										printf("[");
										while (or.regionLength > 0) {
											uint64_t remainingLength = min(current->Length - regionOffset, or.regionLength);

											or.regionLength -= remainingLength;
											if (current->Type != artUnknown) {
												const char *pos = current->Sequence + regionOffset;

												do {
													if (remainingLength >= stepSize) {
														ret = _compute_backwards_edge_count(pos, stepSize, or.MinK, or.MaxK, or.StepK, or.SkipVertices, resultArray);
														printf("%u ", (uint32_t)(pos - seq));
														for (size_t i = 0; i < kCount; ++i)
															printf("%u ", resultArray[i]);
														
														pos += stepSize;
														remainingLength -= stepSize;
													}
													else {
														ret = _compute_backwards_edge_count(pos, remainingLength, or.MinK, or.MaxK, or.StepK, or.SkipVertices, resultArray);
														printf("%u ", (uint32_t)(pos - seq));
														for (size_t i = 0; i < kCount; ++i)
															printf("%u ", resultArray[i]);

														remainingLength = 0;
													}

													printf("\n");
												} while (remainingLength > 0);
											}

											++current;
											regionOffset = 0;
										}

										printf("]\n");
										utils_free(resultArray);
									}
								}

								input_free_regions(regions, regionCount);
							}

							utils_free(seq);
						}
					} else {
						printf("Usage: kgraph [options]\n");
						options_print_help();
					}
				}
			}
		}

		options_module_finit();
	}


	return ret;
}
