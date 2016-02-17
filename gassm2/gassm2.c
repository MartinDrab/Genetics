
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

	if (ret == ERR_SUCCESS) {
		ret = option_get_String(PROGRAM_OPTION_SEQUENCE, &Options->ReferenceSequence);
		if (ret == ERR_SUCCESS && *Options->ReferenceSequence == '\0')
			ret = option_get_String(PROGRAM_OPTION_SEQFILE, &Options->RefSeqFile);
	}

	if (ret == ERR_SUCCESS) {
		char *readFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_READFILE, &readFile);
		if (ret == ERR_SUCCESS && *readFile != '\0') {
			ret = input_get_reads(readFile, "sam", Options->RegionStart, Options->RegionLength, &Options->Reads, &Options->ReadCount);
		}
	}

	return ret;
}


static void _compute_graph(const struct _PROGRAM_OPTIONS *Options, const char *RefSeq, const char *Alternate, const size_t AlternateLen, PPROGRAM_STATISTICS Statistics)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

		ret = kmer_graph_create(Options->KMerSize, &g);
		if (ret == ERR_SUCCESS) {
			ret = kmer_graph_parse_ref_sequence(g, RefSeq, Options->RegionLength, FALSE, Options->Threshold);
			if (ret == ERR_SUCCESS) {
				Statistics->VertexCount = kmer_graph_get_vertex_count(g);
				Statistics->EdgeCount = kmer_graph_get_edge_count(g);
				Statistics->CycleCount = kmer_graph_get_cycle_count(g);
				ret = kmer_graph_parse_reads(g, Options->Reads, Options->ReadCount, FALSE);
				if (ret == ERR_SUCCESS) {
					size_t l = 0;
					char *s = NULL;

					kmer_graph_delete_edges_under_threshold(g, Options->Threshold);
					kmer_graph_delete_trailing_things(g);
						ret = kmer_graph_delete_1to1_vertices(g);
						if (ret == ERR_SUCCESS) {
							ret = kmer_graph_get_seq(g, &s, &l);
							if (ret == ERR_SUCCESS) {
								if ((l == AlternateLen) && memcmp(s, Alternate, l)) {
									fprintf(stderr, "REFSEQ   = %s\n", Alternate);
									fprintf(stderr, "GRAPHSEQ = %s\n", s);
									kmer_graph_print(stderr, g);
								} else {
									if (Options->PrintResults) {
//										printf("SEQ      = %s\n", Alternate);
//										printf("GRAPHSEQ = %s\n", s);
										kmer_graph_print(stderr, g);
									}
								}

								utils_free(s);
							} else printf("kmer_graph_get_seq(): %u\n", ret);
						} else printf("kmer_graph_delete_1to1_vertices(): %u\n", ret);
				} else printf("kmer_graph_parse_reads(): %u\n", ret);
			} else printf("kmer_graph_parse_ref_sequence(): %u\n", ret);

			kmer_graph_destroy(g);
		} else printf("kmer_graph_create(): %u\n", ret);

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


static ERR_VALUE _create_alternatce_sequence(const char *RefSeq, const size_t RefSeqLen, char **Alternate, size_t *AlternateLen)
{
	size_t tmpAlternateLen = RefSeqLen;
	char *rsCopy = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc(tmpAlternateLen + 1, sizeof(char), &rsCopy);
	if (ret == ERR_SUCCESS) {
		memcpy(rsCopy, RefSeq, tmpAlternateLen*sizeof(char));
		rsCopy[tmpAlternateLen] = '\0';
		*Alternate = rsCopy;
		*AlternateLen = tmpAlternateLen;
	}

	return ret;
}


static ERR_VALUE _test_with_reads(PPROGRAM_OPTIONS Options, const char *RefSeq, PPROGRAM_STATISTICS Statistics)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	
	ret = ERR_SUCCESS;
	memset(Statistics, 0, sizeof(PROGRAM_STATISTICS));
	if (Options->ReadCount > 0 && Options->TestReadCycles > 0) {
		for (size_t j = 0; j < Options->TestReadCycles; ++j) {
			char *alternate = NULL;
			size_t alternateLen = 0;

			ret = _create_alternatce_sequence(RefSeq, Options->RegionLength, &alternate, &alternateLen);
			if (ret == ERR_SUCCESS) {
				ret = read_set_generate_from_sequence(alternate, alternateLen, Options->ReadLength, Options->ReadCount, &Options->Reads);
				if (ret == ERR_SUCCESS) {
					PROGRAM_STATISTICS stats;

					_compute_graph(Options, RefSeq, alternate, alternateLen, &stats);
					Statistics->VertexCount += stats.VertexCount;
					Statistics->EdgeCount += stats.EdgeCount;
					Statistics->CycleCount += stats.CycleCount;
					read_set_destroy(Options->Reads, Options->ReadCount);
				}

				utils_free(alternate);
			}

			if (ret != ERR_SUCCESS)
				break;
		}

		if (ret == ERR_SUCCESS) {
			Statistics->VertexCount /= Options->TestReadCycles;
			Statistics->EdgeCount /= Options->TestReadCycles;
			Statistics->CycleCount /= Options->TestReadCycles;
		}
	} else _compute_graph(Options, RefSeq, RefSeq, Options->RegionLength, Statistics);
	
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
									
									if ((i + 1) % (po.TestCount / 100) == 0)
										printf(".");

									memset(rs, 0, po.RegionLength*sizeof(char));
									for (uint32_t j = 0; j < po.RegionLength; ++j)
										rs[j] = _rand_nucleotide();

									po.ReferenceSequence = rs;
									ret = _test_with_reads(&po, rs, &tmpStats);
									if (ret == ERR_SUCCESS) {
										st.CycleCount += tmpStats.CycleCount;
										st.EdgeCount += tmpStats.EdgeCount;
										st.VertexCount += tmpStats.VertexCount;
										st.VertexVariance += (tmpStats.VertexCount*tmpStats.VertexCount);
										st.EdgeVariance += (tmpStats.EdgeCount*tmpStats.EdgeCount);
										st.CycleVariance += (tmpStats.CycleCount*tmpStats.CycleCount);
									} else printf("_test_with_reads(): %u", ret);
								}

								if (ret == ERR_SUCCESS) {
									st.CycleCount /= po.TestCount;
									st.EdgeCount /= po.TestCount;
									st.VertexCount /= po.TestCount;
									st.VertexVariance = (st.VertexVariance / po.TestCount) - st.VertexCount*st.VertexCount;
									st.EdgeVariance = (st.EdgeVariance / po.TestCount) - st.EdgeCount*st.EdgeCount;
									st.CycleVariance = (st.CycleVariance / po.TestCount) - st.CycleCount*st.CycleCount;
									printf("AVG: Vertices: (%" PRIu64 " %lf), Edges: (%" PRIu64 " %lf), Cycles: (%" PRIu64 " %lf)\n", st.VertexCount, sqrt(st.VertexVariance), st.EdgeCount, sqrt(st.EdgeVariance), st.CycleCount, sqrt(st.CycleVariance));
								}

								utils_free(rs);
							}
						} else {
							size_t refSeqLen = 0;
							boolean explicitSequence = (*po.ReferenceSequence != '\0');
							FASTA_FILE seqFile;

							if (!explicitSequence) {
								ret = fasta_load(po.RefSeqFile, &seqFile);
								if (ret == ERR_SUCCESS) {
									ret = fasta_read_seq(&seqFile, &po.ReferenceSequence, &refSeqLen);
									if (ret != ERR_SUCCESS)
										fasta_free(&seqFile);
								}
							} else refSeqLen = strlen(po.ReferenceSequence);

							if (ret == ERR_SUCCESS) {
								do {
									size_t regionCount = 0;
									PACTIVE_REGION regions = NULL;
									uint64_t numberOfAttempts = 0;
									char *origRefSeq = po.ReferenceSequence;

									ret = input_refseq_to_regions(po.ReferenceSequence, refSeqLen, &regions, &regionCount);
									if (ret == ERR_SUCCESS) {
										printf("Going through a reference sequence of length %" PRIu64 " with %u regions...\n", (uint64_t)refSeqLen, regionCount);
										printf("%u test read cycles with %u reads of length %u...\n", po.TestReadCycles, po.ReadCount, po.ReadLength);
										for (size_t i = 0; i < regionCount; ++i) {
											PACTIVE_REGION pa = regions + i;

											printf("Region #%u: Offset: %" PRIu64 ", Length %" PRIu64 "\n", i, pa->Offset, pa->Length);
											if (pa->Type == artValid && pa->Length >= po.RegionLength) {
												int j = 0;
												
												po.ReferenceSequence = pa->Sequence;
#pragma omp parallel for shared(po, st, numberOfAttempts)	
//												for (uint64_t j = 0; j < pa->Length - po.RegionLength; j += po.TestStep) {
												for (j = 0; j < (int)(pa->Length - po.RegionLength); j += (int)po.TestStep) {
												const char *refSeq = pa->Sequence + j;
													PROGRAM_STATISTICS tmpstats;

													ret = _test_with_reads(&po, refSeq, &tmpstats);
													if (ret == ERR_SUCCESS) {
														++numberOfAttempts;
														st.VertexCount += tmpstats.VertexCount;
														st.EdgeCount += tmpstats.EdgeCount;
														st.CycleCount += tmpstats.CycleCount;
														st.VertexVariance += (tmpstats.VertexCount*tmpstats.VertexCount);
														st.EdgeVariance += (tmpstats.EdgeCount*tmpstats.EdgeCount);
														st.CycleVariance += (tmpstats.CycleCount*tmpstats.CycleCount);
													}
												}
											}

											++pa;
										}

										if (ret == ERR_SUCCESS) {
											st.CycleCount /= numberOfAttempts;
											st.EdgeCount /= numberOfAttempts;
											st.VertexCount /= numberOfAttempts;
											st.VertexVariance = (st.VertexVariance / numberOfAttempts) - st.VertexCount*st.VertexCount;
											st.EdgeVariance = (st.EdgeVariance / numberOfAttempts) - st.EdgeCount*st.EdgeCount;
											st.CycleVariance = (st.CycleVariance / numberOfAttempts) - st.CycleCount*st.CycleCount;
											printf("AVG: Vertices: (%" PRIu64 " %lf), Edges: (%" PRIu64 " %lf), Cycles: (%" PRIu64 " %lf)\n", st.VertexCount, sqrt(st.VertexVariance), st.EdgeCount, sqrt(st.EdgeVariance), st.CycleCount, sqrt(st.CycleVariance));
										}

										input_free_regions(regions, regionCount);
									}

									po.ReferenceSequence = origRefSeq;
									if (!explicitSequence) {
										utils_free(po.ReferenceSequence);
										ret = fasta_read_seq(&seqFile, &po.ReferenceSequence, &refSeqLen);
									}
								} while (ret == ERR_SUCCESS && !explicitSequence);

								if (ret == ERR_NO_MORE_ENTRIES)
									ret = ERR_SUCCESS;

								if (!explicitSequence)
									fasta_free(&seqFile);
							}
						}
					} else if (*po.RefSeqFile != '\0') {
						FASTA_FILE seqFile;

						ret = fasta_load(po.RefSeqFile, &seqFile);
						if (ret == ERR_SUCCESS) {
							size_t refSeqLen = 0;

							ret = fasta_read_seq(&seqFile, &po.ReferenceSequence, &refSeqLen);
							if (ret == ERR_SUCCESS) {
								size_t regionCount = 0;
								PACTIVE_REGION regions = NULL;

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

											printf("kmer size: %u\n", po.KMerSize);
											printf("Active region (%" PRIu64 "; %u; %u)...\n", po.RegionStart, po.RegionLength, index);
											_compute_graph(&po, po.ReferenceSequence, po.ReferenceSequence, po.RegionLength, &st);
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
