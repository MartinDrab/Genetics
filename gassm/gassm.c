
#include <stdio.h>
#include <malloc.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "kmer.h"
#include "kmer-table.h"
#include "kmer-edge.h"
#include "kmer-graph.h"
#include "gassm.h"



static ERR_VALUE _set_default_values(void)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(GASSM_OPTION_KMER_SIZE, 3);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(GASSM_OPTION_TEST_INPUT_FILE, "input.txt");

	return ret;
}

static ERR_VALUE _parse_test_data(char *Data, size_t DataLength, char **RefSeq, char ***Reads, size_t *ReadsCount)
{
	char *tmpRefSeq = NULL;
	size_t refSeqLen = 0;

	char *tmp = Data;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	// 1) Capture the reference sequence
	while (DataLength && (*tmp != '\r' && *tmp != '\n' && *tmp != 26)) {
		++tmp;
		--DataLength;
	}

	refSeqLen = tmp - Data;
	if (DataLength > 0 && refSeqLen > 0) {
		ret = utils_preallocate_string(refSeqLen, &tmpRefSeq);
		if (ret == ERR_SUCCESS) {
			memcpy(tmpRefSeq, Data, refSeqLen*sizeof(char));
			// 2) Determine number of reads
			size_t tmpReadsCount = 0;
			char **tmpReadArray = NULL;

			while (DataLength > 0 && (*tmp == '\n' || *tmp == '\r')) {
				++tmp;
				--DataLength;
			}

			if (*tmp != 26 && DataLength > 0) {
				char *readStart = tmp;

				while (DataLength > 0) {
					while (DataLength > 0 && (*tmp != '\r' && *tmp != '\n' && *tmp != 26)) {
						++tmp;
						--DataLength;
					}

					++tmpReadsCount;
					while (DataLength > 0 && (*tmp == '\n' || *tmp == '\r')) {
						++tmp;
						--DataLength;
					}
				}

				// 3) Copy the reads
				tmpReadArray = (char **)calloc(tmpReadsCount, sizeof(char *));
				if (tmpReadArray != NULL) {
					size_t index = 0;
					
					while (ret == ERR_SUCCESS && index < tmpReadsCount) {
						tmp = readStart;
						while ((*tmp != '\r' && *tmp != '\n' && *tmp != 26))
							++tmp;

						char *read = NULL;
						size_t readLen = tmp - readStart;
						ret = utils_preallocate_string(readLen, &read);
						if (ret == ERR_SUCCESS) {
							memcpy(read, readStart, readLen*sizeof(char));
							while ((*tmp == '\n' || *tmp == '\r'))
								++tmp;

							tmpReadArray[index] = read;
							++index;
							readStart = tmp;
						}
					}

					if (ret != ERR_SUCCESS) {
						for (size_t i = 0; i < index; ++i)
							free(tmpReadArray[i]);
					
						free(tmpReadArray);
					}
				} else ret = ERR_OUT_OF_MEMORY;
			} else ret = ERR_SUCCESS;

			if (ret == ERR_SUCCESS) {
				*ReadsCount = tmpReadsCount;
				*Reads = tmpReadArray;
				*RefSeq = tmpRefSeq;
			}

			if (ret != ERR_SUCCESS)
				utils_free_string(tmpRefSeq);
		}
	}

	return ret;
}


void _free_test_data(char *RefSeq, char **Reads, const size_t ReadsCount)
{
	for (size_t i = 0; i < ReadsCount; ++i)
		free(Reads[i]);

	free(Reads);
	free(RefSeq);

	return;
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
				uint32_t kmerSize = 0;
				char *testFile = NULL;

				ret = option_get_UInt32(GASSM_OPTION_KMER_SIZE, &kmerSize);
				if (ret == ERR_SUCCESS)
					ret = option_get_String(GASSM_OPTION_TEST_INPUT_FILE, &testFile);
					
				if (ret == ERR_SUCCESS) {
					char *testData = NULL;
					size_t testDataSize = 0;

					options_print();
					ret = utils_file_read(testFile, &testData, &testDataSize);
					if (ret == ERR_SUCCESS) {
						char *refSeq = NULL;
						char **reads = NULL;
						size_t readCount = 0;

						ret = _parse_test_data(testData, testDataSize, &refSeq, &reads, &readCount);
						if (ret == ERR_SUCCESS) {
							printf("Ref. sequence: %s\n", refSeq);
							printf("Number of reads: %zu\n", readCount);
							for (size_t i = 0; i < readCount; ++i)
								printf("Read #%zu: %s\n", i, reads[i]);
						
							PKMER_GRAPH g = NULL;
						
							ret = kmer_graph_create(kmerSize, &g);
							if (ret == ERR_SUCCESS) {
								PKMER sourceKMer = NULL;
								PKMER destKMer = NULL;

								printf("\nProcessing the reference sequence...");
								sourceKMer = kmer_alloc(kmerSize, refSeq);
								if (sourceKMer != NULL) {
									kmer_back(sourceKMer, '^');
									printf("\nAdding vertext: ");
									kmer_print(sourceKMer);
									ret = kmer_graph_add_vertex(g, sourceKMer);
									if (ret == ERR_SUCCESS) {
										destKMer = kmer_alloc(kmerSize, refSeq);
										if (destKMer != NULL) {
											size_t refSeqLen = strlen(refSeq);
											
											for (size_t i = 0; i <= refSeqLen - kmerSize + 1; ++i) {
												if (i > refSeqLen - kmerSize)
													kmer_advance(destKMer, '$');
												else kmer_init(destKMer, refSeq + i);
												
												printf("\nAdding vertext: ");
												kmer_print(destKMer);
												ret = kmer_graph_add_vertex(g, destKMer);
												if (ret == ERR_ALREADY_EXISTS) {
													ret = ERR_SUCCESS;
													printf("\nAlready exists");
												}

												if (ret == ERR_SUCCESS) {
													printf("\nAdding edge: ");
													kmer_print(sourceKMer);
													printf(" ---> ");
													kmer_print(destKMer);
													ret = kmer_graph_add_edge(g, sourceKMer, destKMer, 0);
													if (ret == ERR_ALREADY_EXISTS) {
														ret = ERR_SUCCESS;
														printf("\nAlready exists");
													}
												
													kmer_init_from_kmer(sourceKMer, destKMer);
												}

												if (ret != ERR_SUCCESS)
													break;
											}

											if (ret == ERR_SUCCESS) {
												for (size_t j = 0; j < readCount; ++j) {
													char *currentRead = reads[j];
													size_t currentReadlen = strlen(currentRead);

													printf("\n\nProcessing read #%zu: %s", j, currentRead);
													kmer_init(sourceKMer, currentRead);
													printf("\nAdding vertext: ");
													kmer_print(sourceKMer);
													ret = kmer_graph_add_vertex(g, sourceKMer);
													if (ret == ERR_ALREADY_EXISTS) {
														ret = ERR_SUCCESS;
														printf("\nAlready exists");
													}

													for (size_t i = 1; i < currentReadlen - kmerSize + 1; ++i) {
														kmer_init(destKMer, currentRead + i);
														printf("\nAdding vertext: ");
														kmer_print(destKMer);
														ret = kmer_graph_add_vertex(g, destKMer);
														if (ret == ERR_ALREADY_EXISTS) {
															ret = ERR_SUCCESS;
															printf("\nAlready exists");
														}

														if (ret == ERR_SUCCESS) {
															printf("\nAdding edge: ");
															kmer_print(sourceKMer);
															printf(" ---> ");
															kmer_print(destKMer);
															ret = kmer_graph_add_edge(g, sourceKMer, destKMer, 1);
															if (ret == ERR_ALREADY_EXISTS) {
																PKMER_EDGE edge = NULL;

																edge = kmer_graph_get_edge(g, sourceKMer, destKMer);
																edge->Weight++;
																ret = ERR_SUCCESS;
																printf("\nAlready exists");
															}

															kmer_init_from_kmer(sourceKMer, destKMer);
														}

														if (ret != ERR_SUCCESS)
															break;
													}
												}

												if (ret == ERR_SUCCESS) {
													printf("\nEdges (and vertices):");
													kmer_edge_table_print(g->EdgeTable);
												}
											}
										} else ret = ERR_OUT_OF_MEMORY;
									}

									kmer_free(sourceKMer);
								} else ret = ERR_OUT_OF_MEMORY;

								kmer_graph_destroy(g);
							}

							_free_test_data(refSeq, reads, readCount);
						}

						free(testData);
					}

				}
			}
		}

		options_module_finit();
	}

	getchar();

	return ret;
}
