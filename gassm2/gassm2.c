
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





UTILS_TYPED_CALLOC_FUNCTION(GEN_ARRAY_ONE_READ)
UTILS_TYPED_CALLOC_FUNCTION(GEN_ARRAY_VARIANT_CALL)


static PUTILS_LOOKASIDE *_vertexLAs;
static PUTILS_LOOKASIDE *_edgeLAs;



static ERR_VALUE _init_default_values()
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = option_add_UInt32(PROGRAM_OPTION_KMERSIZE, 31);
	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_SEQFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt64(PROGRAM_OPTION_SEQSTART, (uint64_t)-1);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_SEQLEN, 2000);


	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_TEST_STEP, 1500);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_THRESHOLD, 2);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_READFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_OUTPUT_DIRECTORY, ".");

	if (ret == ERR_SUCCESS)
		ret = option_add_String(PROGRAM_OPTION_VCFFILE, "\0");

	if (ret == ERR_SUCCESS)
		ret = option_add_Int32(PROGRAM_OPTION_OMP_THREADS, omp_get_num_procs());

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt8(PROGRAM_OPTION_READ_POS_QUALITY, 10);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_CONNECT_REFSEQ, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_CONNECT_READS, FALSE);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_BUBBLE_MERGING, FALSE);
	
	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_LINEAR_SHRINK, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_Boolean(PROGRAM_OPTION_NO_HELPER_VERTICES, FALSE);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_MISSING_EDGE_PENALTY, 3);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_BACKWARD_REFSEQ_PENALTY, 2);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_MAX_PATHS, 1);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_READ_MAX_ERROR_RATE, 20);

	if (ret == ERR_SUCCESS)
		ret = option_add_UInt32(PROGRAM_OPTION_READ_STRIP, 5);

	option_set_description_const(PROGRAM_OPTION_KMERSIZE, PROGRAM_OPTION_KMERSIZE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQFILE, PROGRAM_OPTION_SEQFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQSTART, PROGRAM_OPTION_SEQSTART_DESC);
	option_set_description_const(PROGRAM_OPTION_SEQLEN, PROGRAM_OPTION_SEQLEN_DESC);
	option_set_description_const(PROGRAM_OPTION_TEST_STEP, PROGRAM_OPTION_TEST_STEP_DESC);
	option_set_description_const(PROGRAM_OPTION_THRESHOLD, PROGRAM_OPTION_THRESHOLD_DESC);
	option_set_description_const(PROGRAM_OPTION_READFILE, PROGRAM_OPTION_READFILE_DESC);
	option_set_description_const(PROGRAM_OPTION_OUTPUT_DIRECTORY, PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC);
	option_set_description_const(PROGRAM_OPTION_VCFFILE, PROGRAM_OPTION_VCFFILE_DESC);

	option_set_shortcut(PROGRAM_OPTION_KMERSIZE, 'k');
	option_set_shortcut(PROGRAM_OPTION_SEQFILE, 'f');
	option_set_shortcut(PROGRAM_OPTION_SEQSTART, 'S');
	option_set_shortcut(PROGRAM_OPTION_SEQLEN, 'l');
	option_set_shortcut(PROGRAM_OPTION_TEST_STEP, 'e');
	option_set_shortcut(PROGRAM_OPTION_THRESHOLD, 'w');
	option_set_shortcut(PROGRAM_OPTION_READFILE, 'F');
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
		ret = option_get_UInt32(PROGRAM_OPTION_TEST_STEP, &Options->TestStep);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_THRESHOLD, &Options->Threshold);

	if (ret == ERR_SUCCESS)
		ret = option_get_String(PROGRAM_OPTION_SEQFILE, &Options->RefSeqFile);

	if (ret == ERR_SUCCESS)
		ret = option_get_UInt32(PROGRAM_OPTION_READ_STRIP, &Options->ReadStrip);

	if (ret == ERR_SUCCESS) {
		char *readFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_READFILE, &readFile);
		if (ret == ERR_SUCCESS && *readFile != '\0') {			
			fprintf(stderr, "Loading reads from %s...\n", readFile);
			ret = input_get_reads(readFile, "sam", &Options->Reads, &Options->ReadCount);
			if (ret == ERR_SUCCESS) {
				ret = paired_reads_init();
				if (ret == ERR_SUCCESS) {
					fprintf(stderr, "Filtering out reads with MAPQ less than %u and stripping %u bases from read ends...\n", Options->ReadPosQuality, Options->ReadStrip);
					input_filter_bad_reads(Options->Reads, &Options->ReadCount, Options->ReadPosQuality, Options->ReadStrip);
					input_sort_reads(Options->Reads, Options->ReadCount);
					if (ret == ERR_SUCCESS)
						ret = paired_reads_insert_array(Options->Reads, Options->ReadCount);

					if (ret == ERR_SUCCESS)
						paired_reads_fix_overlaps(FALSE);

					for (size_t i = 0; i < Options->ReadCount; ++i)
						read_shorten(Options->Reads + i, Options->ReadStrip);
				}
			}
		}
	}

	option_get_Boolean(PROGRAM_OPTION_NO_CONNECT_REFSEQ, &b);
	Options->ParseOptions.ConnectRefSeq = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_CONNECT_READS, &b);
	Options->ParseOptions.ConnectReads = !b;
	option_get_Boolean(PROGRAM_OPTION_NO_BUBBLE_MERGING, &b);
	Options->ParseOptions.MergeBubbles = !b;
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



static ERR_VALUE _process_variant_call(const ASSEMBLY_TASK *Task, const size_t RefSeqStart, const size_t RefSeqEnd, const char *AltSeq, const size_t AltSeqLen, const size_t RSWeight, const size_t ReadWeight, const GEN_ARRAY_size_t *RefIndices, const GEN_ARRAY_size_t *AltIndices, PGEN_ARRAY_VARIANT_CALL VCArray)
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
			boolean nothing = TRUE;

			while (ret == ERR_SUCCESS) {
				switch (*opIt) {
					case '\0':
					case 'M':
						if (!nothing) {
							size_t rLen = tmpRS - refSeq;
							size_t aLen = tmpAltS - altSeq;
							
							if (rLen == 0 || aLen == 0 ||
								((rLen > 1 || aLen > 1) && *refSeq != *altSeq)) {
								--rsPos;
								--refSeq;
								--altSeq;
								++rLen;
								++aLen;
							}

							ret = variant_call_init("1", rsPos, ".", refSeq, rLen, altSeq, aLen, 60, RefIndices, AltIndices, &vc);
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
				ret = _process_variant_call(Task, var->RefSeqStart, var->RefSeqEnd, var->Seq1, var->Seq1Len, var->Seq2Weight, var->Seq1Weight, &var->RefReadIndices, &var->ReadIndices, VCArray);

			if (var->Seq2Weight > realThreshold && var->Seq2Type == kmetRead)
				ret = _process_variant_call(Task, var->RefSeqStart, var->RefSeqEnd, var->Seq2, var->Seq2Len, var->Seq1Weight, var->Seq2Weight, &var->RefReadIndices, &var->ReadIndices, VCArray);
		} else {
			printf("VAR-BACK: %u->%u, %s\n", var->RefSeqStart, var->RefSeqEnd, var->Seq1);
		}

		++var;
	}

	return ret;
}


static size_t _compare_alternate_sequences(const PROGRAM_OPTIONS *Options, PKMER_GRAPH Graph, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	boolean notFound = FALSE;
	POINTER_ARRAY_FOUND_SEQUENCE seqArray;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *alternates[2];
	size_t res = 0;
	size_t alternateLens[2];
	GEN_ARRAY_FOUND_SEQUENCE_VARIANT variants;
	FILE *f = NULL;
	char *directory = NULL;
	char graphName[128];

	if (Graph->TypedEdgeCount[kmetRead] > 0) {
		directory = "succ";
#pragma warning (disable : 4996)											
		sprintf(graphName, "%s" PATH_SEPARATOR "%s" PATH_SEPARATOR "%s-k%u.graph", Options->OutputDirectoryBase, directory, Task->Name, kmer_graph_get_kmer_size(Graph));
		unlink(graphName);
		if (Options->VCFFileHandle != NULL) {
			ret = utils_fopen(graphName, FOPEN_MODE_WRITE, &f);
			if (ret == ERR_SUCCESS) {
				kmer_graph_print(f, Graph);
				utils_fclose(f);
			}
		}
	}

	dym_array_init_FOUND_SEQUENCE_VARIANT(&variants, 140);
	ret = kmer_graph_get_variants(Graph, &variants);
	if (ret == ERR_SUCCESS) {
		alternates[0] = Task->Alternate1;
		alternates[1] = Task->Alternate2;
		alternateLens[0] = Task->Alternate1Length;
		alternateLens[1] = Task->Alternate2Length;
		pointer_array_init_FOUND_SEQUENCE(&seqArray, 140);
		if (Options->MaxPaths > 1) {
			ret = kmer_graph_get_seqs(Graph, Task->Reference, Options->MaxPaths, &seqArray);
			if (ret == ERR_SUCCESS)
				res = pointer_array_size(&seqArray);
		} else res = (Graph->TypedEdgeCount[kmetRead] > 0) ? Options->MaxPaths + 1 : 1;

		if (Task->Name != NULL && Graph->TypedEdgeCount[kmetVariant] > 0) {
			/*
			directory = "succ";
#pragma warning (disable : 4996)	
			sprintf(graphName, "%s" PATH_SEPARATOR "%s" PATH_SEPARATOR "%s-k%u-p%Iu.graph", Options->OutputDirectoryBase, directory, Task->Name, kmer_graph_get_kmer_size(Graph), gen_array_size(&seqArray));
			unlink(graphName);
			if (Options->VCFFileHandle != NULL)
			{
				ret = utils_fopen(graphName, FOPEN_MODE_WRITE, &f);
				if (ret == ERR_SUCCESS) {
					kmer_graph_print(f, Graph);
					utils_fclose(f);
				}
			}
			*/
			if (Options->VCFFileHandle != NULL) {
				PFOUND_SEQUENCE *pseq = seqArray.Data;
				PGEN_ARRAY_VARIANT_CALL vcArray = Options->VCSubArrays + omp_get_thread_num();

				if (res <= Options->MaxPaths) {
					_process_variant_calls(vcArray, Task, &variants, Options->Threshold);
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

static ERR_VALUE _compute_graph(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const PARSE_OPTIONS *ParseOptions, const ASSEMBLY_TASK *Task, PPROGRAM_STATISTICS Statistics)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	uint32_t kmerSize = Options->KMerSize;
	size_t pathCount = 0;

	for (uint32_t i = 0; i < 8; ++i) {
		ret = kmer_graph_create(kmerSize, 2500, 6000, &g);
		if (ret == ERR_SUCCESS) {
			if (Allocator != NULL)
//			g->Allocator = *Allocator;
				ret = kmer_graph_parse_ref_sequence(g, Task->Reference, Task->ReferenceLength, ParseOptions);
			if (ret == ERR_SUCCESS) {
				GEN_ARRAY_KMER_EDGE_PAIR ep;

				dym_array_init_KMER_EDGE_PAIR(&ep, 140);
				ret = kmer_graph_parse_reads(g, Task->Reads, Task->ReadCount, Options->Threshold, ParseOptions, &ep);
				if (ret == ERR_SUCCESS) {
					size_t deletedThings = 0;
					GEN_ARRAY_size_t readIndices;
					GEN_ARRAY_size_t refIndices;

					dym_array_init_size_t(&refIndices, 140);
					dym_array_init_size_t(&readIndices, 140);
					
					g->DeleteEdgeCallback = _on_delete_edge;
					g->DeleteEdgeCallbackContext = &ep;
					kmer_graph_delete_edges_under_threshold(g, 0);
					kmer_graph_delete_trailing_things(g, &deletedThings);
					g->DeleteEdgeCallback = NULL;
					if (g->TypedEdgeCount[kmetRead] > 0) {
						size_t changeCount = 0;
						ret = kmer_graph_connect_reads_by_pairs(g, Options->Threshold, &ep, &changeCount);
						if (ret == ERR_SUCCESS) {
							kmer_graph_compute_weights(g);
							if (ParseOptions->LinearShrink)
								kmer_graph_delete_1to1_vertices(g);

							kmer_graph_delete_edges_under_threshold(g, Options->Threshold);
							kmer_graph_delete_trailing_things(g, &deletedThings);
//							kmer_graph_resolve_read_narrowings(g);
							if (ParseOptions->MergeBubbles) {
								boolean changed = FALSE;

								do {
									changed = FALSE;
									ret = kmer_graph_detect_uncertainities(g, Task->Reference, &changed);
								} while (ret == ERR_SUCCESS && changed);
							
//								kmer_graph_resolve_triangles(g, Options->Threshold);
							}

							if (ret == ERR_SUCCESS)
								pathCount = _compare_alternate_sequences(Options, g, Task, Statistics);
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

		if (ret != ERR_SUCCESS || pathCount <= Options->MaxPaths)
			break;

		kmerSize += 10;
	}

	if (ret != ERR_SUCCESS) {
		++Statistics->FailureCount;
		printf("FAILD: %u\n", ret);
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
				printf("%zu ", i);
				if (RefSeq[i] != Alternate[i]) {
					printf("%zu: The position has SNPs but is not covered in any read (%c %c)\n", i, RefSeq[i], Alternate[i]);
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

ERR_VALUE kmer_freq_distribution(const PROGRAM_OPTIONS *Options, const uint32_t KMerSize, const ONE_READ *Reads, const size_t ReadCount)
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
				read_split(r);
				if (r->ReadSequenceLen >= KMerSize) {
					for (size_t j = 0; j < r->ReadSequenceLen - KMerSize + 1; ++j) {
						char *s = NULL;

						memcpy(kmerString, r->ReadSequence + j, KMerSize*sizeof(char));
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

				for (size_t i = 0; i < maxValue; ++i) {
					if (freqArray[i] > 0)
						fprintf(stdout, "%Iu, %Iu, %lf\n", i, freqArray[i], (double)freqArray[i]*100/ (double)kmerCount);
				}

				utils_free(freqArray);
			}
		}

		utils_free(kmerString);
	}

	for (size_t i = kh_begin(table); i < kh_end(table); ++i) {
		if (kh_exist(table, i))
			utils_free(kh_key(table, i));
	}

	kh_destroy(kc, table);

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
					_readBaseCount += fr->ReadSequenceLen;
					baseCount += fr->ReadSequenceLen;
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
			po.ReadThreshold = Options->Threshold;
			po.RegionStart = RegionStart;
			po.RegionLength = Options->RegionLength;
			ret = _compute_graph(Allocator, Options, &po, &task, &tmpstats);
			assembly_task_finit(&task);
		}
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
			uint32_t coverage = 0;
			const ONE_READ *fr = FilteredReads->Data;
			size_t baseCount = 0;

			for (size_t i = 0; i < gen_array_size(FilteredReads); ++i) {
				baseCount += fr->ReadSequenceLen;
				++fr;
			}

			coverage = (uint32_t)(baseCount / Options->RegionLength);

			PARSE_OPTIONS po = Options->ParseOptions;

			po.ReadThreshold = (coverage >= Options->Threshold * 4) ? coverage / 4 : Options->Threshold;
			po.RegionStart = RegionStart;
			po.RegionLength = Options->RegionLength;
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
	const uint32_t realStep = ((Options->RegionLength + Options->TestStep - 1) / Options->TestStep)*Options->TestStep;
	long done = 0;
	PUTILS_LOOKASIDE el = NULL;
	PUTILS_LOOKASIDE vl = NULL;

	for (uint32_t k = 0; k < realStep; k += Options->TestStep) {
		_init_lookasides(Options->KMerSize, &vl, &el);
#pragma omp parallel for shared(k, Options, Contig, _activeRegionProcessed), private(vl, el)
		for (j = k; j < Contig->Length - Options->RegionLength; j += (int)realStep) {
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
			if (done % (_activeRegionCount / 10000) == 0)
				fprintf(stderr, "%u %%\r", done * 10000 / _activeRegionCount);
		}
	}

	KMER_GRAPH_ALLOCATOR ga;

	for (uint32_t k = 0; k < realStep; k += Options->TestStep) {
		_init_lookasides(Options->KMerSize, &vl, &el);
		ga.VertexAllocatorContext = vl;
		ga.VertexAllocator = _lookaside_vertex_alloc;
		ga.VertexFreer = _lookaside_vertex_free;
		ga.EdgeAllocatorContext = el;
		ga.EdgeAllocator = _lookaside_edge_alloc;
		ga.EdgeFreer = _lookaside_edge_free;
		process_active_region(&ga, Options, Contig->Offset + Contig->Length - Options->RegionLength - k, Contig->Sequence + Contig->Length - Options->RegionLength - k, Options->ReadSubArrays);
		done = utils_atomic_increment(&_activeRegionProcessed);
		if (done % (_activeRegionCount / 10000) == 0)
			fprintf(stderr, "%u %%\r", done * 10000 / _activeRegionCount);
	}

	return;
}


static void repair_reads_in_parallel(const ACTIVE_REGION *Contig, const PROGRAM_OPTIONS *Options)
{
	const long iterations = 2;
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

	utils_allocator_init(omp_get_num_procs());
	omp_init_lock(&_readCoverageLock);
#ifdef _MSC_VER
	uint64_t startTime = GetTickCount64();
#endif
	ret = options_module_init(37);
	if (ret == ERR_SUCCESS) {
		ret = _init_default_values();
		if (ret == ERR_SUCCESS) {
			ret = options_parse_command_line(argc - 2, argv + 2);
			if (ret == ERR_SUCCESS) {
				PROGRAM_OPTIONS po;
				PROGRAM_STATISTICS st;

				memset(&st, 0, sizeof(st));
				ret = _capture_program_options(&po);
				if (ret == ERR_SUCCESS) {
					omp_set_num_threads(po.OMPThreads);
					const char *cmd = argv[1];
					if (strncmp(cmd, "help", sizeof("help")) == 0) {
						options_print_help();
					} else if (strncmp(cmd, "repair", sizeof("repair")) == 0) {
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
							ret = utils_calloc_PUTILS_LOOKASIDE(omp_get_num_procs(), &_vertexLAs);
							if (ret == ERR_SUCCESS)
								ret = utils_calloc_PUTILS_LOOKASIDE(omp_get_num_procs(), &_edgeLAs);

							if (ret == ERR_SUCCESS) {
								ret = utils_calloc_GEN_ARRAY_ONE_READ(omp_get_num_procs(), &po.ReadSubArrays);
								if (ret == ERR_SUCCESS) {
									const size_t numThreads = omp_get_num_procs();
									for (size_t i = 0; i < numThreads; ++i) {
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
												_activeRegionCount += (long)(pa->Length / po.TestStep);

											++pa;
										}

										_activeRegionProcessed = 0;
										pa = regions;
										for (size_t i = 0; i < regionCount; ++i) {
											if (pa->Type == artValid && pa->Length >= po.RegionLength)
												repair_reads_in_parallel(pa, &po);
													
											++pa;
										}

										input_free_regions(regions, regionCount);
									}

									PONE_READ r = po.Reads;
									for (size_t i = 0; i < po.ReadCount; ++i) {
										if (r->ReadSequenceLen > 0 && r->NumberOfFixes * 100 / r->ReadSequenceLen <= po.ParseOptions.ReadMaxErrorRate) {
											read_quality_encode(r);											
											read_write_sam(stdout, r);
											read_quality_decode(r);
										}

										++r;
									}

									utils_free(rsFasta);
									int i = 0;
#pragma omp parallel for shared (po)
									for (i = 0; i < numThreads; ++i)
										dym_array_finit_ONE_READ(po.ReadSubArrays + i);

									utils_free(po.ReadSubArrays);
								}
							}

							utils_free(_edgeLAs);
							utils_free(_vertexLAs);
							fasta_free(&seqFile);
						}
					} else if (strncmp(cmd, "rfreq", sizeof("rfreq")) == 0) {
						kmer_freq_distribution(&po, po.KMerSize, po.Reads, po.ReadCount);
					} else if (strncmp(cmd, "paired", sizeof("paired")) == 0) {
						paired_reads_print(stdout);
					} else if (strncmp(cmd, "call", sizeof("call")) == 0) {
						fprintf(stderr, "K-mer size:                 %u\n", po.KMerSize);
						fprintf(stderr, "Active region length:       %u\n", po.RegionLength);
						fprintf(stderr, "Reference:                  %s\n", po.RefSeqFile);
						fprintf(stderr, "Reads:                      %zu\n", po.ReadCount);
						fprintf(stderr, "Read coverage threshold:    %u\n", po.Threshold);
						fprintf(stderr, "Min. read position quality: %u\n", po.ReadPosQuality);
						fprintf(stderr, "OpenMP thread count:        %i\n", po.OMPThreads);
						fprintf(stderr, "Output VCF file:            %s\n", po.VCFFile);

						paired_reads_fix_overlaps(TRUE);
						if (ret == ERR_SUCCESS) {
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
										po.VCFFileHandle = fopen(po.VCFFile, "w");
										ret = (po.VCFFileHandle != NULL) ? ERR_SUCCESS : ERR_NOT_FOUND;
										if (ret == ERR_SUCCESS)
											dym_array_init_VARIANT_CALL(&po.VCArray, 140);
									}

									if (ret == ERR_SUCCESS) {
										ret = utils_calloc_PUTILS_LOOKASIDE(omp_get_num_procs(), &_vertexLAs);
										if (ret == ERR_SUCCESS)
											ret = utils_calloc_PUTILS_LOOKASIDE(omp_get_num_procs(), &_edgeLAs);
										
										ret = utils_calloc_GEN_ARRAY_VARIANT_CALL(omp_get_num_procs(), &po.VCSubArrays);
										if (ret == ERR_SUCCESS) {
											ret = utils_calloc_GEN_ARRAY_ONE_READ(omp_get_num_procs(), &po.ReadSubArrays);
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
															_activeRegionCount += (long)(pa->Length / po.TestStep);

														++pa;
													}
														
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
												fprintf(stderr, "Merging the results...\n");
												ret = vc_array_merge(&po.VCArray, po.VCSubArrays, numThreads);
												int i = 0;
#pragma omp parallel for shared(po)
												for (i = 0; i <(int) numThreads; ++i) {
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
											if (ret == ERR_SUCCESS) {
												VARIANT_GRAPH vg;

												fprintf(stderr, "Creating variant graph...\n");
												ret = vg_graph_init(po.VCArray.Data, gen_array_size(&po.VCArray), po.Threshold, &vg);
												if (ret == ERR_SUCCESS) {
													ret = vg_graph_add_paired(&vg);
													if (ret == ERR_SUCCESS) {
														vg_graph_color(&vg);
														vg_graph_print(stdout, &vg);
														vg_graph_finalize(&vg);
													}

													vg_graph_finit(&vg);
												}

												vc_array_print(po.VCFFileHandle, po.RefSeqFile, "1", &po.VCArray);
											}

											vc_array_finit(&po.VCArray);
											fclose(po.VCFFileHandle);
										}

									}

									fasta_free(&seqFile);
								}
							} else printf("fix_reads(): %u\n", ret);

							fprintf(stderr, "Read coverage: %lf\n", _readBaseCount / _totalRegionLength );
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
	fprintf(stderr, "Time: %I64u s\n", (endTime - startTime) / 1000);
#endif
	omp_destroy_lock(&_readCoverageLock);

	return ret;
}
