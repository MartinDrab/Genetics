
#include <omp.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#include "tinydir.h"
#include "err.h"
#include "utils.h"
#include "kthread.h"
#include "utils-lookaside.h"
#include "file-utils.h"
#include "options.h"
#include "libkmer.h"
#include "input-file.h"
#include "reads.h"
#include "pointer_array.h"
#include "librcorrect.h"
#include "gassm2.h"


UTILS_TYPED_CALLOC_FUNCTION(GEN_ARRAY_ONE_READ)
UTILS_TYPED_CALLOC_FUNCTION(GEN_ARRAY_VARIANT_CALL)


static PUTILS_LOOKASIDE *_vertexLAs;
static PUTILS_LOOKASIDE *_edgeLAs;
static 	EGassm2CommandType _command = gctUnknown;


#define program_option_init(aOptionName, aOptionDescription, aOptionType, aDefaultValue)	\
	{	\
		option_add_##aOptionType((aOptionName), (aDefaultValue));	\
		option_set_description_const((aOptionName), (aOptionDescription));	\
	}	\



static void _init_default_values(int argc, char **argv)
{
	program_option_init(PROGRAM_OPTION_KMERSIZE, PROGRAM_OPTION_KMERSIZE_DESC, UInt32, 21);
	program_option_init(PROGRAM_OPTION_SEQFILE, PROGRAM_OPTION_SEQFILE_DESC, String, "\0");
	program_option_init(PROGRAM_OPTION_SEQLEN, PROGRAM_OPTION_SEQLEN_DESC, UInt32, 2000);
	program_option_init(PROGRAM_OPTION_THRESHOLD, PROGRAM_OPTION_THRESHOLD_DESC, UInt32, 4);
	program_option_init(PROGRAM_OPTION_READFILE, PROGRAM_OPTION_READFILE_DESC, String, "\0");
	program_option_init(PROGRAM_OPTION_VCFFILE, PROGRAM_OPTION_VCFFILE_DESC, String, "-");
	program_option_init(PROGRAM_OPTION_THREADS, PROGRAM_OPTION_THREADS_DESC, Int32, omp_get_num_procs());
	program_option_init(PROGRAM_OPTION_BINOM_THRESHOLD, PROGRAM_OPTION_BINOM_THRESHOLD_DESC, UInt64, 1);
	program_option_init(PROGRAM_OPTION_LOW_QUALITY_VARIANT, PROGRAM_OPTION_LOW_QUALITY_VARIANT_DESC, UInt32, 3);
	program_option_init(PROGRAM_OPTION_OUTPUT_DIRECTORY, PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC, String, "\0");
	program_option_init(PROGRAM_OPTION_READ_POS_QUALITY, PROGRAM_OPTION_READ_POS_QUALITY_DESC, UInt8, 10);
	program_option_init(PROGRAM_OPTION_NO_SHORT_VARIANTS, PROGRAM_OPTION_NO_SHORT_VARIANTS_DESC, Boolean, FALSE);

	option_set_shortcut(PROGRAM_OPTION_KMERSIZE, 'k');
	option_set_shortcut(PROGRAM_OPTION_SEQFILE, 'f');
	option_set_shortcut(PROGRAM_OPTION_SEQLEN, 'l');
	option_set_shortcut(PROGRAM_OPTION_THRESHOLD, 'w');
	option_set_shortcut(PROGRAM_OPTION_READFILE, 'F');
	option_set_shortcut(PROGRAM_OPTION_OUTPUT_DIRECTORY, 'o');
	option_set_shortcut(PROGRAM_OPTION_VCFFILE, 'v');

	_command = gctHelp;
	if (argc > 1) {
		if (strcmp(argv[1], "call") == 0)
			_command = gctCall;
		else if (strcmp(argv[1], "correct") == 0)
			_command = gctCorrect;
	}

	return;
}


static ERR_VALUE _capture_program_options(PPROGRAM_OPTIONS Options)
{
	boolean b = FALSE;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Options, 0, sizeof(PROGRAM_OPTIONS));
	ret = ERR_SUCCESS;
	if (_command == gctCall) {
		ret = option_get_UInt64(PROGRAM_OPTION_BINOM_THRESHOLD, &Options->ParseOptions.BinomThreshold);
		if (ret != ERR_SUCCESS || !in_range(0, 101, Options->ParseOptions.BinomThreshold)) {
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_BINOM_THRESHOLD);
			ret = ERR_INTERNAL_ERROR;
		}
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_UInt32(PROGRAM_OPTION_LOW_QUALITY_VARIANT, &Options->ParseOptions.LQVariant);
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_LOW_QUALITY_VARIANT);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_UInt32(PROGRAM_OPTION_KMERSIZE, &Options->KMerSize);
		if (ret != ERR_SUCCESS || Options->KMerSize > KMER_MAXIMUM_SIZE) {
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_KMERSIZE);
			ret = ERR_INTERNAL_ERROR;
		}
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_Int32(PROGRAM_OPTION_THREADS, &Options->OMPThreads);
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_THREADS);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_UInt8(PROGRAM_OPTION_READ_POS_QUALITY, &Options->ReadPosQuality);
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_READ_POS_QUALITY);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		char *outputDirectory = NULL;

		ret = option_get_String(PROGRAM_OPTION_OUTPUT_DIRECTORY, &outputDirectory);
		if (ret == ERR_SUCCESS) {
			size_t len = strlen(outputDirectory);

			if (outputDirectory[len] == '\\')
				outputDirectory[len] = '\0';

			Options->OutputDirectoryBase = outputDirectory;
		} else fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_OUTPUT_DIRECTORY);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		char *vcfFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_VCFFILE, &vcfFile);
		if (ret == ERR_SUCCESS)
			Options->VCFFile = vcfFile;
	
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_VCFFILE);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_UInt32(PROGRAM_OPTION_SEQLEN, &Options->RegionLength);
		if (ret != ERR_SUCCESS || !in_range(500, 50000, Options->RegionLength)) {
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_SEQLEN);
			ret = ERR_INTERNAL_ERROR;
		}
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_UInt32(PROGRAM_OPTION_THRESHOLD, &Options->Threshold);
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_THRESHOLD);
	}

	if (_command == gctCall && ret == ERR_SUCCESS) {
		ret = option_get_String(PROGRAM_OPTION_SEQFILE, &Options->RefSeqFile);
		if (ret != ERR_SUCCESS)
			fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_SEQFILE);
	}

	if (ret == ERR_SUCCESS) {
		char *readFile = NULL;

		ret = option_get_String(PROGRAM_OPTION_READFILE, &readFile);
		if (ret == ERR_SUCCESS && *readFile != '\0') {			
			fprintf(stderr, "Loading reads from %s...\n", readFile);
			ret = input_get_reads(readFile, "sam", &Options->Reads, &Options->ReadCount);
			if (ret != ERR_SUCCESS)
				fprintf(stderr, "Error during read loading: %u\n", ret);
		} else fprintf(stderr, "Invalid value for the \"%s\" parameter\n", PROGRAM_OPTION_READFILE);
	}

	if (ret == ERR_SUCCESS) {
		if (_command == gctCorrect) {
			if (Options->ReadCount == 0) {
				fprintf(stderr, "No reads found\n");
				ret = ERR_INTERNAL_ERROR;
			}
		} else {
			option_get_Boolean(PROGRAM_OPTION_NO_SHORT_VARIANTS, &b);
			Options->ParseOptions.OptimizeShortVariants = !b;
			Options->ParseOptions.PlotOptions.PlotFlags = GRAPH_PRINT_ALL;
			Options->ParseOptions.ConnectReads = TRUE;
			Options->ParseOptions.ConnectRefSeq = TRUE;
			Options->ParseOptions.HelperVertices = TRUE;
			Options->ParseOptions.LinearShrink = TRUE;
			Options->ParseOptions.MergeBubbles = TRUE;
			Options->ParseOptions.BackwardRefseqPenalty = 8;
			Options->ParseOptions.MissingEdgePenalty = 3;
			Options->TestStep = Options->RegionLength * 3 / 4;
			Options->ReadStrip = 5;
			if (Options->ReadCount == 0) {
				fprintf(stderr, "No reads found\n");
				ret = ERR_INTERNAL_ERROR;
			}

			if (*Options->RefSeqFile == '\0') {
				fprintf(stderr, "No reference sequence file specified\n");
				ret = ERR_INTERNAL_ERROR;
			}
		}
	}

	return ret;
}


static ERR_VALUE _print_graph(const KMER_GRAPH *Graph, const PROGRAM_OPTIONS *Options, const ASSEMBLY_TASK *Task, uint16_t Flag)
{
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char directory[261];
	char graphName[261];

	ret = ERR_SUCCESS;
	if (*Options->OutputDirectoryBase != '\0' && (Options->ParseOptions.PlotOptions.PlotFlags & Flag)) {
		const char *flagNames[] = {
			"ref", "reads", "helper", "long", "th1",
			"conn", "th2", "shrink", "vars", "res"
		};
		uint32_t flagPos = 0;

		while ((Flag & (1 << flagPos)) == 0)
			++flagPos;

		memset(directory, 0, sizeof(directory));
#pragma warning (disable : 4996)											
		snprintf(directory, sizeof(directory), "%s" PATH_SEPARATOR "%s", Options->OutputDirectoryBase, Task->Name);
#ifdef _WIN32
		int err = mkdir(directory);
#else
		int err = mkdir(directory, 0x777);
#endif

		if (err == 0 || errno == EEXIST) {
			if (Graph->TypedEdgeCount[kmetRead] + Graph->TypedEdgeCount[kmetVariant] > 0)
				snprintf(graphName, sizeof(graphName), "%s" PATH_SEPARATOR "%s-k%u-%s.graph", directory, Task->Name, kmer_graph_get_kmer_size(Graph), flagNames[flagPos]);
			else snprintf(graphName, sizeof(graphName), "%s" PATH_SEPARATOR "NR-%s-k%u-%s.graph", directory, Task->Name, kmer_graph_get_kmer_size(Graph), flagNames[flagPos]);

			unlink(graphName);
			if (Options->VCFFileHandle != NULL) {
				ret = utils_fopen(graphName, FOPEN_MODE_WRITE, &f);
				if (ret == ERR_SUCCESS) {
					kmer_graph_print(f, Graph);
					utils_fclose(f);
				}
			}
		}
	}

	return ret;
}


static void _on_delete_edge(const KMER_GRAPH *Graph, const KMER_EDGE *Edge, void *Context)
{
	PGEN_ARRAY_KMER_EDGE_PAIR pairs = (PGEN_ARRAY_KMER_EDGE_PAIR)Context;
	PKMER_EDGE_PAIR p = pairs->Data;

	for (size_t i = 0; i < gen_array_size(pairs); ++i) {
		boolean remove = FALSE;

		remove = (p->U == Edge || p->V == Edge || p->ConnectingEdge == Edge);
		if (!remove) {
			for (size_t j = 0; j < p->EdgeCount; ++j) {
				remove = (p->Edges[j] == Edge);
				if (remove)
					break;
			}
		}
		
		if (remove) {
			if (p->Edges != NULL)
				utils_free(p->Edges);

			memset(p, 0, sizeof(KMER_EDGE_PAIR));
		}

		++p;
	}

	return;
}


static ERR_VALUE _compute_graph(uint32_t KMerSize, const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const PARSE_OPTIONS *ParseOptions, const ASSEMBLY_TASK *Task, PGEN_ARRAY_VARIANT_CALL VCArray)
{
	PKMER_GRAPH g = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	ASSEMBLY_STATE state;

	ret = kmer_graph_create(KMerSize, 2500, 6000, &g);
	if (ret == ERR_SUCCESS) {
		ret = assembly_state_init(g, ParseOptions, Task->Reads, Task->ReadCount, &state);
		if (ret == ERR_SUCCESS) {
			size_t deletedThings = 0;
			GEN_ARRAY_KMER_EDGE_PAIR ep;

			dym_array_init_KMER_EDGE_PAIR(&ep, 140);
			ret = assembly_parse_reference(&state);
			_print_graph(g, Options, Task, GRAPH_PRINT_REFERENCE);
			if (ret == ERR_REF_REPEATS)
				ret = ERR_SUCCESS;
			
			if (ret == ERR_SUCCESS)
				ret = assembly_parse_reads(&state);

			_print_graph(g, Options, Task, GRAPH_PRINT_RAW_READS);
			if (ret == ERR_SUCCESS)
				ret = assembly_add_helper_vertices(&state);
				
			_print_graph(g, Options, Task, GRAPH_PRINT_HELPER);
			if (ret == ERR_SUCCESS)
				ret = assembly_create_long_edges(&state, &ep);

			_print_graph(g, Options, Task, GRAPH_PRINT_LONG_EDGES);
			if (ret == ERR_SUCCESS) {
				g->DeleteEdgeCallback = _on_delete_edge;
				g->DeleteEdgeCallbackContext = &ep;
				kmer_graph_delete_edges_under_threshold(g, 0);
				kmer_graph_delete_trailing_things(g, &deletedThings);
				g->DeleteEdgeCallback = NULL;
			}

			_print_graph(g, Options, Task, GRAPH_PRINT_THRESHOLD_1);
			if (ret == ERR_SUCCESS && g->TypedEdgeCount[kmetRead] > 0) {
				size_t changeCount = 0;

				ret = kmer_graph_connect_reads_by_pairs(g, ParseOptions->ReadThreshold, &ep, &changeCount);
				_print_graph(g, Options, Task, GRAPH_PRINT_CONNECT);
				if (ret == ERR_SUCCESS) {
					kmer_graph_compute_weights(g);
					kmer_graph_delete_edges_under_threshold(g, ParseOptions->ReadThreshold);
					kmer_graph_delete_trailing_things(g, &deletedThings);
				}
				
				_print_graph(g, Options, Task, GRAPH_PRINT_THRESHOLD_2);
				if (ret == ERR_SUCCESS && ParseOptions->LinearShrink)
					kmer_graph_delete_1to1_vertices(g);

				_print_graph(g, Options, Task, GRAPH_PRINT_SHRINK);
				if (ret == ERR_SUCCESS) {
					boolean changed = FALSE;

					do {
						changed = FALSE;
						ret = kmer_graph_detect_variant(g, VCArray, Options->RefSeq.Name, ParseOptions, &changed);
					} while (ret == ERR_SUCCESS && changed);
				}

				if (ret == ERR_SUCCESS)
					ret = assembly_variants_to_edges(&state, VCArray);

				_print_graph(g, Options, Task, GRAPH_PRINT_VARIANTS);
				if (g->TypedEdgeCount[kmetRead] > 0)
					ret = ERR_TOO_COMPLEX;
			}
		
			PKMER_EDGE_PAIR p = ep.Data;

			for (size_t i = 0; i < gen_array_size(&ep); ++i) {
				if (p->Edges != NULL)
					utils_free(p->Edges);

				++p;
			}

			dym_array_finit_KMER_EDGE_PAIR(&ep);
			assembly_state_finit(&state);
		}

		kmer_graph_destroy(g);
	} else printf("kmer_graph_create(): %u\n", ret);

	return ret;
}


static ERR_VALUE _compute_graphs(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const PARSE_OPTIONS *ParseOptions, const ASSEMBLY_TASK *Task, PGEN_ARRAY_VARIANT_CALL VCArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	uint32_t kmerSize = Options->KMerSize;
	GEN_ARRAY_VARIANT_CALL lowerArray;
	GEN_ARRAY_VARIANT_CALL higherArray;
	const uint32_t step = 10;

	dym_array_init_VARIANT_CALL(&lowerArray, 140);
	dym_array_init_VARIANT_CALL(&higherArray, 140);
	for (uint32_t i = 0; i < 8; ++i) {
		if (kmerSize > KMER_MAXIMUM_SIZE)
			break;
		
		ret = _compute_graph(kmerSize, Allocator, Options, ParseOptions, Task, &lowerArray);
		if (ret == ERR_SUCCESS ||
			(ret == ERR_TOO_COMPLEX && kmerSize + step > KMER_MAXIMUM_SIZE)) {
			vc_array_intersection(&lowerArray, &lowerArray, VCArray);			
			break;
		}

		vc_array_clear(&lowerArray);
		vc_array_clear(&higherArray);
		if (ret != ERR_REF_REPEATS && ret != ERR_TOO_COMPLEX)
			break;

		kmerSize += step;
	}

	dym_array_finit_VARIANT_CALL(&higherArray);
	dym_array_finit_VARIANT_CALL(&lowerArray);

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


static double _readBaseCount = 0;
static size_t _totalRegionLength = 0;
static omp_lock_t _readCoverageLock;


ERR_VALUE process_active_region(const KMER_GRAPH_ALLOCATOR *Allocator, const PROGRAM_OPTIONS *Options, const uint64_t RegionStart, const char *RefSeq, PGEN_ARRAY_ONE_READ FilteredReads, PGEN_ARRAY_VARIANT_CALL VCArray)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = input_filter_reads(Options->KMerSize, Options->Reads, Options->ReadCount, RegionStart, Options->RegionLength, FilteredReads);
	if (ret == ERR_SUCCESS) {
		if (gen_array_size(FilteredReads) > 0) {
			char taskName[128];
			ASSEMBLY_TASK task;
			
			PARSE_OPTIONS po = Options->ParseOptions;
			uint32_t coverage = 0;

			memset(po.ReadQualityDistribution, 0, sizeof(po.ReadQualityDistribution));
			{
				const ONE_READ *fr = FilteredReads->Data;
				uint32_t baseCount = 0;

				omp_set_lock(&_readCoverageLock);
				_totalRegionLength += Options->RegionLength;
				for (size_t i = 0; i < gen_array_size(FilteredReads); ++i) {
					_readBaseCount += fr->ReadSequenceLen;
					baseCount += fr->ReadSequenceLen;
					for (size_t j = 0; j < fr->ReadSequenceLen; ++j)
						++po.ReadQualityDistribution[fr->Quality[j]];

					++fr;
				}

				omp_unset_lock(&_readCoverageLock);
				coverage = baseCount / Options->RegionLength;
			}
			
			sprintf(taskName, "%08" PRIu64 "-r%Iu", (uint64_t)RegionStart, gen_array_size(FilteredReads));
			assembly_task_init(&task, RefSeq, Options->RegionLength, NULL, 0, NULL, 0, FilteredReads->Data, gen_array_size(FilteredReads));
			assembly_task_set_name(&task, taskName);
			task.RegionStart = RegionStart;

			po.ReadCoverage = coverage;
			po.ReadThreshold = Options->Threshold;
			po.RegionStart = RegionStart;
			po.RegionLength = Options->RegionLength;
			po.Reference = RefSeq;
			ret = _compute_graphs(Allocator, Options, &po, &task, VCArray);
			assembly_task_finit(&task);
		}

		input_free_filtered_reads(FilteredReads->Data, gen_array_size(FilteredReads));
	}

	dym_array_clear_ONE_READ(FilteredReads);

	return ret;
}


static long _activeRegionCount = 0;
static volatile long _activeRegionProcessed = 0;


static ERR_VALUE _init_lookasides(const uint32_t KmerSize, PUTILS_LOOKASIDE *VA, PUTILS_LOOKASIDE *EA, size_t ThreadIndex)
{
	ERR_VALUE ret = ERR_SUCCESS;

	if (_vertexLAs[ThreadIndex] == NULL) {
		ret = utils_malloc(sizeof(UTILS_LOOKASIDE), _vertexLAs + ThreadIndex);
		if (ret == ERR_SUCCESS) {
			ret = utils_lookaside_init(_vertexLAs[ThreadIndex], sizeof(KMER_VERTEX) + KmerSize, 3000);
			if (ret != ERR_SUCCESS) {
				utils_free(_vertexLAs[ThreadIndex]);
				_vertexLAs[ThreadIndex] = NULL;
			}
		}
	}

	if (ret == ERR_SUCCESS) {
		*VA = _vertexLAs[ThreadIndex];
		if (_edgeLAs[ThreadIndex] == NULL) {
			ret = utils_malloc(sizeof(UTILS_LOOKASIDE), _edgeLAs + ThreadIndex);
			if (ret == ERR_SUCCESS) {
				ret = utils_lookaside_init(_edgeLAs[ThreadIndex], sizeof(KMER_EDGE), 5000);
				if (ret != ERR_SUCCESS) {
					utils_free(_edgeLAs[ThreadIndex]);
					_edgeLAs[ThreadIndex] = NULL;
				}
			}
		}

		if (ret == ERR_SUCCESS)
			*EA = _edgeLAs[ThreadIndex];
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


typedef struct _AR_WRAPPER_CONTEXT{
	const char *Reference;
	uint64_t RegionStart;
	const PROGRAM_OPTIONS *Options;
} AR_WRAPPER_CONTEXT, *PAR_WRAPPER_CONTEXT;

GEN_ARRAY_TYPEDEF(AR_WRAPPER_CONTEXT);
GEN_ARRAY_IMPLEMENTATION(AR_WRAPPER_CONTEXT)


static GEN_ARRAY_AR_WRAPPER_CONTEXT _assemblyTasks;

static void _ar_wrapper(PAR_WRAPPER_CONTEXT Context, long WorkIndex, size_t ThreadNo)
{
	long done = 0;
	PUTILS_LOOKASIDE el = NULL;
	PUTILS_LOOKASIDE vl = NULL;
	KMER_GRAPH_ALLOCATOR ga;
	PAR_WRAPPER_CONTEXT task = Context + WorkIndex;

	_init_lookasides(task->Options->KMerSize, &vl, &el, ThreadNo);
	ga.VertexAllocatorContext = vl;
	ga.VertexAllocator = _lookaside_vertex_alloc;
	ga.VertexFreer = _lookaside_vertex_free;
	ga.EdgeAllocatorContext = el;
	ga.EdgeAllocator = _lookaside_edge_alloc;
	ga.EdgeFreer = _lookaside_edge_free;
	process_active_region(&ga, task->Options, task->RegionStart, task->Reference, task->Options->ReadSubArrays + ThreadNo, task->Options->VCSubArrays + ThreadNo);
	done = utils_atomic_increment(&_activeRegionProcessed);
	if (done % (_activeRegionCount / 10000) == 0)
		fprintf(stderr, "%u %%\r", done * 10000 / _activeRegionCount);

	return;
}


static void process_active_region_in_parallel(const ACTIVE_REGION *Contig, const PROGRAM_OPTIONS *Options)
{
	PUTILS_LOOKASIDE el = NULL;
	PUTILS_LOOKASIDE vl = NULL;
	
	_init_lookasides(Options->KMerSize, &vl, &el, 0);

	for (uint64_t i = 0; i < Contig->Length - Options->RegionLength; i += Options->TestStep) {
		AR_WRAPPER_CONTEXT arCtx;

		arCtx.Options = Options;
		arCtx.Reference = Contig->Sequence + i;
		arCtx.RegionStart = Contig->Offset + i;
		dym_array_push_back_no_alloc_AR_WRAPPER_CONTEXT(&_assemblyTasks, arCtx);
//		kt_for(Options->OMPThreads, _ar_wrapper, &arCtx, MaxWorkIndex);
	}

	KMER_GRAPH_ALLOCATOR ga;

	_init_lookasides(Options->KMerSize, &vl, &el, 0);
	ga.VertexAllocatorContext = vl;
	ga.VertexAllocator = _lookaside_vertex_alloc;
	ga.VertexFreer = _lookaside_vertex_free;
	ga.EdgeAllocatorContext = el;
	ga.EdgeAllocator = _lookaside_edge_alloc;
	ga.EdgeFreer = _lookaside_edge_free;
	process_active_region(&ga, Options, Contig->Offset + Contig->Length - Options->RegionLength, Contig->Sequence + Contig->Length - Options->RegionLength, Options->ReadSubArrays, Options->VCSubArrays);
		
	long done = utils_atomic_increment(&_activeRegionProcessed);
	if (done % (_activeRegionCount / 10000) == 0)
		fprintf(stderr, "%u %%\r", done * 10000 / _activeRegionCount);

	return;
}


int main(int argc, char *argv[])
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

		utils_allocator_init(omp_get_num_procs());
		omp_init_lock(&_readCoverageLock);
		ret = options_module_init(37);
		if (ret == ERR_SUCCESS) {
			_init_default_values(argc, argv);
			if (argc > 1)
				ret = options_parse_command_line(argc - 2, argv + 2);
			
			if (ret == ERR_SUCCESS) {
				PROGRAM_OPTIONS po;

				if (_command != gctHelp)
					ret = _capture_program_options(&po);
				
				if (ret == ERR_SUCCESS) {
					const char *cmd = argv[1];

					omp_set_num_threads(po.OMPThreads);
					if (_command == gctHelp) {
						fprintf(stdout, "Usage: gassm2 call -f <reference.fa> -F <reads.sam> [OPTIONS]\n");
						fprintf(stdout, "Usage: gassm2 correct -F <reads.sam>\n");
						fprintf(stdout, "\nOptions:\n");
						options_print_help();
					} else if (_command == gctCorrect) {
						LIBRCORRECT_STATISTICS correctStats;
						LIBCORRECT_STATE correctState;

						fprintf(stderr, "Correcting reads (%u iterations)...\n", po.Threshold);
						ret = libcorrect_state_init(&correctState, po.Reads, po.ReadCount);
						if (ret == ERR_SUCCESS) {
							fprintf(stderr, "Correcting...\n");
							libcorrect_correct(&correctState);
							ret = libcorrect_correct_stats(&correctState, &correctStats);
							if (ret == ERR_SUCCESS) {
								libcorrect_stats_print(stderr, &correctStats);
								librcorrect_stats_free(&correctStats);
							}

							ret = librcorrect_state_finit(&correctState);
							if (ret == ERR_SUCCESS) {
								for (size_t i = 0; i < po.ReadCount; ++i) {
									if (po.Reads[i].ReadSequenceLen > 0) {
										read_quality_encode(po.Reads + i);
										read_write_sam(stdout, po.Reads + i);
										read_quality_decode(po.Reads + i);
									}
								}
							}
						}
					} else if (_command == gctCall) {
						fprintf(stderr, "K-mer size:                 %u\n", po.KMerSize);
						fprintf(stderr, "Active region length:       %u\n", po.RegionLength);
						fprintf(stderr, "Reference:                  %s\n", po.RefSeqFile);
						fprintf(stderr, "Reads:                      %zu\n", po.ReadCount);
						fprintf(stderr, "Read coverage threshold:    %u\n", po.Threshold);
						fprintf(stderr, "Min. read position quality: %u\n", po.ReadPosQuality);
						fprintf(stderr, "OpenMP thread count:        %i\n", po.OMPThreads);
						fprintf(stderr, "Output VCF file:            %s\n", po.VCFFile);
						fprintf(stderr, "Read end strip:             %u\n", po.ReadStrip);
						fprintf(stderr, "Step size:                  %u\n", po.TestStep);
						fprintf(stderr, "backward penalty:           %u\n", po.ParseOptions.BackwardRefseqPenalty);
						fprintf(stderr, "Missing edge penaly:        %u\n", po.ParseOptions.MissingEdgePenalty);
						fprintf(stderr, "Short variant optimization: %u\n", po.ParseOptions.OptimizeShortVariants);
						fprintf(stderr, "Connect reference:          %u\n", po.ParseOptions.ConnectReads);
						fprintf(stderr, "Connect reads:              %u\n", po.ParseOptions.ConnectRefSeq);
						fprintf(stderr, "Helper vertices:            %u\n", po.ParseOptions.HelperVertices);
						fprintf(stderr, "Linear shrink:              %u\n", po.ParseOptions.LinearShrink);
						fprintf(stderr, "Read threshold:             %u\n", po.Threshold);
						fprintf(stderr, "Binomic threshold:          %" PRIu64 "\n", po.ParseOptions.BinomThreshold);
						fprintf(stderr, "Low q. variant threshold:   %u\n", po.ParseOptions.LQVariant);

						fprintf(stderr, "Filtering out reads with MAPQ less than %u and stripping %u bases from read ends...\n", po.ReadPosQuality, po.ReadStrip);
						BAD_READS_STATISTICS badStats;
						read_set_stats(po.Reads, po.ReadCount, po.ReadPosQuality, &badStats);
						read_set_stats_print(stderr, &badStats);
						input_filter_bad_reads(po.Reads, &po.ReadCount, po.ReadPosQuality, TRUE);
						input_sort_reads(po.Reads, po.ReadCount);

						ret = paired_reads_init();
						if (ret == ERR_SUCCESS) {
							ret = paired_reads_insert_array(po.Reads, po.ReadCount);
							if (ret == ERR_SUCCESS) {
								paired_reads_fix_overlaps(FALSE);

								for (size_t i = 0; i < po.ReadCount; ++i)
									read_shorten(po.Reads + i, po.ReadStrip);

								paired_reads_fix_overlaps(TRUE);
							}

							if (ret == ERR_SUCCESS) {
								if (ret == ERR_SUCCESS) {
									FASTA_FILE seqFile;

									ret = fasta_load(po.RefSeqFile, &seqFile);
									if (ret == ERR_SUCCESS) {
										ret = fasta_read_seq(&seqFile, &po.RefSeq);
										if (ret != ERR_SUCCESS)
											fasta_free(&seqFile);
									}

									if (ret == ERR_SUCCESS) {
										po.VCFFileHandle = NULL;
										if (*po.VCFFile != '\0') {
											if (strcmp(po.VCFFile, "-") != 0) {
												po.VCFFileHandle = fopen(po.VCFFile, "w");
												ret = (po.VCFFileHandle != NULL) ? ERR_SUCCESS : ERR_NOT_FOUND;
											} else po.VCFFileHandle = stdout;

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

													dym_array_init_AR_WRAPPER_CONTEXT(&_assemblyTasks, 140);
													ret = dym_array_reserve_AR_WRAPPER_CONTEXT(&_assemblyTasks, po.RefSeq.Length / po.TestStep);
													size_t regionCount = 0;
													PACTIVE_REGION regions = NULL;

													ret = input_refseq_to_regions(po.RefSeq.Sequence, po.RefSeq.Length, &regions, &regionCount);
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

														kt_for(po.OMPThreads, _ar_wrapper, _assemblyTasks.Data, (long)gen_array_size(&_assemblyTasks));
														input_free_regions(regions, regionCount);
													}

													dym_array_finit_AR_WRAPPER_CONTEXT(&_assemblyTasks);
													fasta_free_seq(&po.RefSeq);
													fprintf(stderr, "Merging the results...\n");
													vc_array_merge(&po.VCArray, po.VCSubArrays, numThreads);

													int i = 0;
#pragma omp parallel for shared(po)
													for (i = 0; i < (int)numThreads; ++i) {
														dym_array_finit_ONE_READ(po.ReadSubArrays + i);
														vc_array_finit(po.VCSubArrays + i);
													}

													utils_free(po.ReadSubArrays);
												}

												utils_free(po.VCSubArrays);
											}

											utils_free(_edgeLAs);
											utils_free(_vertexLAs);

											if (po.VCFFileHandle != NULL && gen_array_size(&po.VCArray) > 0) {
												if (ret == ERR_SUCCESS) {
													VARIANT_GRAPH vg;

													fprintf(stderr, "Creating variant graph...\n");
													ret = vg_graph_init(po.VCArray.Data, gen_array_size(&po.VCArray), po.Threshold, &vg);
													if (ret == ERR_SUCCESS) {
														ret = vg_graph_add_paired(&vg);
														if (ret == ERR_SUCCESS) {
															vg_graph_color(&vg);
//															vg_graph_print(stdout, &vg);
															vg_graph_finalize(&vg);
														}

														vg_graph_finit(&vg);
													}

													vc_array_print(po.VCFFileHandle, po.RefSeqFile, &po.VCArray);
												}

												vc_array_finit(&po.VCArray);
												if (po.VCFFileHandle != NULL && po.VCFFileHandle != stdout)
													fclose(po.VCFFileHandle);
											}

										}

										fasta_free(&seqFile);
									}
								}

							}

							fprintf(stderr, "Read coverage: %lf\n", _readBaseCount / _totalRegionLength);
							paired_reads_finit();
						}
					}
				}
			}

			options_module_finit();
		}

		omp_destroy_lock(&_readCoverageLock);

	return ret;
}
