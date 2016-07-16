
#ifndef __GASSM2_H__
#define __GASSM2_H__


#include "err.h"
#include "utils.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "found-sequence.h"

/************************************************************************/
/*                   PROGRAM OTPIONS                                    */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE							"kmer-size"
#define PROGRAM_OPTION_SEQUENCE							"seq"
#define PROGRAM_OPTION_SEQFILE							"seq-file"
#define PROGRAM_OPTION_SEQSTART							"seq-start"
#define PROGRAM_OPTION_SEQLEN							"len"
#define PROGRAM_OPTION_TEST								"test"
#define PROGRAM_OPTION_TEST_COUNT						"test-count"
#define PROGRAM_OPTION_TEST_STEP						"test-step"
#define PROGRAM_OPTION_TEST_READ_CYCLES					"read-cycles"
#define PROGRAM_OPTION_READ_COUNT						"read-count"
#define PROGRAM_OPTION_READ_LENGTH						"read-length"
#define PROGRAM_OPTION_HELP								"help"
#define PROGRAM_OPTION_PRINT_RESULTS					"print-results"
#define PROGRAM_OPTION_THRESHOLD						"threshold"
#define PROGRAM_OPTION_READFILE							"read-file"
#define PROGRAM_OPTION_SNP_RATIO						"snp-ratio"
#define PROGRAM_OPTION_INSERT_RATIO						"insert-ratio"
#define PROGRAM_OPTION_DELETE_RATIO						"delete-ratio"
#define PROGRAM_OPTION_TESTFILE							"test-file"
#define PROGRAM_OPTION_OUTPUT_DIRECTORY					"output-directory"
#define PROGRAM_OPTION_VCFFILE							"vcf-file"
#define PROGRAM_OPTION_OMP_THREADS						"omp-threads"
#define PROGRAM_OPTION_NO_READ_FIXING					"no-read-fixing"
#define PROGRAM_OPTION_BASE_QUALITY_THRESHOLD			"base-quality"
#define PROGRAM_OPTION_BASE_QUALITY_MULTIPLIER			"base-multiplier"
#define PROGRAM_OPTION_READ_POS_QUALITY					"pos-quality"

#define PROGRAM_OPTION_ALT1_SEQ							"alternate1-seq"
#define PROGRAM_OPTION_ALT2_SEQ							"alternate2-seq"

#define PROGRAM_OPTION_DISTINCT_PASSES					"distinct-passes"
#define PROGRAM_OPTION_CONNECT_READS					"connect-reads"
#define PROGRAM_OPTION_RESOLVE_BUBBLES					"resolve-bubbles"
#define PROGRAM_OPTION_MERGE_UNBRANCHED					"merge-unbranched"


/************************************************************************/
/*                  OPTION DESCRIPTION                                  */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE_DESC					"Size of a kmer"
#define PROGRAM_OPTION_SEQUENCE_DESC					"Reference sequence"
#define PROGRAM_OPTION_SEQFILE_DESC						"File (FASTA) containing a reference sequence"
#define PROGRAM_OPTION_SEQSTART_DESC					"zero-based offset to the start of the active region."
#define PROGRAM_OPTION_SEQLEN_DESC						"Length of a reference sequence or an active region"
#define PROGRAM_OPTION_TEST_DESC						"Test mode"
#define PROGRAM_OPTION_TEST_COUNT_DESC					"Number of tests"
#define PROGRAM_OPTION_TEST_STEP_DESC					"Determines the number of bases the active region is moved forward when making tests on a reference sequence"
#define PROGRAM_OPTION_HELP_DESC						"This help"
#define PROGRAM_OPTION_PRINT_RESULTS_DESC				"Print results even if no error occurs"
#define PROGRAM_OPTION_TEST_READ_CYCLES_DESC			"Number of times the test reads will be generated and assembled to the graph with a given reference sequence"
#define PROGRAM_OPTION_READ_COUNT_DESC					"Number of reads to generate for a given active region"
#define PROGRAM_OPTION_READ_LENGTH_DESC					"Length of the generated reads"
#define PROGRAM_OPTION_THRESHOLD_DESC					"Weight threshold"
#define PROGRAM_OPTION_READFILE_DESC					"Name of a file (SAM) that contains reads. Valid only for non-test mode."
#define PROGRAM_OPTION_SNP_RATIO_DESC					"A probability that a base in an alternative sequence is chosen randomly rather that copied from the reference sequence."						"snp-ratio"
#define PROGRAM_OPTION_INSERT_RATIO_DESC				"A probability that a base in an alternate sequence starts an insertion"
#define PROGRAM_OPTION_DELETE_RATIO_DESC				"A probability that a base in an alternate sequence is deleted"
#define PROGRAM_OPTION_TESTFILE_DESC					"Load a reference sequence, alternate sequences and reads from a given file"
#define PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC			"Base output directory"
#define PROGRAM_OPTION_VCFFILE_DESC						"VCF file name"

#define PROGRAM_OPTION_ALT1_SEQ_DESC					"alternate1-seq"
#define PROGRAM_OPTION_ALT2_SEQ_DESC					"alternate2-seq"

#define PROGRAM_OPTION_DISTINCT_PASSES_DESC				"distinct-passes"
#define PROGRAM_OPTION_CONNECT_READS_DESC				"connect-reads"
#define PROGRAM_OPTION_RESOLVE_BUBBLES_DESC				"resolve-bubbles"
#define PROGRAM_OPTION_MERGE_UNBRANCHED_DESC			"merge-unbranched"

/************************************************************************/
/*                                                                      */
/************************************************************************/

typedef struct _PROGRAM_OPTIONS {
	boolean Help;
	boolean Test;
	char *TestFile;
	const char *OutputDirectoryBase;
	boolean PrintResults;
	uint32_t KMerSize;
	const char *ReferenceSequence;
	char *RefSeqFile;
	uint64_t RegionStart;
	uint32_t RegionLength;
	uint32_t TestCount;
	uint32_t TestReadCycles;
	uint32_t TestStep;
	uint32_t Threshold;
	uint32_t ReadLength;
	uint32_t ReadCount;
	PONE_READ Reads;
	double SNPRatio;
	double InsertRatio;
	double DeleteRatio;
	char *AltenrateSequence1;
	char *AlternateSequence2;
	boolean MakeDistinctPasses;
	boolean ConnectReads;
	boolean ResolveBubbles;
	boolean MergeUnbranched;
	const char *VCFFile;
	int32_t OMPThreads;
	uint8_t BaseQualityThreshold;
	uint32_t BaseQualityMultiplier;
	uint8_t ReadPosQuality;
	boolean NoFixReads;
	FILE *VCFFileHandle;
	GEN_ARRAY_VARIANT_CALL *VCSubArrays;
	GEN_ARRAY_ONE_READ *ReadSubArrays;
	GEN_ARRAY_VARIANT_CALL VCArray;
} PROGRAM_OPTIONS, *PPROGRAM_OPTIONS;

typedef struct _PROGRAM_STATISTICS {
	uint64_t SuccessCount;
	uint64_t FailureCount;
	uint64_t CannotSucceed;
} PROGRAM_STATISTICS, *PPROGRAM_STATISTICS;



#endif 
