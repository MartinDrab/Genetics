
#ifndef __GASSM2_H__
#define __GASSM2_H__


#include "err.h"
#include "utils.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "input-file.h"
#include "found-sequence-types.h"
#include "libkmer.h"

/************************************************************************/
/*                   PROGRAM OTPIONS                                    */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE							"kmer-size"
#define PROGRAM_OPTION_SEQFILE							"seq-file"
#define PROGRAM_OPTION_SEQSTART							"seq-start"
#define PROGRAM_OPTION_SEQLEN							"len"
#define PROGRAM_OPTION_TEST_STEP						"test-step"
#define PROGRAM_OPTION_THRESHOLD						"threshold"
#define PROGRAM_OPTION_READFILE							"read-file"
#define PROGRAM_OPTION_OUTPUT_DIRECTORY					"output-directory"
#define PROGRAM_OPTION_VCFFILE							"vcf-file"
#define PROGRAM_OPTION_OMP_THREADS						"omp-threads"
#define PROGRAM_OPTION_READ_POS_QUALITY					"pos-quality"

#define PROGRAM_OPTION_NO_CONNECT_REFSEQ				"no-connect-refseq"
#define PROGRAM_OPTION_NO_CONNECT_READS					"no-connect-reads"
#define PROGRAM_OPTION_NO_BUBBLE_MERGING				"no-bubble-merging"
#define PROGRAM_OPTION_NO_LINEAR_SHRINK					"no-linear-shrink"
#define PROGRAM_OPTION_NO_HELPER_VERTICES				"no-helper-vertices"
#define PROGRAM_OPTION_MISSING_EDGE_PENALTY				"missing-edge-penalty"
#define PROGRAM_OPTION_BACKWARD_REFSEQ_PENALTY			"backward-refseq-penalty"
#define PROGRAM_OPTION_READ_STRIP						"read-strip"
#define PROGRAM_OPTION_NO_SHORT_VARIANTS				"no-short-variants"

#define PROGRAM_OPTION_PLOT_START						"plot-start"
#define PROGRAM_OPTION_PLOT_END							"plot-end"
#define PROGRAM_OPTION_PLOT_STEP						"plot-step"


/************************************************************************/
/*                  OPTION DESCRIPTION                                  */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE_DESC					"Size of a kmer"
#define PROGRAM_OPTION_SEQFILE_DESC						"File (FASTA) containing a reference sequence"
#define PROGRAM_OPTION_SEQSTART_DESC					"zero-based offset to the start of the active region."
#define PROGRAM_OPTION_SEQLEN_DESC						"Length of a reference sequence or an active region"
#define PROGRAM_OPTION_TEST_STEP_DESC					"Determines the number of bases the active region is moved forward when making tests on a reference sequence"
#define PROGRAM_OPTION_THRESHOLD_DESC					"Weight threshold"
#define PROGRAM_OPTION_READFILE_DESC					"Name of a file (SAM) that contains reads. Valid only for non-test mode."
#define PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC			"Base output directory"
#define PROGRAM_OPTION_VCFFILE_DESC						"VCF file name"
#define PROGRAM_OPTION_NO_SHORT_VARIANTS_DESC			"Do not Optimize for short variants"

/************************************************************************/
/*                                                                      */
/************************************************************************/

typedef struct _PROGRAM_OPTIONS {
	const char *OutputDirectoryBase;
	uint32_t KMerSize;
	REFSEQ_DATA RefSeq;
	char *RefSeqFile;
	uint64_t RegionStart;
	uint32_t RegionLength;
	uint32_t TestStep;
	uint32_t Threshold;
	size_t ReadCount;
	PONE_READ Reads;
	const char *VCFFile;
	int32_t OMPThreads;
	uint8_t ReadPosQuality;
	FILE *VCFFileHandle;
	GEN_ARRAY_VARIANT_CALL *VCSubArrays;
	GEN_ARRAY_ONE_READ *ReadSubArrays;
	GEN_ARRAY_VARIANT_CALL VCArray;
	uint32_t ReadStrip;
	PARSE_OPTIONS ParseOptions;
} PROGRAM_OPTIONS, *PPROGRAM_OPTIONS;

typedef struct _PROGRAM_STATISTICS {
	uint64_t SuccessCount;
	uint64_t FailureCount;
	uint64_t CannotSucceed;
} PROGRAM_STATISTICS, *PPROGRAM_STATISTICS;



#endif 
