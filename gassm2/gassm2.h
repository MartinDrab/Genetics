
#ifndef __GASSM2_H__
#define __GASSM2_H__


#include "err.h"
#include "utils.h"
#include "reads.h"
#include "gen_dym_array.h"
#include "input-file.h"
#include "variant-types.h"
#include "libkmer.h"

/************************************************************************/
/*                   PROGRAM OTPIONS                                    */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE							"kmer-size"
#define PROGRAM_OPTION_SEQFILE							"seq-file"
#define PROGRAM_OPTION_SEQLEN							"len"
#define PROGRAM_OPTION_THRESHOLD						"threshold"
#define PROGRAM_OPTION_READFILE							"read-file"
#define PROGRAM_OPTION_OUTPUT_DIRECTORY					"output-directory"
#define PROGRAM_OPTION_VCFFILE							"vcf-file"
#define PROGRAM_OPTION_OMP_THREADS						"omp-threads"
#define PROGRAM_OPTION_READ_POS_QUALITY					"pos-quality"
#define PROGRAM_OPTION_LOW_QUALITY_VARIANT				"low-quality-variant"
#define PROGRAM_OPTION_BINOM_THRESHOLD					"binom-threshold"
#define PROGRAM_OPTION_NO_SHORT_VARIANTS				"no-short-variants"



/************************************************************************/
/*                  OPTION DESCRIPTION                                  */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE_DESC					"Initial k-mer size (21...100)"
#define PROGRAM_OPTION_SEQFILE_DESC						"File (FASTA) containing the reference sequence"
#define PROGRAM_OPTION_SEQLEN_DESC						"Length of a reference sequence or an active region (500..50000)"
#define PROGRAM_OPTION_THRESHOLD_DESC					"Global threshold"
#define PROGRAM_OPTION_READFILE_DESC					"Name of a file (SAM) that contains reads. Valid only for non-test mode."
#define PROGRAM_OPTION_OUTPUT_DIRECTORY_DESC			"Directory for debug outputs"
#define PROGRAM_OPTION_VCFFILE_DESC						"VCF file name"
#define PROGRAM_OPTION_NO_SHORT_VARIANTS_DESC			"Do not optimize for short variants"
#define PROGRAM_OPTION_LOW_QUALITY_VARIANT_DESC			"Maximum number of supporting reads for low quality variants"
#define PROGRAM_OPTION_BINOM_THRESHOLD_DESC				"Binomial threshold (0..100)"
#define PROGRAM_OPTION_READ_POS_QUALITY_DESC			"Minimal mapping quality of accepted reads"

/************************************************************************/
/*                                                                      */
/************************************************************************/

#define GRAPH_PRINT_REFERENCE			0x1
#define GRAPH_PRINT_RAW_READS			0x2
#define GRAPH_PRINT_HELPER				0x4
#define GRAPH_PRINT_LONG_EDGES			0X8
#define GRAPH_PRINT_THRESHOLD_1			0x10
#define GRAPH_PRINT_CONNECT				0x20
#define GRAPH_PRINT_THRESHOLD_2			0x40
#define GRAPH_PRINT_SHRINK				0x80
#define GRAPH_PRINT_VARIANTS			0x100
#define GRAPH_PRINT_RESULT				0x200

#define GRAPH_PRINT_ALL	(														\
	GRAPH_PRINT_LONG_EDGES | GRAPH_PRINT_THRESHOLD_1 | GRAPH_PRINT_CONNECT |	\
	GRAPH_PRINT_THRESHOLD_2 | GRAPH_PRINT_SHRINK | GRAPH_PRINT_VARIANTS)		\


/************************************************************************/
/*                                                                      */
/************************************************************************/


typedef enum _EGassm2CommandType {
	gctUnknown,
	gctHelp,
	gctCall,
	gctCorrect,
} EGassm2CommandType, *PEGassm2CommandType;

typedef struct _PROGRAM_OPTIONS {
	const char *OutputDirectoryBase;
	uint32_t KMerSize;
	REFSEQ_DATA RefSeq;
	char *RefSeqFile;
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



#endif 
