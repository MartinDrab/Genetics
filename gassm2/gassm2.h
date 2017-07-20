
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

#define PROGRAM_OPTION_MISSING_EDGE_PENALTY				"missing-edge-penalty"
#define PROGRAM_OPTION_BACKWARD_REFSEQ_PENALTY			"backward-refseq-penalty"
#define PROGRAM_OPTION_NO_SHORT_VARIANTS				"no-short-variants"



/************************************************************************/
/*                  OPTION DESCRIPTION                                  */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE_DESC					"Size of a kmer"
#define PROGRAM_OPTION_SEQFILE_DESC						"File (FASTA) containing a reference sequence"
#define PROGRAM_OPTION_SEQLEN_DESC						"Length of a reference sequence or an active region"
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
