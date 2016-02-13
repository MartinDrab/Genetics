
#ifndef __GASSM2_H__
#define __GASSM2_H__


#include "err.h"
#include "utils.h"

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
#define PROGRAM_OPTION_HELP								"help"
#define PROGRAM_OPTION_PRINT_RESULTS					"print-results"

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
#define PROGRAM_OPTION_HELP_DESC						"This help"
#define PROGRAM_OPTION_PRINT_RESULTS_DESC				"Print results even if no error occurs"

/************************************************************************/
/*                                                                      */
/************************************************************************/

typedef struct _PROGRAM_OPTIONS {
	boolean Help;
	boolean Test;
	boolean PrintResults;
	uint32_t KMerSize;
	char *ReferenceSequence;
	uint64_t RegionStart;
	uint32_t RegionLength;
	uint32_t TestCount;
	uint32_t ReadLength;
	uint32_t ReadCount;
	char **Reads;
} PROGRAM_OPTIONS, *PPROGRAM_OPTIONS;




#endif 
