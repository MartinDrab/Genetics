
#ifndef __GASSM2_H__
#define __GASSM2_H__


#include "err.h"
#include "utils.h"

/************************************************************************/
/*                   PROGRAM OTPIONS                                    */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE							"kmer-size"
#define PROGRAM_OPTION_SEQUENCE							"seq"
#define PROGRAM_OPTION_RANDSEQLEN						"randseq-len"
#define PROGRAM_OPTION_TEST								"test"
#define PROGRAM_OPTION_TEST_COUNT						"test-count"
#define PROGRAM_OPTION_HELP								"help"
#define PROGRAM_OPTION_PRINT_RESULTS					"print-results"

/************************************************************************/
/*                  OPTION DESCRIPTION                                  */
/************************************************************************/

#define PROGRAM_OPTION_KMERSIZE_DESC					"Size of a kmer"
#define PROGRAM_OPTION_SEQUENCE_DESC					"Reference sequence"
#define PROGRAM_OPTION_RANDSEQLEN_DESC					"Length of a randomly generated reference sequence"
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
	uint32_t TestSeqLen;
	uint32_t TestCount;
} PROGRAM_OPTIONS, *PPROGRAM_OPTIONS;




#endif 
