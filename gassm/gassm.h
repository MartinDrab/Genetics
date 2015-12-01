
#ifndef __GASSM_H__
#define __GASSM_H__


/************************************************************************/
/*                 PROGRAM OPTION NAMES                                 */
/************************************************************************/

#define GASSM_OPTION_KMER_SIZE			"kmer-size"
#define GASSM_OPTION_REFSEQ_INPUT_FILE	"refseq-input-file"
#define GASSM_OPTION_REFSEQ_INPUT_TYPE	"refseq-input-type"
#define GASSM_OPTION_REFSEQ_SKIP_VERT	"refseq-skip-vertices"
#define GASSM_OPTION_READS_INPUT_FILE	"reads-input-file"
#define GASSM_OPTION_READS_INPUT_TYPE	"reads-input-type"
#define GASSM_OPTION_READS_SKIP_VERT	"reads-skip-vertices"
#define GASSM_OPTION_OUTPUT_FILE		"output-file"
#define GASSM_OPTION_ACTION				"action"
#define GASSM_OPTION_REGION_START		"region-start"
#define GASSM_OPTION_REGION_LENGTH		"region-length"
#define GASSM_OPTION_WEIGHT_THRESHOLD	"weight-threshold"
#define GASSM_OPTION_NUM_PATHS			"num-paths"

#define GASSM_OPTION_HELP				"help"
#define GASSM_OPTION_VERBOSE			"verbose"


/************************************************************************/
/*                PROGRAM OPTIONS DESCRIPTION                           */
/************************************************************************/

#define GASSM_OPTION_KMER_SIZE_DESC				"Size of K-mer to use"
#define GASSM_OPTION_REFSEQ_INPUT_FIL_DESC	"Name of a file with the reference sequence"
#define GASSM_OPTION_REFSEQ_INPUT_TYPE_DESC	"Format of the sref. sequence input file"
#define GASSM_OPTION_REFSEQ_SKIP_VERT_DESC	"When constructing a graph for the ref. sequence, skip vertices to avoid cycles"
#define GASSM_OPTION_READS_INPUT_FILE_DESC	"Name of a file that stores reads to be processed"
#define GASSM_OPTION_READS_INPUT_TYPE_DESC	"Format of the input file with reads"
#define GASSM_OPTION_READS_SKIP_VERT_DESC	"When adding reads to the graph, attempt to skip vertices to avoid cycles"
#define GASSM_OPTION_OUTPUT_FILE_DESC		"Specifies the place where the output should be stored"
#define GASSM_OPTION_ACTION_DESC			"Action to take"
#define GASSM_OPTION_REGION_START_DESC		"Offset to the start of the active region (relative to the start of the reference sequence)"
#define GASSM_OPTION_REGION_LENGTH_DESC		"Length of the active region"
#define GASSM_OPTION_WEIGHT_THRESHOLD_DESC	"Edges with weights below this value will be deleted from the graph"
#define GASSM_OPTION_NUM_PATHS_DESC			"Number of best paths to select as variants candidates"

#define GASSM_OPTION_HELP_DESC					"See this help"
#define GASSM_OPTION_VERBOSE_DESC				"Produce verbose output"

#endif 
