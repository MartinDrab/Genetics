
#ifndef __KGRAPH_H__
#define __KGRAPH_H__



/************************************************************************/
/*                      PROGRAM OPTIONS                                 */
/************************************************************************/

#define KGRAPH_OPTION_INPUT_FILE			"input-file"
#define KGRAPH_OPTION_INPUT_TYPE			"input-type"
#define KGRAPH_OPTION_REGION_START			"region-start"
#define KGRAPH_OPTION_REGION_LENGTH			"region-length"
#define KGRAPH_OPTION_UNIT_SIZE				"unit-size"
#define KGRAPH_OPTION_UNIT_OVERLAP			"unit-overlap"
#define KGRAPH_OPTION_MIN_K					"min-k"
#define KGRAPH_OPTION_MAX_K					"max-k"
#define KGRAPH_OPTION_STEP_K				"step-k"
#define KGRAPH_OPTION_SKIP_VERTICES			"skip-vert"

#define KGRAPH_OPTION_HELP					"help"
#define KGRAPH_OPTION_VERBOSE				"verbose"

/************************************************************************/
/*                 PROGRAM OPTION DESCRIPTIONS                          */
/************************************************************************/

#define KGRAPH_OPTION_INPUT_FILE_DESC				"input-file"
#define KGRAPH_OPTION_INPUT_TYPE_DESC				"input-type"
#define KGRAPH_OPTION_REGION_START_DESC				"region-start"
#define KGRAPH_OPTION_REGION_LENGTH_DESC			"region-length"
#define KGRAPH_OPTION_UNIT_SIZE_DESC				"unit-size"
#define KGRAPH_OPTION_UNIT_OVERLAP_DESC				"unit-overlap"
#define KGRAPH_OPTION_MIN_K_DESC					"min-k"
#define KGRAPH_OPTION_MAX_K_DESC					"max-k"
#define KGRAPH_OPTION_STEP_K_DESC					"step-k"
#define KGRAPH_OPTION_SKIP_VERTICES_DESC			"skip-vert"

#define KGRAPH_OPTION_HELP_DESC						"help"
#define KGRAPH_OPTION_VERBOSE_DESC					"verbose"

/************************************************************************/
/*                OPTION RECORD                                         */
/************************************************************************/

typedef struct _KGRAPH_OPTIONS_RECORD {
	char *InputFile;
	char *Inputtype;
	uint64_t Regionstart;
	uint64_t regionLength;
	uint32_t UnitSize;
	float UnitOverlap;
	uint32_t MinK;
	uint32_t MaxK;
	uint32_t StepK;
	boolean SkipVertices;
	boolean Verbose;
	boolean Help;
} KGRAPH_OPTIONS_RECORD, *PKGRAPH_OPTIONS_RECORD;




#endif 
