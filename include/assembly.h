
#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__

#include "err.h"
#include "utils.h"
#include "reads.h"


typedef struct _PLOT_OPTIONS {
	uint64_t PlotRefStart;
	uint64_t PlotRefEnd;
	uint16_t PlotStep;
	uint16_t PlotFlags;
} PLOT_OPTIONS, *PPLOT_OPTIONS;

typedef struct _PARSE_OPTIONS {
	uint32_t ReadCoverage;
	uint32_t ReadQualityDistribution[256];
	boolean OptimizeShortVariants;
	boolean ConnectRefSeq;
	boolean ConnectReads;
	boolean MergeBubbles;
	boolean LinearShrink;
	boolean HelperVertices;
	uint32_t MissingEdgePenalty;
	uint32_t BackwardRefseqPenalty;
	uint32_t ReadThreshold;
	const char *Reference;
	uint64_t RegionStart;
	uint32_t RegionLength;
	PLOT_OPTIONS PlotOptions;
} PARSE_OPTIONS, *PPARSE_OPTIONS;

typedef struct _ASSEMBLY_STATE {
	PKMER_GRAPH Graph;
	PONE_READ Reads;
	size_t ReadCount;
	PARSE_OPTIONS ParseOptions;
	PKMER_VERTEX **Paths;
	PKMER_EDGE **EdgePaths;
	uint8_t **FlagPaths;
	size_t *PathLengths;
} ASSEMBLY_STATE, *PASSEMBLY_STATE;


ERR_VALUE assembly_parse_reference(PASSEMBLY_STATE State);
ERR_VALUE assembly_parse_reads(PASSEMBLY_STATE State);
ERR_VALUE assembly_add_helper_vertices(PASSEMBLY_STATE State);
ERR_VALUE assembly_create_long_edges(PASSEMBLY_STATE State, PGEN_ARRAY_KMER_EDGE_PAIR PairArray);
ERR_VALUE assembly_variants_to_edges(PASSEMBLY_STATE State, const GEN_ARRAY_VARIANT_CALL *VCArray);

ERR_VALUE assembly_state_init(PKMER_GRAPH Graph, const PARSE_OPTIONS *ParseOptions, PONE_READ Reads, size_t ReadCount, PASSEMBLY_STATE State);
void assembly_state_finit(PASSEMBLY_STATE State);



#endif 
