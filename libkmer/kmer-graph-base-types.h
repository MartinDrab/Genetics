
#ifndef __KMER_GRAPH_BASE_TYPES_H__
#define __KMER_GRAPH_BASE_TYPES_H__



typedef enum _EKMerVertexType {
	kmvtRefSeqStart,
	kmvtRefSeqMiddle,
	kmvtRefSeqEnd,
	kmvtRead,
	kmvtMax
} EKMerVertexType, PEKMerVertexType;


typedef enum _EKMerEdgeType {
	kmetReference,
	kmetRead,
	kmetVariant,
	kmetNone,
	kmetMax,
} EKMerEdgeType, *PEKMerEdgeType;







#endif
