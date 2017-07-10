
#ifndef __LIBCORRECT_H__
#define __LIBCORRECT_H__


#include "err.h"
#include "utils.h"


typedef struct _LIBRCORRECT_STATISTICS {
	uint64_t TotalReads;
	uint64_t TotalBases;
	uint64_t TotalRepairs;
	uint64_t ReadsRemoved;
	uint64_t ReadsShortened;
	uint64_t *RepairCountDistribution;
	uint64_t *RepairBasePositionDistribution;
	uint32_t RepairCountDistributionCount;
	int K;
} LIBRCORRECT_STATISTICS, *PLIBRCORRECT_STATISTICS;

typedef struct _LIBCORRECT_STATE {
	PONE_READ Reads;
	size_t Count;
	void *ConvertedReads;
	int K;
} LIBCORRECT_STATE, *PLIBCORRECT_STATE;


ERR_VALUE libcorrect_state_init(PLIBCORRECT_STATE State, const uint32_t NumberOfIterations, PONE_READ Reads, const size_t Count);
ERR_VALUE libcorrect_correct(PLIBCORRECT_STATE State);
ERR_VALUE libcorrect_correct_stats(const LIBCORRECT_STATE *State, PLIBRCORRECT_STATISTICS Stats);
void librcorrect_state_finit(PLIBCORRECT_STATE State);





#endif
