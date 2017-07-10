
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


ERR_VALUE libcorrect_correct(PONE_READ Reads, size_t Count, const uint32_t Iterations, PLIBRCORRECT_STATISTICS Stats);






#endif
