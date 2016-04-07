
#include <random>
#include "librandom.h"




static std::mt19937 mt1;
static std::knuth_b knuthb;
static std::minstd_rand0 minstdrand0;
static std::minstd_rand minstdrand;



extern "C" unsigned int rand_C(void)
{
	return rand();
}

extern "C" unsigned int rand_mt(void)
{
	return mt1();
}

extern "C" unsigned int rand_knuthb(void)
{
	return knuthb();
}

extern "C" unsigned int rand_minstd(void)
{
	return minstdrand();
}

extern "C" unsigned int rand_minstd0(void)
{
	return minstdrand0();
}

