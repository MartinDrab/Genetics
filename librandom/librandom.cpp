
#include <random>
#include <iostream>
#include "librandom.h"




static std::mt19937 mt1;
static std::knuth_b knuthb;
static std::minstd_rand0 minstdrand0;
static std::minstd_rand minstdrand;
static std::ranlux24_base ranlux24base;
static std::ranlux48_base ranlux48base;
static std::ranlux24 ranlux24;
static std::ranlux48 ranlux48;
// static std::random_device randomDevice;


extern "C" unsigned int rand_C(void)
{
	return rand();
}

extern "C" unsigned int rand_mt(void)
{
	return mt1() - mt1.min();
}

extern "C" unsigned int rand_knuthb(void)
{
	return knuthb() - knuthb.min();
}

extern "C" unsigned int rand_minstd(void)
{
	return minstdrand() - minstdrand.min();
}


extern "C" unsigned int rand_minstd0(void)
{
	return minstdrand0() - minstdrand0.min();
}


extern "C" unsigned int rand_ranlux48base(void)
{
	return ranlux48base() - ranlux48base.min();
}


extern "C" unsigned int rand_ranlux24base(void)
{
	return ranlux24base() - ranlux24base.min();
}

extern "C" unsigned int rand_ranlux48(void)
{
	return ranlux48() - ranlux48.min();
}


extern "C" unsigned int rand_ranlux24(void)
{
	return ranlux24base() - ranlux24base.min();
}


extern "C" void rand_stats(void)
{
	std::cout << "mt1:         " << mt1.min() << " " << mt1.max() << std::endl;
	std::cout << "knuthb:      " << knuthb.min() << " " << knuthb.max() << std::endl;
	std::cout << "minstdrand0: " << minstdrand0.min() << " " << minstdrand0.max() << std::endl;
	std::cout << "minstdrand:  " << minstdrand.min() << " " << minstdrand.max() << std::endl;
	std::cout << "ranlux24base " << ranlux24base.min() << " " << ranlux24base.max() << std::endl;
	std::cout << "ranlux48base " << ranlux48base.min() << " " << ranlux48base.max() << std::endl;
	std::cout << "ranlux24     " << ranlux24.min() << " " << ranlux24.max() << std::endl;
	std::cout << "ranlux48     " << ranlux48.min() << " " << ranlux48.max() << std::endl;
//	std::cout << "device       " << randomDevice.min() << " " << randomDevice << std::endl;


	return;
}
