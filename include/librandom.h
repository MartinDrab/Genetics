
#ifndef __LIBRANDOM_H__
#define __LIBRANDOM_H__




#ifdef __cplusplus
extern "C" {
#endif
	unsigned int rand_C(void);
	unsigned int rand_mt(void);
	unsigned int rand_knuthb(void);
	unsigned int rand_minstd(void);
	unsigned int rand_minstd0(void);
	unsigned int rand_ranlux48base(void);
	unsigned int rand_ranlux24base(void);
	unsigned int rand_ranlux48(void);
	unsigned int rand_ranlux24(void);
	void rand_stats(void);
#ifdef __cplusplus
}
#endif 



#endif 
