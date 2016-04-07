
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
#ifdef __cplusplus
}
#endif 



#endif 
