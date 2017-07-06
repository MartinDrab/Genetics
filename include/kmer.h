

#ifndef __GASSM_KMER_H__
#define __GASSM_KMER_H__


#include <stdint.h>
#include <malloc.h>
#include <stdio.h>
#ifndef _MSC_VER
#include <alloca.h>
#endif
#include "err.h"
#include "utils.h"
#include "kmer-debug.h"
#include "kmer-short.h"


#define kmer_get_size				kmer_debug_get_size
#define kmer_set_size				kmer_debug_set_size
#define kmer_get_base				kmer_debug_get_base
#define kmer_set_base				kmer_debug_set_base
#define kmer_get_last_base			kmer_debug_get_last_base
#define kmer_get_number				kmer_debug_get_number
#define kmer_set_number				kmer_debug_set_number

#define	kmer_seq_init_raw				kmer_debug_seq_init_raw
#define	kmer_seq_init_by_sequence		kmer_debug_seq_init_by_sequence
#define	kmer_advance					kmer_debug_advance
#define	kmer_back						kmer_debug_back
#define	kmer_seq_equal					kmer_debug_seq_equal
#define	kmer_print						kmer_debug_print
#define	kmer_hash						kmer_debug_hash


#include "kmer-base.h"




#endif
