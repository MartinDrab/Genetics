

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


#if 1

#define KMER_MAXIMUM_SIZE						KMER_DEBUG_MAXIMUM_SIZE
#define KMER_BYTES(aKMerSize)					KMER_DEBUG_BYTES(aKMerSize)
#define KMER_BYTES_EXTRA(aKMerSize, aExtra)		KMER_DEBUG_BYTES_EXTRA(aKMerSize, aExtra)


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

#else

#define KMER_MAXIMUM_SIZE						KMER_SHORT_MAXIMUM_SIZE
#define KMER_BYTES(aKMerSize)					KMER_SHORT_BYTES(aKMerSize)
#define KMER_BYTES_EXTRA(aKMerSize, aExtra)		KMER_SHORT_BYTES_EXTRA(aKMerSize, aExtra)

#define kmer_get_size				kmer_short_get_size
#define kmer_set_size				kmer_short_set_size
#define kmer_get_base				kmer_short_get_base
#define kmer_set_base				kmer_short_set_base
#define kmer_get_last_base			kmer_short_get_last_base
#define kmer_get_number				kmer_short_get_number
#define kmer_set_number				kmer_short_set_number

#define	kmer_seq_init_raw				kmer_short_seq_init_raw
#define	kmer_seq_init_by_sequence		kmer_short_seq_init_by_sequence
#define	kmer_advance					kmer_short_advance
#define	kmer_back						kmer_short_back
#define	kmer_seq_equal					kmer_short_seq_equal
#define	kmer_print						kmer_short_print
#define	kmer_hash						kmer_short_hash

#endif


#include "kmer-base.h"




#endif
