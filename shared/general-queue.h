
#ifndef __GENERAL_QUEUE_H__
#define __GENERAL_QUEUE_H__


#include "utils.h"



typedef struct _GQUEUE_HEADER {
	struct _GQUEUE_HEADER *Next;
} GQUEUE_HEADER, *PGQUEUE_HEADER;

typedef struct _GQUEUE {
	PGQUEUE_HEADER First;
	PGQUEUE_HEADER Last;
} GQUEUE, *PGQUEUE;



void gqueue_init(PGQUEUE Queue);
void gqueue_insert(PGQUEUE Queue, PGQUEUE_HEADER Item);
PGQUEUE_HEADER gqueue_remove(PGQUEUE Queue);
boolean gqueue_empty(PGQUEUE Queue);




#endif 
