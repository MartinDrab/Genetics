
#include "err.h"
#include "utils.h"
#include "general-queue.h"




/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/


void gqueue_init(PGQUEUE Queue)
{
	memset(Queue, 0, sizeof(GQUEUE));

	return;
}


void gqueue_insert(PGQUEUE Queue, PGQUEUE_HEADER Item)
{
	Item->Next = NULL;
	if (Queue->Last != NULL)
		Queue->Last->Next = Item;

	Queue->Last = Item;
	if (Queue->First == NULL)
		Queue->First = Queue->Last;

	return;
}


PGQUEUE_HEADER gqueue_remove(PGQUEUE Queue)
{
	PGQUEUE_HEADER ret = Queue->First;

	Queue->First = ret->Next;
	if (Queue->First == NULL)
		Queue->Last = NULL;

	return ret;
}


boolean gqueue_empty(PGQUEUE Queue)
{
	return (Queue->First == NULL && Queue->Last == NULL);
}
