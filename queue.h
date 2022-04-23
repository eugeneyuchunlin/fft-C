// Implementation of queue based on linux-list api
#ifndef __QUEUE_H__
#define __QUEUE_H__

#include "list.h"

#include <pthread.h>
#include <complex.h>

#define complex _Complex // for stupid compiler


#define enqueue(node, head) list_add(node, head)

static inline void dequeue(struct list_head **node, struct list_head *head){

    if(list_empty(head)){
        *node = NULL;
        return;
    }

    *node = head->next;
    list_del(*node);
}

#endif
