// Implementation of queue based on linux-list api
#ifndef __QUEUE_H__
#define __QUEUE_H__

#include "list.h"

#include <complex.h>
#include <pthread.h>

#define complex _Complex  // for stupid compiler


#define enqueue(node, head) list_add(node, head)

static inline void dequeue(struct list_head **node, struct list_head *head)
{
    if (list_empty(head)) {
        *node = NULL;
        return;
    }

    *node = head->next;
    list_del(*node);
}

static inline void dequeue_no_del(struct list_head **node,
                                  struct list_head *head)
{
    if (list_empty(head)) {
        *node = NULL;
        return;
    }

    *node = head->next;
}

#endif
