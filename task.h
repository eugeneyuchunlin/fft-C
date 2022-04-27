#ifndef __TASK_H__
#define __TASK_H__

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

enum TASK_TYPE { DIVIDE, MERGE };

struct __task_t {
    struct list_head node;
    int size;
    complex double *in;
    complex double *out;
    complex double **entries_out;
    int dim_x, dim_y;
    enum TASK_TYPE type;
};
typedef struct __task_t fft_task_t;

struct __task_queue_t {
    struct list_head list;
    int size;
    // pthread_mutex_t _queue_mutex;
};
typedef struct __task_queue_t task_queue_t;

// user-level api

static inline void init_task_queue(task_queue_t *queue)
{
    INIT_LIST_HEAD(&queue->list);
    queue->size = 0;
}

static inline void add_task(fft_task_t *task, task_queue_t *queue)
{
    enqueue(&task->node, &queue->list);
    ++queue->size;
}

static inline void pop_task(fft_task_t **task, task_queue_t *queue)
{
    struct list_head *node;
    dequeue(&node, &queue->list);
    *task = list_entry(node, fft_task_t, node);
    --queue->size;
}

static inline void list_task(task_queue_t *queue)
{
    struct list_head *node;
    fft_task_t *t;
    for (node = queue->list.next; node != &queue->list; node = node->next) {
        t = list_entry(node, fft_task_t, node);
        printf("\tT size = % 2d\n", t->size);
    }
}

static inline void free_task(fft_task_t **t)
{
    if (t) {
        fft_task_t *_t = *t;
        if (_t->entries_out) {
            for (int i = 0; i < _t->dim_y; ++i)
                free(_t->entries_out[i]);
            free(_t->entries_out);
        }
        if (_t->in)
            free(_t->in);
        free(_t);
        *t = NULL;
    }
}

#endif  // __TASK_H__
