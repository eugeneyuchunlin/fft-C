#ifndef __TASK_H__
#define __TASK_H__

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "queue.h"


enum TASK_TYPE { DIVIDE = 0, CROSS = 1, MERGE = 2 };

struct __memory_segment_node_t {
    struct list_head node;
    void *_address;
    size_t bytes;
};
typedef struct __memory_segment_node_t mem_seg_t;

typedef void *(*thread_function)(void *);
struct __task_t {
    struct list_head node;
    int size;
    complex double *in;
    complex double *out;
    complex double **entries_out;
    int dim_x, dim_y;
    enum TASK_TYPE type;

    int children_cnt;
    _Atomic int finished_children_acnt;
    _Atomic int *parent_finished_acnt;
    thread_function func;
};
typedef struct __task_t fft_task_t;


enum POOL_TYPE { started = 0, shutdown = 1, soft = 2, hard = 3 };
struct __task_queue_t {
    struct list_head list;
    struct list_head _free_nodes;
    struct list_head _allocated_memories;
    int size;
    int __size_of_free_nodes;

    complex double const *const *omegas;

    // TODO: functions:
    // be like FFT1, IFFT1

    // mutexes
    pthread_mutex_t queue_lock;
    pthread_mutex_t mem_bank_lock;
    pthread_cond_t notify;
    enum POOL_TYPE flag;
    pthread_t *threads;
    int thread_count;
};
typedef struct __task_queue_t task_queue_t;


struct __thread_data_t {
    fft_task_t *task;
    task_queue_t *queue;
};
typedef struct __thread_data_t thread_data_t;

// user-level api

static inline void init_task_queue(task_queue_t *queue)
{
    INIT_LIST_HEAD(&queue->list);
    INIT_LIST_HEAD(&queue->_free_nodes);
    INIT_LIST_HEAD(&queue->_allocated_memories);
    queue->size = 0;
    queue->__size_of_free_nodes = 0;
    queue->omegas = NULL;

    queue->queue_lock = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    queue->mem_bank_lock = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
    queue->notify = (pthread_cond_t) PTHREAD_COND_INITIALIZER;
    queue->flag = started;
    queue->threads = NULL;
    queue->thread_count = 0;
}

static inline void premalloc_tasks(task_queue_t *queue, int number)
{
    mem_seg_t *seg =
        (mem_seg_t *) malloc(sizeof(mem_seg_t) + sizeof(fft_task_t) * number);
    if (seg == NULL) {
        perror("Malloc Failed");
        exit(EXIT_FAILURE);
    }
    seg->bytes = sizeof(mem_seg_t) + sizeof(fft_task_t) * number;
    INIT_LIST_HEAD(&seg->node);
    enqueue(&seg->node, &queue->_allocated_memories);
    fft_task_t *__memory = (fft_task_t *) ((char *) seg + sizeof(mem_seg_t));
    for (int i = 0; i < number; ++i) {
        fft_task_t *t = &__memory[i];
        INIT_LIST_HEAD(&t->node);
        enqueue(&t->node, &queue->_free_nodes);
        ++queue->__size_of_free_nodes;
    }
}

static inline void destroy_task_queue(task_queue_t *queue)
{
    // free the memory
    struct list_head *node;
    mem_seg_t *sg;
    while (!list_empty(&queue->_allocated_memories)) {
        dequeue(&node, &queue->_allocated_memories);
        sg = list_entry(node, mem_seg_t, node);
        free(sg);
    }
}

static inline fft_task_t *get_free_task(task_queue_t *queue)
{
    fft_task_t *t;
    struct list_head *node;
    // mem bank critical section begins
    pthread_mutex_lock(&queue->mem_bank_lock);
    if (!queue->__size_of_free_nodes) {
        premalloc_tasks(queue, 1000);
    }
    dequeue(&node, &queue->_free_nodes);
    --queue->__size_of_free_nodes;
    // mem bank critical section end
    pthread_mutex_unlock(&queue->mem_bank_lock);

    t = list_entry(node, fft_task_t, node);
    return t;
}

static inline void recycle_task(fft_task_t *task, task_queue_t *queue)
{
    INIT_LIST_HEAD(&task->node);

    // mem bank critical section begins
    pthread_mutex_lock(&queue->mem_bank_lock);
    enqueue(&task->node, &queue->_free_nodes);
    ++queue->__size_of_free_nodes;
    if (task) {
        int y;
        y = (task->type == MERGE) * (task->dim_x) +
            (task->type != MERGE) * (task->dim_y);
        if (task->entries_out) {
            for (int i = 0; i < y; ++i)
                free(task->entries_out[i]);
            free(task->entries_out);
        }
        if (task->in)
            free(task->in);
    }
    // mem bank critical section end
    pthread_mutex_unlock(&queue->mem_bank_lock);
}

static inline void add_task(fft_task_t *task, task_queue_t *queue)
{
    // critical section begin
    pthread_mutex_lock(&queue->queue_lock);
    enqueue(&task->node, &queue->list);
    ++queue->size;
    // critical section end
    pthread_mutex_unlock(&queue->queue_lock);
    // notify
    pthread_cond_signal(&queue->notify);
}

static inline void pop_task(fft_task_t **task, task_queue_t *queue)
{
    // critical section begin
    pthread_mutex_lock(&queue->queue_lock);
    struct list_head *node;
    dequeue(&node, &queue->list);
    *task = list_entry(node, fft_task_t, node);
    --queue->size;
    // critical section end
    pthread_mutex_unlock(&queue->queue_lock);
}

static inline bool has_ready_task(task_queue_t *queue)
{
    struct list_head *node;

    // critical section begin
    // pthread_mutex_lock(&queue->queue_lock);
    fft_task_t *t;
    for (node = queue->list.next; node != &queue->list; node = node->next) {
        t = list_entry(node, fft_task_t, node);
        if (t->finished_children_acnt == t->children_cnt) {
            return true;
        }
    }
    return false;
}

static inline void get_ready_task(fft_task_t **task, task_queue_t *queue)
{
    struct list_head *node;

    // critical section begin
    // pthread_mutex_lock(&queue->queue_lock);
    fft_task_t *t;
    for (node = queue->list.next; node != &queue->list; node = node->next) {
        t = list_entry(node, fft_task_t, node);
        if (t->finished_children_acnt == t->children_cnt) {
            *task = t;
            list_del(node);
            --queue->size;
            return;
        }
    }
    *task = NULL;
    // critical section end
    // pthread_mutex_unlock(&queue->queue_lock);
}

static inline void list_task(task_queue_t *queue)
{
    struct list_head *node;
    fft_task_t *t;
    for (node = queue->list.next; node != &queue->list; node = node->next) {
        t = list_entry(node, fft_task_t, node);
        printf(
            "\tT size = % 2d, type = %d, acnt = % 2d, parent = %p, "
            "finished_children_acnt_address = %p\n",
            t->size, t->type, t->finished_children_acnt,
            t->parent_finished_acnt, &t->finished_children_acnt);
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
