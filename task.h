#ifndef __TASK_H__
#define __TASK_H__

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

enum TASK_TYPE { 
    DIVIDE, 
    CROSS,
    MERGE 
};

struct __memory_segment_node_t{
    struct list_head node;
    void *_address;
    size_t bytes;
};
typedef struct __memory_segment_node_t mem_seg_t;

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
    struct list_head _free_nodes;
    struct list_head _allocated_memories;
    int size;
    int __size_of_free_nodes;
};
typedef struct __task_queue_t task_queue_t;

// user-level api

static inline void init_task_queue(task_queue_t *queue)
{
    INIT_LIST_HEAD(&queue->list);
    INIT_LIST_HEAD(&queue->_free_nodes);
    INIT_LIST_HEAD(&queue->_allocated_memories);
    queue->size = 0;
    queue->__size_of_free_nodes = 0;
}

static inline void premalloc_tasks(task_queue_t *queue, int number){
    mem_seg_t * seg = (mem_seg_t *)malloc(sizeof(mem_seg_t) +sizeof(fft_task_t)*number);
    INIT_LIST_HEAD(&seg->node);
    enqueue(&seg->node, &queue->_allocated_memories);

    fft_task_t *__memory = (fft_task_t*)(seg+sizeof(mem_seg_t));
    for(int i = 0; i < number; ++i){
        fft_task_t *t = &__memory[i];
        INIT_LIST_HEAD(&t->node);
        enqueue(&t->node, &queue->_free_nodes);
        ++queue->__size_of_free_nodes;
    }
}

static inline void destroy_task_queue(task_queue_t *queue){
    // free the memory
    struct list_head *node;
    mem_seg_t *sg;
    while(!list_empty(&queue->_allocated_memories)){
        dequeue(&node, &queue->_allocated_memories);
        sg = list_entry(node, mem_seg_t, node);
        free(sg);
    }
}

static inline fft_task_t *get_free_task(task_queue_t *queue){
    fft_task_t *t;
    struct list_head *node;
    if(!queue->__size_of_free_nodes){
        premalloc_tasks(queue, 1000);
    }
    dequeue(&node, &queue->_free_nodes);
    --queue->__size_of_free_nodes;
    t = list_entry(node, fft_task_t, node);
    return t;
}

static inline void recycle_task(fft_task_t *task, task_queue_t *queue){
    enqueue(&task->node, &queue->_free_nodes);
    ++queue->__size_of_free_nodes;
    if (task) {
        if (task->entries_out) {
            for (int i = 0; i < task->dim_y; ++i)
                free(task->entries_out[i]);
            free(task->entries_out);
        }
        if (task->in)
            free(task->in);
    }
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
