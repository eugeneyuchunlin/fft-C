#ifndef __LIST_H__
#define __LIST_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GNUC__)
#define __LIST_HAVE_TYPEOF 1
#endif

// #ifndef container_of
// #ifdef __LIST_HAVE_TYPEOF
// #define container_of(ptr, type, member)                            \
//     __extension__({                                                \
//         const __typeof__(((type *) 0)->member) *__pmember = (ptr); \
//         (type *) ((char *) __pmember - offsetof(type, member));    \
//     })
// #else
// #define container_of(ptr, type, member) \
//     ((type *) ((char *) (ptr) -offsetof(type, member)))
// #endif

#ifndef container_of
#ifdef __LIST_HAVE_TYPEOF
#define container_of(ptr, type, member)                            \
    __extension__({                                                \
        const __typeof__(((type *) 0)->member) *__pmember = (ptr); \
        (type *) ((char *) __pmember - offsetof(type, member));    \
    })
#else
#define container_of(ptr, type, member) \
    ((type *) ((char *) (ptr) -offsetof(type, member)))
#endif
#endif

struct list_head {
    struct list_head *prev;
    struct list_head *next;
};

#define LIST_HEAD(head) struct list_head head = {&head, &head}


static inline void INIT_LIST_HEAD(struct list_head *head)
{
    head->next = head;
    head->prev = head;
}

static inline void list_add(struct list_head *node, struct list_head *head)
{
    struct list_head *next = head->next;

    next->prev = node;
    node->next = next;
    node->prev = head;
    head->next = node;
}

static inline void list_add_tail(struct list_head *node, struct list_head *head)
{
    struct list_head *prev = head->prev;

    prev->next = node;
    node->prev = prev;
    head->prev = node;
    node->next = head;
}

static inline void list_del(struct list_head *node)
{
    struct list_head *next = node->next;
    struct list_head *prev = node->prev;

    next->prev = prev;
    prev->next = next;

    node->next = NULL;
    node->prev = NULL;
}

static inline void list_del_init(struct list_head *node)
{
    list_del(node);
    INIT_LIST_HEAD(node);
}


static inline int list_empty(const struct list_head *head)
{
    return (head->next == head);
}

static inline int list_is_singular(const struct list_head *head)
{
    return (list_empty(head) && head->prev == head->next);
}

#define list_entry(node, type, member) container_of(node, type, member)


#define list_first_entry(head, type, member) \
    container_of((head)->next, type, member)

#define list_last_entry(head, type, memebr) \
    container_of((head)->prev, type, member)



#ifdef __cplusplus
}
#endif

#endif  //__LIST_H__
