#include "fft.h"
#include "queue.h"
#include "task.h"

#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define SQRT3 1.7320508075688771931766041
#define SQRT3_DIV_2 0.8660254037844385965883021
#define M_PI_MUL_2 6.2831853071795862319959269
#define ONE_OVER_THREE 0.333333333333333333333333

#define W51_REAL_ABS 0.3090169943749474512628694
#define W51_IMG_ABS 0.9510565162951535311819384
#define W52_REAL_ABS 0.8090169943749473402405670
#define W52_IMG_ABS 0.5877852522924732481257593

#ifndef _GNU_SOURCE
#ifdef __APPLE__
#define sincos __sincos
#else
#define sincos(deg, _s, _c) \
    do {                    \
        *(_s) = sin((deg)); \
        *(_c) = cos((deg)); \
    } while (0)
#endif
#endif

#define OMEGA_WITH_DEG(deg)          \
    __extension__({                  \
        double __cos, __sin;         \
        sincos(deg, &__sin, &__cos); \
        __cos - I *__sin;            \
    })

#define CONJ_OMEGA_WITH_DEG(deg)     \
    __extension__({                  \
        double __cos, __sin;         \
        sincos(deg, &__sin, &__cos); \
        __cos + I *__sin;            \
    })

#define OMEGA_WITH_2PI_DIV_N(ij, PI_2_DIV_N) OMEGA_WITH_DEG(ij *PI_2_DIV_N)
#define CONJ_OMEGA_WITH_2PI_DIV_N(ij, PI_2_DIV_N) \
    CONJ_OMEGA_WITH_DEG(ij *PI_2_DIV_N)

static inline complex double omega_with_2pi_div_n(int ij, double cst)
{
    return OMEGA_WITH_2PI_DIV_N(ij, cst);
}

static inline complex double conj_omega_with_2pi_div_n(int ij, double cst)
{
    return CONJ_OMEGA_WITH_2PI_DIV_N(ij, cst);
}


#define OMEGA(ij, N) OMEGA_WITH_2PI_DIV_N(ij, (M_PI_MUL_2 / (N)))
#define CONJ_OMEGA(ij, N) CONJ_OMEGA_WITH_2PI_DIV_N(ij, (M_PI_MUL_2 / (N)))


// TODO: I need to check two types of a, b
#define SWAP(a, b)          \
    __extension__({         \
        __typeof__(a) temp; \
        temp = a;           \
        a = b;              \
        b = temp;           \
    })


#define _FFT1_CONTENT   \
    {                   \
        out[0] = in[0]; \
    }

#define _INV_FFT1_CONTENT \
    {                     \
        out[0] = in[0];   \
    }


// In this function I want to perform this operation
//
//         +--     --+
// [out] = |  1,  1  | * [in]
//         |  1, -1  |
//         +--     --+
#define _FFT2_CONTENT           \
    {                           \
        out[0] = in[0] + in[1]; \
        out[1] = in[0] - in[1]; \
    }


#define _INV_FFT2_CONTENT               \
    {                                   \
        out[0] = 0.5 * (in[0] + in[1]); \
        out[1] = 0.5 * (in[0] - in[1]); \
    }



// In this function I want to perform FFT3 operation
// and the FFT3 operation is defined as
//           0          1              2
//         +--                                     --+
//  out1   | 1,         1,             1             |  in1
//         |                                         |
//  out2 = | 1, (-1-sqrt(3))/2    (-1+sqrt(3))/2     |  in2
//         |                                         |
//  out3   | 1  (-1+sqrt(3))/2    (-1-sqrt(3))/2     |  in3
//         +--                                     --+
//
#define _FFT3_CONTENT                                                      \
    {                                                                      \
        static const complex double sq3_ps_1_d_2 = -0.5 + SQRT3_DIV_2 * I; \
        static const complex double sq3_mi_1_d_2 = -0.5 - SQRT3_DIV_2 * I; \
        out[0] = in[0] + in[1] + in[2];                                    \
        out[1] = in[0] + in[1] * sq3_mi_1_d_2 + in[2] * sq3_ps_1_d_2;      \
        out[2] = in[0] + in[1] * sq3_ps_1_d_2 + in[2] * sq3_mi_1_d_2;      \
    }

#define _INV_FFT3_CONTENT                                                  \
    {                                                                      \
        static const complex double sq3_ps_1_d_2 = -0.5 + SQRT3_DIV_2 * I; \
        static const complex double sq3_mi_1_d_2 = -0.5 - SQRT3_DIV_2 * I; \
        out[0] = in[0] + in[1] + in[2];                                    \
        out[1] = in[0] + in[1] * sq3_ps_1_d_2 + in[2] * sq3_mi_1_d_2;      \
        out[2] = in[0] + in[1] * sq3_mi_1_d_2 + in[2] * sq3_ps_1_d_2;      \
        out[0] *= ONE_OVER_THREE;                                          \
        out[1] *= ONE_OVER_THREE;                                          \
        out[2] *= ONE_OVER_THREE;                                          \
    }


#define _FFT5_CONTENT                                                      \
    {                                                                      \
        static const complex double w51 = W51_REAL_ABS - I * W51_IMG_ABS;  \
        static const complex double w52 = -W52_REAL_ABS - I * W52_IMG_ABS; \
        static const complex double w53 = -W52_REAL_ABS + I * W52_IMG_ABS; \
        static const complex double w54 = W51_REAL_ABS + I * W51_IMG_ABS;  \
        out[0] = in[0] + in[1] + in[2] + in[3] + in[4];                    \
        out[1] =                                                           \
            in[0] + w51 * in[1] + w52 * in[2] + w53 * in[3] + w54 * in[4]; \
        out[2] =                                                           \
            in[0] + w52 * in[1] + w54 * in[2] + w51 * in[3] + w53 * in[4]; \
        out[3] =                                                           \
            in[0] + w53 * in[1] + w51 * in[2] + w54 * in[3] + w52 * in[4]; \
        out[4] =                                                           \
            in[0] + w54 * in[1] + w53 * in[2] + w52 * in[3] + w51 * in[4]; \
    }

#define _INV_FFT5_CONTENT                                                    \
    {                                                                        \
        static const complex double w51 = (W51_REAL_ABS - I * W51_IMG_ABS);  \
        static const complex double w52 = (-W52_REAL_ABS - I * W52_IMG_ABS); \
        static const complex double w53 = (-W52_REAL_ABS + I * W52_IMG_ABS); \
        static const complex double w54 = (W51_REAL_ABS + I * W51_IMG_ABS);  \
        out[0] = 0.2 * (in[0] + in[1] + in[2] + in[3] + in[4]);              \
        out[1] = 0.2 * (in[0] + w54 * in[1] + w53 * in[2] + w52 * in[3] +    \
                        w51 * in[4]);                                        \
        out[2] = 0.2 * (in[0] + w53 * in[1] + w51 * in[2] + w54 * in[3] +    \
                        w52 * in[4]);                                        \
        out[3] = 0.2 * (in[0] + w52 * in[1] + w54 * in[2] + w51 * in[3] +    \
                        w53 * in[4]);                                        \
        out[4] = 0.2 * (in[0] + w51 * in[1] + w52 * in[2] + w53 * in[3] +    \
                        w54 * in[4]);                                        \
    }



// clang-format off
#define DECLARE_FFTX_FUNC(num)                                             \
    static inline void FFT##num(                                           \
            complex double in[],                                           \
            complex double out[])                                          \
                    _FFT##num##_CONTENT                                    \
                                                                           \
    static inline void _FFT##num(                                          \
            complex double in[],                                           \
            complex double out[], int dum)                                 \
                    _FFT##num##_CONTENT                                    \
                                                                           \
    static inline void __FFT##num##_MAT_MUL(                               \
                complex double in[],                                       \
                complex double out[],                                      \
                complex double const *const *dum_mat,                      \
                int dum_size)                                              \
                    _FFT##num##_CONTENT                                    

#define DECLARE_IFFTX_FUNC(num)                                                    \
    static inline void IFFT##num(                                          \
            complex double in[],                                           \
            complex double out[])                                          \
                    _INV_FFT##num##_CONTENT                                \
                                                                           \
    static inline void _IFFT##num(                                         \
            complex double in[],                                           \
            complex double out[], int dum)                                 \
                    _INV_FFT##num##_CONTENT                                \
                                                                           \
    static inline void __IFFT##num##_MAT_MUL(                              \
                complex double in[],                                       \
                complex double out[],                                      \
                complex double const *const *dum_mat,                      \
                int dum_size)                                              \
                    _INV_FFT##num##_CONTENT
// clang-format on

DECLARE_FFTX_FUNC(1)
DECLARE_FFTX_FUNC(2)
DECLARE_FFTX_FUNC(3)
DECLARE_FFTX_FUNC(5)

DECLARE_IFFTX_FUNC(1)
DECLARE_IFFTX_FUNC(2)
DECLARE_IFFTX_FUNC(3)
DECLARE_IFFTX_FUNC(5)


#define DETERMINE_FFT_MUL(num, func) \
    do {                             \
        switch ((num)) {             \
        case 1:                      \
            (func) = __FFT1_MAT_MUL; \
            break;                   \
        case 2:                      \
            (func) = __FFT2_MAT_MUL; \
            break;                   \
        case 3:                      \
            (func) = __FFT3_MAT_MUL; \
            break;                   \
        case 5:                      \
            (func) = __FFT5_MAT_MUL; \
            break;                   \
        default:                     \
            (func) = __FFT;          \
            break;                   \
        }                            \
    } while (0)

#define DETERMINE_IFFT_MUL(num, func) \
    do {                              \
        switch ((num)) {              \
        case 1:                       \
            (func) = __IFFT1_MAT_MUL; \
            break;                    \
        case 2:                       \
            (func) = __IFFT2_MAT_MUL; \
            break;                    \
        case 3:                       \
            (func) = __IFFT3_MAT_MUL; \
            break;                    \
        case 5:                       \
            (func) = __IFFT5_MAT_MUL; \
            break;                    \
        default:                      \
            (func) = __IFFT;          \
            break;                    \
        }                             \
    } while (0)


typedef void (*fft_native_func_t)(complex double *, complex double *);

typedef void (*fft_mat_mul_func_t)(complex double *,
                                   complex double *,
                                   complex double const *const *,
                                   int);
typedef void (*mat_mult_func_t)(complex double *,
                                complex double *,
                                complex double **,
                                int);

typedef complex double (*omega_with_2pi_div_n_func_t)(int, double);

void init_fft_task(complex double in[],
                   complex double out[],
                   int size,
                   fft_task_t *task)
{
    *task = (fft_task_t){.in = in,
                         .out = out,
                         .size = size,
                         .type = DIVIDE,
                         .entries_out = NULL,
                         .dim_x = 0,
                         .dim_y = 0};
}

int determine_p(int size)
{
    int p = 1;
    if (size % 2 == 0) {
        p = 2;
    } else if (size % 3 == 0) {
        p = 3;
    } else {
        for (int i = 5; i <= size; ++i) {
            if (size % i == 0) {
                p = i;
                break;
            }
        }
    }
    return p;
}

void _FFT_iter(task_queue_t *queue, complex double const *const *omegas)
{
    int p = 2;
    fft_task_t *nt, *t;
    while (queue->size > 0) {
        // pop node
        pop_task(&t, queue);
        if (t->size == 1) {
            __FFT1_MAT_MUL(t->in, t->out, NULL, 0);
            recycle_task(t, queue);
            continue;
        } else if (t->size == 2) {
            __FFT2_MAT_MUL(t->in, t->out, NULL, 0);
            recycle_task(t, queue);
            continue;
        } else if (t->size == 3) {
            __FFT3_MAT_MUL(t->in, t->out, NULL, 0);
            recycle_task(t, queue);
            continue;
        } else if (t->size == 5) {
            __FFT5_MAT_MUL(t->in, t->out, NULL, 0);
            recycle_task(t, queue);
            continue;
        } else {
            p = determine_p(t->size);
            if (p == t->size) {
                for (int i = 0; i < t->size; ++i) {
                    complex double sum = 0;
                    for (int j = 0; j < t->size; ++j) {
                        sum += t->in[j] * omegas[t->size][(i * j) % t->size];
                    }
                    t->out[i] = sum;
                }

                recycle_task(t, queue);
                continue;
            }
        }
        if (t->type == DIVIDE) {
            int new_size = t->size / p;
            // reentrant
            add_task(t, queue);
            t->type = CROSS;
            complex double **entries_out =
                (complex double **) malloc(sizeof(complex double *) * p);
            for (int i = 0; i < p; ++i) {
                // divide
                complex double *entries_in = (complex double *) malloc(
                    sizeof(complex double) * new_size);
                entries_out[i] = (complex double *) malloc(
                    sizeof(complex double) * new_size);

                // place data
                for (int j = 0; j < new_size; ++j) {
                    entries_in[j] = t->in[j * p + i];
                }

                // nt = malloc(sizeof(fft_task_t));
                nt = get_free_task(queue);
                init_fft_task(entries_in, entries_out[i], new_size, nt);
                add_task(nt, queue);
            }
            t->entries_out = entries_out;
            t->dim_y = p;
            t->dim_x = new_size;
            // printf("[DIVIDE]List task : \n");
            // list_task(&queue);
            // printf("");
        } else if (t->type == CROSS) {  // type is CROSS
            complex double **out_temp =
                malloc(sizeof(complex double *) * t->dim_x);
            add_task(t, queue);
            for (int i = 0; i < t->dim_x; ++i) {
                complex double *in_temp = (complex double *) malloc(
                    sizeof(complex double) * t->dim_y);
                out_temp[i] = malloc(sizeof(complex double) * t->dim_y);
                for (int j = 0; j < t->dim_y; ++j) {
                    in_temp[j] = t->entries_out[j][i] * OMEGA(i * j, t->size);
                }
                // printf("t->size = %d\n", t->size);
                nt = get_free_task(queue);
                init_fft_task(in_temp, out_temp[i], t->dim_y, nt);
                add_task(nt, queue);
            }
            for (int i = 0; i < t->dim_y; ++i) {
                free(t->entries_out[i]);
            }
            free(t->entries_out);
            t->entries_out = out_temp;
            t->type = MERGE;
            // printf("[CROSS]List task : \n");
            // list_task(&queue);
            // printf("");
        } else {  // MERGE
            for (int i = 0; i < t->dim_x; ++i) {
                for (int j = 0; j < t->dim_y; ++j) {
                    t->out[i + t->dim_x * j] = t->entries_out[i][j];
                }
                // free(t->entries_out[i]);
            }
            // free(t->entries_out);
            recycle_task(t, queue);
        }
    }
}

void FFT_iter(complex double in[], complex double out[], int size)
{
    int p = 2;
    complex double **omegas =
        (complex double **) malloc(sizeof(complex double *) * 1000);
    memset(omegas, 0, sizeof(complex double *) * 1000);
    int n_size = size;
    while (p <= n_size) {
        if (n_size % p == 0) {
            n_size /= p;
            if (omegas[p] == NULL) {
                omegas[p] =
                    (complex double *) malloc(sizeof(complex double) * (p + 1));
                int half_p = p >> 1;
                for (int i = 0; i <= half_p; ++i) {
                    omegas[p][i] = OMEGA_WITH_DEG(i * M_PI_MUL_2 / p);
                    omegas[p][p - i] = conj(omegas[p][i]);
                }
                omegas[p][0] = 1;
            }
        } else {
            ++p;
        }
    }

    complex double *_in =
        (complex double *) malloc(sizeof(complex double) * size);
    memcpy(_in, in, sizeof(complex double) * size);
    task_queue_t queue;
    init_task_queue(&queue);
    premalloc_tasks(&queue, 10000);

    fft_task_t *t = get_free_task(&queue);
    init_fft_task(_in, out, size, t);
    add_task(t, &queue);

    _FFT_iter(&queue, (complex double const *const *) omegas);
    destroy_task_queue(&queue);
    for (int i = 0; i < 1000; ++i) {
        if (omegas[i] != NULL)
            free(omegas[i]);
    }
    free(omegas);
}



static void __FFT(complex double in[],
                  complex double out[],
                  complex double const *const *OMEGA_MAT,
                  fft_native_func_t functions[],
                  omega_with_2pi_div_n_func_t omega_func,
                  int size)
{
    if (size == 0) {
        return;
    } else if (size == 1) {
        out[0] = in[0];
    } else if (size == 2) {
        functions[0](in, out);
    } else if (size == 3) {
        functions[1](in, out);
    } else if (size == 5) {
        functions[2](in, out);
    } else {
        int p = 1;
        if (size % 2 == 0) {
            p = 2;
        } else if (size % 3 == 0) {
            p = 3;
        } else {
            for (int i = 5; i <= size; ++i) {
                if (size % i == 0) {
                    p = i;
                    break;
                }
            }
            if (p == size) {
                for (int i = 0; i < size; ++i) {
                    complex double sum = 0;
                    for (int j = 0; j < size; ++j) {
                        sum += in[j] * OMEGA_MAT[size][(i * j) % size];
                    }
                    out[i] = sum;
                }
                return;
            }
        }

        int new_size = size / p;
        complex double **entries_out, *__entries_out_content;
        complex double *entries_in = malloc(sizeof(complex double) * new_size);
        __entries_out_content =
            (complex double *) malloc(sizeof(complex double) * p * new_size);
        entries_out = (complex double **) malloc(sizeof(complex double *) * p);
        for (int i = 0; i < p; ++i) {
            entries_out[i] = __entries_out_content + (i * new_size);
        }

        // place the data;
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < new_size; ++j) {
                entries_in[j] = in[j * p + i];
            }
            __FFT(entries_in, entries_out[i], OMEGA_MAT, functions, omega_func,
                  new_size);
        }


        complex double *in_temp =
            (complex double *) malloc(sizeof(complex double) * (p << 1));
        complex double *out_temp = in_temp + p;

        for (int i = 0; i < new_size; ++i) {
            double m_pi_mul_2_sz = M_PI_MUL_2 / size;
            for (int j = 0; j < p; ++j) {
                in_temp[j] =
                    entries_out[j][i] * omega_func(i * j, m_pi_mul_2_sz);
            }
            __FFT(in_temp, out_temp, OMEGA_MAT, functions, omega_func, p);

            for (int j = 0; j < p; ++j) {
                out[i + new_size * j] = out_temp[j];
            }
        }
        free(in_temp);
        free(entries_in);
        free(entries_out);
        free(__entries_out_content);
    }
}

void _FFT(complex double in[], complex double out[], int size, bool inverse)
{
    // from https://en.wikipedia.org/wiki/List_of_prime_numbers
    // There are 1000 prime numbers under 7919

    // omegas is a table to compute the FFT matrix first to reduct the redundant
    // computation. The initialization process would base on the prime divisor
    // number of @size.
    complex double **omegas = malloc(sizeof(complex double *) * 1000);
    memset(omegas, 0, sizeof(complex double *) * 1000);

    // defraction of size to init the omegas table
    int _tmp_size = size;
    int p = 2;
    while (p <= _tmp_size) {
        bool _in = false;
        while (_tmp_size % p == 0) {
            _tmp_size /= p;
            _in = true;
        }
        int _exp_conj = inverse * (-1) + !inverse;
        double __one_over_p = (inverse * (1.0 / p) + !inverse);
        if (_in) {
            omegas[p] = malloc(sizeof(complex double) * (p + 1));
            int half_p = p >> 1;
            for (int i = 0; i <= half_p; ++i) {
                omegas[p][i] = OMEGA(i * _exp_conj, p) * __one_over_p;
                omegas[p][p - i] = conj(omegas[p][i]);
            }
            omegas[p][0] = __one_over_p;
        }
        ++p;
    }

    fft_native_func_t functions[3];
    if (inverse) {
        functions[0] = IFFT2;
        functions[1] = IFFT3;
        functions[2] = IFFT5;
    } else {
        functions[0] = FFT2;
        functions[1] = FFT3;
        functions[2] = FFT5;
    }

    __FFT(in, out, (complex double const *const *) omegas, functions,
          (omega_with_2pi_div_n_func_t) ((unsigned long)
                                                 conj_omega_with_2pi_div_n *
                                             inverse +
                                         (unsigned long) omega_with_2pi_div_n *
                                             !inverse),
          size);

    for (int i = 0; i < 1000; ++i) {
        if (omegas[i] != 0) {
            free(omegas[i]);
        }
    }
    free(omegas);
}

void FFT(complex double in[], complex double out[], int size)
{
    _FFT(in, out, size, false);
}

void IFFT(complex double in[], complex double out[], int size)
{
    _FFT(in, out, size, true);
}


void legacy_FFT(complex double in[], complex double out[], int size)
{
    double pi2_div_size = 2 * M_PI / size;
    for (int k = 0; k < size; ++k) {
        complex double sum = 0;
        double angle = pi2_div_size * k;
        for (int j = 0; j < size; ++j) {
            sum += in[j] * (cos(angle * j) - I * sin(angle * j));
        }
        out[k] = sum;
    }
}
