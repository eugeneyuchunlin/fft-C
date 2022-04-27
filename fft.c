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

#define W51_REAL_ABS 0.3090169943749474512628694
#define W51_IMG_ABS 0.9510565162951535311819384
#define W52_REAL_ABS 0.8090169943749473402405670
#define W52_IMG_ABS 0.5877852522924732481257593

#ifdef __linux__
#define _GNU_SOURCE
#endif
#include <math.h>

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

#define OMEGA_WITH_2PI_DIV_N(ij, PI_2_DIV_N) OMEGA_WITH_DEG(ij *PI_2_DIV_N)

#define OMEGA(ij, N) OMEGA_WITH_2PI_DIV_N(ij, (M_PI_MUL_2 / (N)))


#define _FFT1_CONTENT   \
    {                   \
        out[0] = in[0]; \
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

#define DECLARE_FFTX_FUNC(num)                                             \
    static inline void FFT##num(complex double in[], complex double out[]) \
        _FFT##num##_CONTENT static inline void _FFT##num(                  \
            complex double in[], complex double out[], int dum)            \
            _FFT##num##_CONTENT static inline void __FFT##num##_MAT_MUL(   \
                complex double in[], complex double out[],                 \
                complex double const *const *dum_mat, int dum_size)        \
                _FFT##num##_CONTENT

DECLARE_FFTX_FUNC(1)
DECLARE_FFTX_FUNC(2)
DECLARE_FFTX_FUNC(3)
DECLARE_FFTX_FUNC(5)

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
        default:                     \
            (func) = __FFT;          \
            break;                   \
        }                            \
    } while (0)

typedef void (*fft_mat_mul_func_t)(complex double *,
                                   complex double *,
                                   complex double const *const *,
                                   int);


void mat_multiplication(complex double in[],
                        complex double out[],
                        complex double const *const *mat,
                        int N)
{
    for (int i = 0; i < N; ++i) {
        complex double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += mat[i][j] * in[j];
        }
        out[i] = sum;
    }
}
typedef void (*mat_mult_func_t)(complex double *,
                                complex double *,
                                complex double **,
                                int);

complex double **create_mat(int N)
{
    complex double *_mat;
    complex double **mat;
    _mat = malloc(sizeof(complex double) * N * N);
    mat = malloc(sizeof(complex double *) * N);

    for (int i = 0; i < N; ++i) {
        mat[i] = _mat + (i * N);
    }

    return mat;
}

static void __init_OMEGA_MAT(complex double **mat, int N)
{
    double pi_2_DIV_N = M_PI_MUL_2 / N;
    double degree;

    int half_n = N >> 1;
    double _cos, _sin;
    for (int i = 1; i <= half_n; ++i) {
        degree = pi_2_DIV_N * i;
        _cos = cos(degree);
        _sin = sin(degree);
        mat[1][i] = _cos - I * _sin;
        mat[1][N - i] = _cos + I * _sin;
    }
    int idx;
    for (int i = 2; i < N; ++i) {
        for (int j = 1; j < N; ++j) {
            idx = i * j % N;
            mat[i][j] = mat[1][idx];  // cexp in slower than cos/sin
        }
    }
}

static void ___init_OMEGA_MAT(complex double **mat,
                              complex double *omegas,
                              int N)
{
    for (int i = 1; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            int idx = i * j % N;
            mat[j][i] = mat[i][j] = omegas[idx];
        }
    }
}

void FFTN(complex double in[], complex double out[], int N)
{
    complex double **mat = create_mat(N);
    __init_OMEGA_MAT(mat, N);
    for (int i = 0; i < N; ++i) {
        mat[0][i] = mat[i][0] = 1;
    }
    mat_multiplication(in, out, (complex double const *const *) mat, N);
    free(mat[0]);
    free(mat);
}



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

int greatest_prime(int number)
{
    int p = 2;
    while (p <= number) {
        if (number % p == 0) {
            number /= p;
        } else {
            ++p;
        }
    }
    return p;
}


void FFT_iter(complex double in[], complex double out[], int size)
{
    // find the greatest prime
    int gp = greatest_prime(size);
    int p, mat_dim;
    complex double **mat = create_mat(gp);
    for (int i = 0; i < gp; ++i) {
        mat[0][i] = mat[i][0] = 1;
    }

    complex double *_in =
        (complex double *) malloc(sizeof(complex double) * size);
    complex double *in_temp =
        (complex double *) malloc(sizeof(complex double) * gp);
    complex double *out_temp =
        (complex double *) malloc(sizeof(complex double) * gp);

    complex double *omegas[10000] = {0};
    p = 2;
    int n_size = size;
    double _sin, _cos;
    while (p <= n_size) {
        if (n_size % p == 0) {
            n_size /= p;
            if (omegas[p] == NULL) {
                omegas[p] =
                    (complex double *) malloc(sizeof(complex double) * p);
                omegas[p][0] = 1;
                for (int i = 1; i < p; ++i) {
                    sincos(i * M_PI_MUL_2 / p, &_sin, &_cos);
                    omegas[p][i] = _cos - I * _sin;
                }
            }
        } else {
            ++p;
        }
    }


    memcpy(_in, in, sizeof(complex double) * size);
    task_queue_t queue;
    init_task_queue(&queue);

    fft_task_t *t = (fft_task_t *) malloc(sizeof(fft_task_t));
    init_fft_task(_in, out, size, t);
    add_task(t, &queue);

    mat_mult_func_t mat_mul;

    fft_task_t *nt;
    while (queue.size > 0) {
        // pop node
        pop_task(&t, &queue);
        if (t->size == 2) {
            FFT2(t->in, t->out);
            free_task(&t);
            continue;
        } else if (t->size == 3) {
            FFT3(t->in, t->out);
            free_task(&t);
            continue;
        } else if (t->size == 5) {
            FFT5(t->in, t->out);
            free_task(&t);
            continue;
        } else {
            p = determine_p(t->size);
            if (p == t->size) {
                if (p != mat_dim) {
                    ___init_OMEGA_MAT(mat, omegas[p], p);
                    mat_dim = p;
                }
                mat_multiplication(t->in, t->out,
                                   (complex double const *const *) mat, p);
                free_task(&t);
                continue;
            }
        }
        if (t->type == DIVIDE) {
            int new_size = t->size / p;
            // reentrant
            add_task(t, &queue);
            t->type = MERGE;
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

                nt = malloc(sizeof(fft_task_t));
                init_fft_task(entries_in, entries_out[i], new_size, nt);
                add_task(nt, &queue);
            }
            t->entries_out = entries_out;
            t->dim_y = p;
            t->dim_x = new_size;
            // printf("List task : \n");
            // list_task(&queue);
            // printf("");
        } else {  // type is MERGE
            mat_mul =
                (mat_mult_func_t) ((p == 2) * (unsigned long) __FFT2_MAT_MUL +
                                   (p == 3) * (unsigned long) __FFT3_MAT_MUL +
                                   (p == 5) * (unsigned long) __FFT5_MAT_MUL +
                                   (p != 2 && p != 3 && p != 5) *
                                       (unsigned long) mat_multiplication);
            if (p != mat_dim) {
                ___init_OMEGA_MAT(mat, omegas[p], p);
                mat_dim = p;
            }

            for (int i = 0; i < t->dim_x; ++i) {
                for (int j = 0; j < p; ++j) {
                    int idx = i * j;
                    double degree = idx * M_PI_MUL_2 / t->size;
                    sincos(degree, &_sin, &_cos);
                    in_temp[j] = t->entries_out[j][i] * (_cos - I * _sin);
                }
                mat_mul(in_temp, out_temp, mat, p);
                for (int j = 0; j < p; ++j) {
                    t->out[i + t->dim_x * j] = out_temp[j];
                }
            }
            free_task(&t);
        }
    }
    free(mat[0]);
    free(mat);
    free(in_temp);
    free(out_temp);
}



static void __FFT(complex double in[],
                  complex double out[],
                  complex double const *const *OMEGA_MAT,
                  int size)
{
    fft_mat_mul_func_t __fft_mat_mul_func;
    if (size == 0) {
        return;
    } else if (size == 1) {
        out[0] = in[0];
    } else if (size == 2) {
        FFT2(in, out);
    } else if (size == 3) {
        FFT3(in, out);
    } else if (size == 5) {
        FFT5(in, out);
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
        DETERMINE_FFT_MUL(new_size, __fft_mat_mul_func);
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
            __fft_mat_mul_func(entries_in, entries_out[i], OMEGA_MAT, new_size);
        }


        complex double *in_temp =
            (complex double *) malloc(sizeof(complex double) * (p << 1));
        complex double *out_temp = in_temp + p;

        DETERMINE_FFT_MUL(p, __fft_mat_mul_func);
        for (int i = 0; i < new_size; ++i) {
            double m_pi_mul_2_sz = M_PI_MUL_2 / size;
            for (int j = 0; j < p; ++j) {
                in_temp[j] = entries_out[j][i] *
                             OMEGA_WITH_2PI_DIV_N(i * j, m_pi_mul_2_sz);
            }
            __fft_mat_mul_func(in_temp, out_temp, OMEGA_MAT, p);

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

void FFT(complex double in[], complex double out[], int size)
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
    double _cos, _sin;
    while (p <= _tmp_size) {
        bool _in = false;
        while (_tmp_size % p == 0) {
            _tmp_size /= p;
            _in = true;
        }
        if (_in) {
            omegas[p] = malloc(sizeof(complex double) * (p + 1));
            int half_p = p >> 1;
            for (int i = 0; i <= half_p; ++i) {
                omegas[p][i] = OMEGA(i, p);
                omegas[p][p - i] = conj(omegas[p][i]);
            }
            omegas[p][0] = 1;
        }
        ++p;
    }

    __FFT(in, out, (complex double const *const *) omegas, size);
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
