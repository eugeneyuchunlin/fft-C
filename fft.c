#include "fft.h"
#include "queue.h"
#include "task.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQRT3 1.7320508075688771931766041
#define SQRT3_DIV_2 0.8660254037844385965883021
#define M_PI_MUL_2 6.2831853071795862319959269

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


static inline void FFT1(complex double in[], complex double out[]) _FFT1_CONTENT
    static inline void _FFT1(complex double in[],
                             complex double out[],
                             int dum) _FFT1_CONTENT

    void FFT2(complex double in[], complex double out[]) _FFT2_CONTENT
    static inline void _FFT2(complex double in[],
                             complex double out[],
                             int dum) _FFT2_CONTENT

    void FFT3(complex double in[], complex double out[]) _FFT3_CONTENT
    static inline void _FFT3(complex double in[],
                             complex double out[],
                             int dum) _FFT3_CONTENT

    void mat_multiplication(complex double in[],
                            complex double out[],
                            complex double **mat,
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

static void __init_omega_mat(complex double **mat, int N)
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

static void ___init_omega_mat(complex double **mat,
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
    __init_omega_mat(mat, N);
    for (int i = 0; i < N; ++i) {
        mat[0][i] = mat[i][0] = 1;
    }
    mat_multiplication(in, out, mat, N);
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
                         .size_of_entries_out = 0};
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

    complex double **omegas =
        (complex double **) malloc(sizeof(complex double *) * (gp + 1));
    memset(omegas, 0, sizeof(complex double *) * (gp + 1));
    p = 2;
    int n_size = size;
    while (p <= n_size) {
        if (n_size % p == 0) {
            n_size /= p;
            if (omegas[p] == NULL) {
                omegas[p] =
                    (complex double *) malloc(sizeof(complex double) * p);
                omegas[p][0] = 1;
                for (int i = 1; i < p; ++i) {
                    omegas[p][i] =
                        cos(i * M_PI_MUL_2 / p) - I * sin(i * M_PI_MUL_2 / p);
                }
            }

            // if(omegas[n_size] == NULL){
            //     omegas[n_size] = (complex double*)malloc(sizeof(complex
            //     double)*n_size); omegas[n_size][0] = 1; for(int i = 1; i <
            //     n_size; ++i){
            //         omegas[n_size][i] = cos(i*M_PI_MUL_2 / n_size) -
            //         I*sin(i*M_PI_MUL_2/ n_size);
            //     }
            // }
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

    fftn_function_t fft_function;

    fft_task_t *nt;
    double _cos, _sin;
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
        } else {
            p = determine_p(t->size);
            if (p == t->size) {
                if (p != mat_dim) {
                    ___init_omega_mat(mat, omegas[p], p);
                    mat_dim = p;
                }
                mat_multiplication(t->in, t->out, mat, p);
                free_task(&t);
                continue;
            }
        }
        int new_size = t->size / p;
        if (t->type == DIVIDE) {
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
            t->size_of_entries_out = p;
            // printf("List task : \n");
            // list_task(&queue);
            // printf("");
        } else {  // type is MERGE
            fft_function =
                (fftn_function_t) ((p == 2) * (unsigned long) _FFT2 +
                                   (p == 3) * (unsigned long) _FFT3 +
                                   (p != 2 && p != 3) * (unsigned long) FFTN);
            if (p != mat_dim) {
                ___init_omega_mat(mat, omegas[p], p);
                mat_dim = p;
            }

            for (int i = 0; i < new_size; ++i) {
                for (int j = 0; j < p; ++j) {
                    int idx = i * j;
                    double degree = idx * M_PI_MUL_2 / t->size;
                    in_temp[j] =
                        t->entries_out[j][i] * (cos(degree) - I * sin(degree));
                }
                mat_multiplication(in_temp, out_temp, mat, p);
                for (int j = 0; j < p; ++j) {
                    t->out[i + new_size * j] = out_temp[j];
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


void FFT(complex double in[], complex double out[], int size)
{
    if (size <= 0) {
        return;
    }
    fftn_function_t fft_function;
    if (size == 1) {
        out[0] = in[0];
    } else if (size == 2) {
        FFT2(in, out);
    } else if (size == 3) {
        FFT3(in, out);
    } else if (size == 5) {
        FFTN(in, out, 5);
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
        }

        int new_size = size / p;
        complex double **entries_in;
        complex double **entries_out;
        entries_in = (complex double **) malloc(sizeof(complex double *) * p);
        entries_out = (complex double **) malloc(sizeof(complex double *) * p);
        for (int i = 0; i < p; ++i) {
            entries_in[i] =
                (complex double *) malloc(sizeof(complex double) * new_size);
            entries_out[i] =
                (complex double *) malloc(sizeof(complex double) * new_size);
        }



        // place the data;
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < new_size; ++j) {
                entries_in[i][j] = in[j * p + i];
            }
        }


        fft_function = (fftn_function_t) ((p == size) * (unsigned long) FFT1 +
                                          (p != size) * (unsigned long) FFT);
        // invoke FFT
        for (int i = 0; i < p; ++i) {
            fft_function(entries_in[i], entries_out[i], new_size);
        }

        complex double *in_temp =
            (complex double *) malloc(sizeof(complex double) * p);
        complex double *out_temp =
            (complex double *) malloc(sizeof(complex double) * p);
        // collect the answer

        fft_function =
            (fftn_function_t) ((p == 2) * (unsigned long) _FFT2 +
                               (p == 3) * (unsigned long) _FFT3 +
                               (p != 2 && p != 3) * (unsigned long) FFTN);
        for (int i = 0; i < new_size; ++i) {
            for (int j = 0; j < p; ++j) {
                double degree = i * j * M_PI_MUL_2 / size;
                complex double w = cos(degree) - I * sin(degree);
                in_temp[j] = entries_out[j][i] * w;
            }
            fft_function(in_temp, out_temp, p);

            for (int j = 0; j < p; ++j) {
                out[i + new_size * j] = out_temp[j];
            }
        }
        free(in_temp);
        free(out_temp);
        for (int i = 0; i < p; ++i) {
            free(entries_in[i]);
            free(entries_out[i]);
        }
        free(entries_in);
        free(entries_out);
    }
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
