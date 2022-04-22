#include "fft.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SQRT3       1.7320508075688771931766041
#define SQRT3_DIV_2 0.8660254037844385965883021
#define M_PI_MUL_2  6.2831853071795862319959269

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
typedef void(*fftn_function_t)(complex double *, complex double *, int);


void FFT2(complex double in[], complex double out[]) _FFT2_CONTENT
void _FFT2(complex double in[], complex double out[], int dum) _FFT2_CONTENT
void FFT3(complex double in[], complex double out[]) _FFT3_CONTENT
void _FFT3(complex double in[], complex double out[], int dum) _FFT3_CONTENT

void FFTN(complex double in[], complex double out[], int N)
{
    complex double mat[N][N];
    double pi_2_DIV_N = M_PI_MUL_2 / N;
    for (int i = 1; i < N; ++i) {
        for (int j = 1; j < N; ++j) {
            double degree = pi_2_DIV_N * i * j;
            mat[i][j] = cos(degree) - I*sin(degree); // cexp in slower than cos/sin
        }
    }
    for (int i = 0; i < N; ++i) {
        mat[0][i] = mat[i][0] = 1;
    }

    for (int i = 0; i < N; ++i) {
        complex double sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += mat[i][j] * in[j];
        }
        out[i] = sum;
    }
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

        fft_function = (fftn_function_t)((p == size)*(unsigned long)FFTN + 
                    (p != size)*(unsigned long)FFT);
        // invoke FFT
        for (int i = 0; i < p; ++i) {
            fft_function(entries_in[i], entries_out[i], new_size);
        }

        complex double *in_temp =
            (complex double *) malloc(sizeof(complex double) * p);
        complex double *out_temp =
            (complex double *) malloc(sizeof(complex double) * p);
        // collect the answer

        fft_function = (fftn_function_t)((p == 2) * (unsigned long)_FFT2 + 
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
