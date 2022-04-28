#ifndef __FFT_H__
#define __FFT_H__

#undef _GNU_SOURCE
#if defined(__linux__) || defined(linux) || defined(__linux) || defined(__gnu_linux__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#include <math.h>
#include <complex.h>
#include "task.h"


typedef void (*fftn_function_t)(complex double *, complex double *, int);

void init_fft_task(complex double in[],
                   complex double out[],
                   int size,
                   fft_task_t *task);



void FFT_iter(complex double in[], complex double out[], int size);

void FFT2N(complex double in[], complex double out[], int size);
/**
 * FFT() - A general FFT for 2^q elements of data
 */
void FFT(complex double in[], complex double out[], int size);


void FFT3N(complex double in[], complex double out[], int size);

void FFTN(complex double in[], complex double out[], int N);

/**
 * legacy_FFT() - A legacy fft which performs on O(N^2) complexity
 */
void legacy_FFT(complex double in[], complex double out[], int size);

#endif
