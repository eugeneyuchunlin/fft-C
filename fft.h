#ifndef __FFT_H__
#define __FFT_H__

#include <complex.h>

/**
 * FFT2() -  A simple FFT for 2 elements
 */
void FFT2(complex double in[], complex double out[]);
void FFT2N(complex double in[], complex double out[], int size);

/**
 * FFT3() - A simple FFT for 3 elements
 */
void FFT3(complex double in[], complex double out[]);

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
