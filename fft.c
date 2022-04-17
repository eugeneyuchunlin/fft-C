#include "fft.h"

#include <complex.h>
#include <math.h>

void FFT2(complex double in[], complex double out[])
{
    // In this function I want to perform this operation
    //
    //         +--     --+
    // [out] = |  1,  1  | * [in]
    //         |  1, -1  |
    //         +--     --+
    out[0] = in[0] + in[1];
    out[1] = in[0] - in[1];
}

void FFT(complex double in[], complex double out[], int size)
{
    // in this function I wanna seperate the array to 2 parts
    // even part and odd part
    // out = even_part + W_{n}^-k

    if (size == 2) {
        FFT2(in, out);
    } else {
        int half_size = size >> 1;
        // seperate in to two parts
        complex double even[half_size];
        complex double even_result[half_size];

        complex double odd[half_size];
        complex double odd_result[half_size];

        // place the elements to right place
        for (int i = 0; i < size; ++i) {
            even[i] = in[2 * i];
            odd[i] = in[2 * i + 1];
        }

        // FFT seperately
        FFT(even, even_result, half_size);
        FFT(odd, odd_result, half_size);


        for (int i = 0; i < half_size; ++i) {
            complex double w =
                cos(i * 2 * M_PI / size) - I * sin(i * 2 * M_PI / size);
            out[i] = even_result[i] + odd_result[i] * w;
            out[half_size + i] = even_result[i] + (-1) * odd_result[i] * w;
        }
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
