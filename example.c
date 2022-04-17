#include <complex.h>  // for complex number
#include <math.h>     // for cos, sin
#include <stdio.h>    // for srandard i/o
#include <stdlib.h>   // for memory mangement

#include "fft.h"

int main()
{
    const int size = 8;

    complex double array[size];
    complex double out[size];

    for (int i = 0; i < size; ++i) {
        array[i] = i;
    }

    FFT(array, out, size);

    for (int i = 0; i < size; ++i) {
        printf("%.3f + %.3fi\n", creal(out[i]), cimag(out[i]));
    }
}
