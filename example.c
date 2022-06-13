#include <assert.h>
#include <complex.h>  // for complex number
#include <math.h>     // for cos, sin
#include <stdio.h>    // for srandard i/o
#include <stdlib.h>   // for memory mangement
#include <time.h>
#include "fft.h"



int main()
{
    const int size = 12;
    clock_t t1, t2;
    complex double *array =
        (complex double *) malloc(sizeof(complex double) * size);
    complex double *out =
        (complex double *) malloc(sizeof(complex double) * size);

    for (int i = 0; i < size; ++i) {
        array[i] = i;
    }
    t1 = clock();
    FFT_iter(array, out, size);
    t2 = clock();
    printf("Iter : T2 - T2 = %.7f\n", (double) (t2 - t1) / CLOCKS_PER_SEC);
    for (int i = 0; i < size; ++i) {
        printf("%.3f %s %.3fi\n", creal(out[i]), cimag(out[i]) > 0 ? "+" : "-",
               fabs(cimag(out[i])));
    }


    free(array);
    free(out);
}
