#include <complex.h>  // for complex number
#include <math.h>     // for cos, sin
#include <stdio.h>    // for srandard i/o
#include <stdlib.h>   // for memory mangement
#include <assert.h>

#include "fft.h"

int main()
{
    const int size =  7*17*11*30*30;

    // complex double array[4] = {1, 3, 5, 7};
    complex double *array = (complex double *)malloc(sizeof(complex double)*size);
    complex double *out = (complex double *)malloc(sizeof(complex double)*size);
    // complex double out_legacy[size];
    
    // array[0] = 1;
    for (int i = 1; i < size; ++i) {
        array[i - 1] = i;
    }
    FFT(array, out, size);
    // legacy_FFT(array, out_legacy, size);
    // for(int i = 0; i < size; ++i){
    //     assert(fabs(creal(out[i]) - creal(out_legacy[i])) < 0.001);
    //     assert(fabs(cimag(out[i]) - cimag(out_legacy[i])) < 0.001);
    // } 

    // for (int i = 0; i < size; ++i) {
    //     printf("%.3f %s %.3fi\n", creal(out[i]), cimag(out[i]) > 0 ? "+" : "-",
    //             fabs(cimag(out[i])));
    // }
}
