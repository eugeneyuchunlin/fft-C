#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "fft.h"

void test(int test_size);

int main(){
    int test_cases[] = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
    30, 900, 1024, 1600, 2160, 3840, 4320, 7680, 7*11*10*13*31};
    int size = sizeof(test_cases)/sizeof(int);
    for(int i = 0; i < size; ++i){
        fprintf(stdout, "Running test case % 3d for size %- 10d...", i+1, test_cases[i]);
        fflush(stdout);
        test(test_cases[i]);
        fprintf(stdout, "Pass!\n");
    }
    return 0;
}

void test(int test_size){
    complex double *array = (complex double *)malloc(sizeof(complex double)*test_size);
    complex double *out = (complex double *)malloc(sizeof(complex double)*test_size);
    complex double *out_legacy = (complex double *)malloc(sizeof(complex double)*test_size);

   for (int i = 0; i < test_size; ++i) {
        array[i] = i;
    }
    FFT(array, out, test_size);
    fprintf(stdout, "[FFT]");
    fflush(stdout);
    legacy_FFT(array, out_legacy, test_size);
    fprintf(stdout, "[Legacy]");
    fflush(stdout); 
    for(int i = 0; i < test_size; ++i){
        assert(fabs(creal(out[i]) - creal(out_legacy[i])) < 0.001);
        assert(fabs(cimag(out[i]) - cimag(out_legacy[i])) < 0.001);
    } 

    free(array);
    free(out);
    free(out_legacy);
}
