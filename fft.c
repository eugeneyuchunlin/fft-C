#include "fft.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define SQRT3 1.7320508075688771931766041
#define SQRT3_DIV_2 0.8660254037844385965883021

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

void FFT3(complex double in[], complex double out[]){
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
    
    static const complex double sq3_ps_1_d_2 = -0.5 + SQRT3_DIV_2 * I;
    static const complex double sq3_mi_1_d_2 = -0.5 - SQRT3_DIV_2 * I;
    out[0] = in[0] + in[1] + in[2];
    out[1] = in[0] + in[1] * sq3_mi_1_d_2 + in[2] * sq3_ps_1_d_2;
    out[2] = in[0] + in[1] * sq3_ps_1_d_2 + in[2] * sq3_mi_1_d_2;
}

void FFTN(complex double in[], complex double out[], int N){
    complex double mat[N][N];
    double pi_2_DIV_N = M_PI * 2 / N;
    for(int i = 1; i < N; ++i){
        for(int j = 1; j < N; ++j){
            // exp(-I*2pi/N * ij)
            mat[i][j] = cos(pi_2_DIV_N * i * j) - I*sin(pi_2_DIV_N * i * j);
        }
    }

    for(int i = 0; i < N; ++i){
        mat[0][i] = mat[i][0] = 1;
    }
    
    for(int i = 0 ; i < N; ++i){
        complex double sum = 0;
        for(int j = 0; j < N; ++j){
            sum += mat[i][j] * in[j];
        }
        out[i] = sum;
    }
}

void FFT3N(complex double in[], complex double out[], int size){
    if(size == 3){
        FFT3(in, out);
    }else{
        int next_size = size / 3;
        complex double n3[next_size];
        complex double n3_1[next_size];
        complex double n3_2[next_size];

        complex double out3[next_size];
        complex double out3_1[next_size];
        complex double out3_2[next_size];
        
        for(int i = 0; i < next_size; ++i){
            n3[i]   = in[3*i];
            n3_1[i] = in[3*i + 1];
            n3_2[i] = in[3*i + 2];
        }

        FFT3N(n3, out3, next_size);
        FFT3N(n3_1, out3_1, next_size);
        FFT3N(n3_2, out3_2, next_size);

        complex double in_temp[3];
        complex double out_temp[3];

        for(int i = 0; i < next_size; ++i){
            complex double w1 = 
                cos(i * 2 * M_PI / size) - I * sin(i * 2 * M_PI / size);
            complex double w2 = 
                cos(i * 4 * M_PI / size) - I * sin(i * 4 * M_PI / size);

            in_temp[0] = out3[i];
            in_temp[1] = out3_1[i] * w1;
            in_temp[2] = out3_2[i] * w2;

            FFT3(in_temp, out_temp);
            out[i] = out_temp[0];
            out[i + next_size] = out_temp[1];
            out[i + next_size * 2] = out_temp[2];
        }
    }
}


void FFT2N(complex double in[], complex double out[], int size)
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
        for (int i = 0; i < half_size; ++i) {
            even[i] = in[2 * i];
            odd[i] = in[2 * i + 1];
        }

        // FFT seperately
        FFT2N(even, even_result, half_size);
        FFT2N(odd, odd_result, half_size);


        for (int i = 0; i < half_size; ++i) {
            complex double w =
                cos(i * 2 * M_PI / size) - I * sin(i * 2 * M_PI / size);
            out[i] = even_result[i] + odd_result[i] * w;
            out[half_size + i] = even_result[i] + (-1) * odd_result[i] * w;
        }
    }
}

void FFT(complex double in[], complex double out[], int size){
    bool special_flag = false;
    void (*fft_function)(complex double *, complex double *, int) = FFT;
    if(size == 1){
        out[0] = in[0];
    }else if(size == 2){
        FFT2(in, out);
    }else if (size == 3){
        FFT3(in, out);
    }else if (size == 5){
        FFTN(in, out, 5);
    }else{
        int p = 1;
        if(size % 2 == 0){
            p = 2;
        }else if (size % 3 == 0){
            p = 3;
        }else{
            for(int i = 5; i <= size; ++i){
                if(size % i == 0){
                    p = i;
                    break;
                }
            }
        }
        // printf("p = %d\n", p);
        int new_size = size / p;
        complex double **entries_in;
        complex double **entries_out;

        entries_in = (complex double **)malloc(sizeof(complex double*) * p);
        entries_out  = (complex double **)malloc(sizeof(complex double *)*p);

        for(int i = 0; i < p; ++i){
            entries_in[i] = 
                (complex double *)malloc(sizeof(complex double)*new_size);
            entries_out[i] = 
                (complex double *)malloc(sizeof(complex double)*new_size);
        }

      

        // place the data;
        for(int i = 0; i < p; ++i){
            for(int j = 0; j < new_size; ++j){
                entries_in[i][j] = in[j * p + i];
            }
        }

    
        if(p == size) fft_function = FFTN;
        // invoke FFT
        for(int i = 0; i < p; ++i){
            fft_function(entries_in[i], entries_out[i], new_size);
        }

        complex double *in_temp = 
            (complex double *)malloc(sizeof(complex double)*p);
        complex double *out_temp = 
            (complex double *)malloc(sizeof(complex double)*p);
        // collect the answer
         
        for(int i = 0; i < new_size; ++i){
            for(int j = 0; j < p; ++j){
                complex double w = 
                    cos(i * j * 2 * M_PI / size) - 
                    I * sin(i * j * 2 * M_PI / size);
                in_temp[j] = entries_out[j][i] * w;
            }
            fft_function(in_temp, out_temp, p);
            
            for(int j = 0; j < p; ++j){
                out[i + new_size * j] = out_temp[j];
            }
        }
        free(in_temp);
        free(out_temp);
        for(int i = 0; i < p; ++i){
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
