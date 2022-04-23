#include <complex.h>  // for complex number
#include <math.h>     // for cos, sin
#include <stdio.h>    // for srandard i/o
#include <stdlib.h>   // for memory mangement
#include <assert.h>
#include <time.h>
#include "fft.h"



int main()
{
    const int size =  7*11*17*30*30;
    
//     fft_task_queue_t queue;
//     INIT_LIST_HEAD(&queue.list);
// 
//     fft_task_t t1;
//     t1.radix = 2;
//     enqueue(&t1.node, &queue.list);
// 
//     struct list_head *n;
//     dequeue(&n, &queue.list);
//     dequeue(&n, &queue.list); 
//     fft_task_t *task1 = container_of(n, fft_task_t, node);
//     printf("Task %d\n", task1->radix);
// 
    
    clock_t t1, t2; 
    // complex double array[4] = {1, 3, 5, 7};
    complex double *array = (complex double *)malloc(sizeof(complex double)*size);
    complex double *out = (complex double *)malloc(sizeof(complex double)*size);
    complex double *out_legacy = (complex double *)malloc(sizeof(complex double)*size);
    
    array[0] = 1;
    for (int i = 0; i < size; ++i) {
        array[i] = i;
    }
    t1 = clock();
    FFT_iter(array, out, size);
    t2 = clock();
    printf("Iter : T2 - T2 = %.7f\n", (double)(t2 - t1) / CLOCKS_PER_SEC);

    t1 = clock();
    FFT(array, out_legacy, size);
    t2 = clock();
    printf("Recursive : T2 - T2 = %.7f\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
    for(int i = 0; i < size; ++i){
        assert(fabs(creal(out[i]) - creal(out_legacy[i])) < 0.001);
        assert(fabs(cimag(out[i]) - cimag(out_legacy[i])) < 0.001);
    } 
    // for (int i = 0; i < size; ++i) {
    //     printf("%.3f %s %.3fi\n", creal(out[i]), cimag(out[i]) > 0 ? "+" : "-",
    //             fabs(cimag(out[i])));
    // }
    free(array);
    free(out);
    free(out_legacy);
}
