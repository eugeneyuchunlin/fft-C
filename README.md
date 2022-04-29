# fft-C
Fast Fouier Transform in C

This is introduction for the repository and I would like to describe what I have done for implementing fast fourier transform in C language. 
This is the course work in **National Cheng Kung University**. The professor tought me all these magic is [ycshu](https://github.com/ycshu). Mathmetician is megician 

I wrote two versions to implement fast fourier transform algorithm one is using recursive function call the other is using loop to iterate. 

## Usage

It's welcome to use the algorithm but please note the LICENSE.

```c=
#include <complex.h>
#include "fft.h"


int main(){
    ...
    complex double in[200];
    complex double out[200];
    
    // call FFT recursive version which is faster
    FFT(in, out, 200);
    
    // or call FFT iteration version
    FFT_iter(in, out, 200);
    ...
    
    return 0
}
```

## Contribution Guide

It's welcome to contribute the code. Please feel free to have modification. After you done, you can use `verify.c` to check the answer is correct or not. In the pull request I'll review the GitHub action.

## Run the Code

### Compile

If you have `make`
```shell=
make
```

If you don't have `make`
```shell=
gcc -Ofast fft.c example.c -o example
gcc -Ofast fft.c verify.c -o verify
```

## Recursive

TO BE DONE

## Iterative
