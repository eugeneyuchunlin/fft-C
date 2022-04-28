CC := gcc
CFLAGS := -O3 -Wall


all : example verify

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	CFLAGS += -lm
endif

%.o : %.c
	$(CC)  -c $< -o $@ $(CFLAGS)

example : example.o fft.o
	$(CC) example.o fft.o -o example $(CFLAGS)

verify: fft.o verify.o 
	$(CC) fft.o verify.o -o verify $(CFLAGS) 

clean:
	-rm example
	-rm verify
	-rm *.o

.PHONY: clean
