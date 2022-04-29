CC := gcc
CFLAGS := -Ofast -Wall


all : example verify

VALGRIND := valgrind

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

valgrind: example
	$(VALGRIND) --leak-check=full ./example	

clean:
	-rm example
	-rm verify
	-rm *.o

.PHONY: clean
