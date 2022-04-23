CC := gcc
CFLAGS := -O3

all : example verify

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ -lm

example : example.o fft.o
	$(CC) $(CFLAGS) example.o fft.o -o example -lm

verify: fft.o verify.o 
	$(CC) $(CFLAGS)  fft.o verify.o -o verify -lm

clean:
	-rm example
	-rm verify
	-rm *.o

.PHONY: clean

