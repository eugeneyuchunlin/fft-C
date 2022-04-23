CC := gcc
CFLAGS := -O3 

all : example verify

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ -lm

example : fft.o example.o
	$(CC) $(CFLAGS) fft.o example.o -o example -lm

verify: fft.o verify.o
	$(CC) $(CFLAGS) fft.o verify.o -o verify -lm

clean:
	-rm example
	-rm *.o

