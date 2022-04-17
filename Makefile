CC := gcc
CFLAGS := -g

all : example

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

example : fft.o example.o
	$(CC) $(CFLAGS) fft.o example.o -o example
