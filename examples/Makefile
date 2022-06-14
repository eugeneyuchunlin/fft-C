
include ../common.mk

CFLAGS += -I../include

all: example fast_mul

%.o: %.c
	$(VECHO) "  CC\t$@\n"
	$(Q)$(CC)  -c $< -o $@ $(CFLAGS)

example: example.o ../src/fft.o
	$(VECHO) "  LD\t$@\n"
	$(Q)$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)

fast_mul: fast_mult.o ../src/fft.o
	$(VECHO) "  LD\t$@\n"
	$(Q)$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)

