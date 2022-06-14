CFLAGS += -I../include
OBJS = fft.o

include ../common.mk

all: $(OBJS)

%.o : %.c
	$(VECHO) "  CC\t$@\n"
	$(Q)$(CC) -c $< -o $@ $(CFLAGS)

clean:
	-rm *.o
