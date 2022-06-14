include ../common.mk

CFLAGS += -I../include
OBJS = fft.o


all: $(OBJS)

%.o : %.c
	$(VECHO) "  CC\t$@\n"
	$(Q)$(CC) -c $< -o $@ $(CFLAGS)

clean:
	-rm *.o
