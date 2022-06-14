include ../common.mk

CFLAGS += -I../include

all: verify

%.o : %.c
	$(VECHO) "  CC\t$@\n"
	$(Q)$(CC) -c $^ -o $@ $(CFLAGS)

verify: verify.o ../src/fft.o
	$(VECHO) "  LD\t$@\n"
	$(Q)$(CC) $^ -o $@ $(CFALGS) $(LDFLAGS)
