CC := gcc
CFLAGS := -Ofast -pthread -Wall
LDFLAGS := -pthread
VALGRIND := valgrind

include common.mk

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	CFLAGS += -lm
endif


all: 
	$(MAKE) -C src
	$(MAKE) -C examples
	# $(MAKE) -C tests

valgrind: example
	$(VALGRIND) --leak-check=full ./example	

clean:
	$(MAKE) -C src clean

.PHONY: clean
