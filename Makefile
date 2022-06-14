VALGRIND := valgrind

include common.mk

all: 
	@$(MAKE) -C src
	@$(MAKE) -C examples
	@$(MAKE) -C tests

valgrind: example
	$(VALGRIND) --leak-check=full ./example	

clean:
	$(MAKE) -C src clean

.PHONY: clean
