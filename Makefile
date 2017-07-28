INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

CC = gcc
LDFLAGS = -lm

TARGETS = veff

all: $(TARGETS) croptest

veff: veff.c vec.o par.o
	gcc -Wall -Wextra -fsanitize=address -o $@ $^ -I$(CSPICE_DIR)/include -L$(CSPICE_DIR)/lib -lasan -lcspice -lm

test:
	$(MAKE) -C test_example

croptest: croptest.c ss.o
	gcc -o $@ $^ -lm

parabfit: pf.c
	gcc -o $@ $< -lm -lgsl -lgslcblas

bin2asc: bin2asc.c
	gcc -o $@ $<

#vec.o: vec.c
#	gcc -c -o $@ $<

install:
	cp veff $(INSTALL_DIR)
	cp parabfit $(INSTALL_DIR)

clean:
	$(RM) veff *.o

