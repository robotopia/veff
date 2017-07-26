INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

CC = gcc
LDFLAGS = -lm

TARGETS = veff parabfit bin2asc

all: $(TARGETS)

veff: veff.c vec.o par.o
	gcc -Wall -Wextra -fsanitize=address -o $@ $^ -I$(CSPICE_DIR)/include -L$(CSPICE_DIR)/lib -lasan -lcspice -lm

test: 0437.par
	./veff -p $< -e 56559.878 -v -s /usr/local/share/jpl/de430.bsp

0437.par:
	psrcat -e2 J0437-4715 > $@

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
	$(RM) veff

