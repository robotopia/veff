INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

TARGETS = veff parabfit bin2asc

all: $(TARGETS)

veff: veff.c
	gcc -Wall -Wextra -fsanitize=address -o $@ $< -I$(CSPICE_DIR)/include -L$(CSPICE_DIR)/lib -lasan -lpsrcat -lcspice -lm

test: veff
	./$< -p J0437-4715 -e 56559.878 -v -s /usr/local/share/jpl/de430.bsp

parabfit: pf.c
	gcc -o $@ $< -lm -lgsl -lgslcblas

bin2asc: bin2asc.c
	gcc -o $@ $<

install:
	cp veff $(INSTALL_DIR)
	cp parabfit $(INSTALL_DIR)

clean:
	$(RM) veff

