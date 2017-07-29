INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

CC = gcc
LDFLAGS = -lm

TARGETS = veff

all: $(TARGETS)

veff: veff.c vec.o par.o ss.o hough.o
	gcc -Wall -Wextra -fsanitize=address -o $@ $^ -I$(CSPICE_DIR)/include -L$(CSPICE_DIR)/lib -lasan -lcspice -lm

test:
	$(MAKE) -C test_example

parabfit: parabfit.c ss.o hough.o
	gcc -o $@ $^ -lm

bin2asc: bin2asc.c
	gcc -o $@ $<

#vec.o: vec.c
#	gcc -c -o $@ $<

install:
	cp veff $(INSTALL_DIR)
	cp parabfit $(INSTALL_DIR)

clean:
	$(RM) veff *.o

