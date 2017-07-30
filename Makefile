INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

CC = gcc
LDFLAGS = -lm
CFLAGS = -Wall -Wextra

TARGETS = veff parabfit

all: $(TARGETS)

veff: veff.c vec.o par.o ss.o hough.o
	gcc -Wall -Wextra -o $@ $^ -I$(CSPICE_DIR)/include -L$(CSPICE_DIR)/lib -lcspice -lm

test:
	$(MAKE) -C test_example

parabfit: parabfit.c ss.o hough.o

documentation:
	$(MAKE) -C doc

install:
	cp $(TARGETS) $(INSTALL_DIR)

clean:
	$(RM) $(TARGETS) *.o

