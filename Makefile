INSTALL_DIR = /usr/local/bin
CSPICE_DIR = $(HOME)/src/cspice/cspice

CC = gcc
LDFLAGS = -L$(CSPICE_DIR)/lib
LDLIBS = -lcspice -lm
CFLAGS = -Wall -Wextra -I$(CSPICE_DIR)/include

TARGETS = veff parabfit

all: $(TARGETS)

veff: veff.o vec.o par.o ss.o hough.o

test:
	$(MAKE) -C test_example

parabfit: parabfit.c ss.o hough.o

documentation:
	$(MAKE) -C doc

install:
	cp $(TARGETS) $(INSTALL_DIR)

clean:
	$(RM) $(TARGETS) *.o

