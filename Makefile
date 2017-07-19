INSTALL_DIR = ~/bin
TARGETS = parabfit bin2asc
LIBS = -lm -lgsl -lgslcblas

all: $(TARGETS)

parabfit: pf.c
	gcc -o $@ $< $(LIBS)

bin2asc: bin2asc.c
	gcc -o $@ $<

install:
	mv parabfit $(INSTALL_DIR)

clean:
	-rm $(TARGETS)
