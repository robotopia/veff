TARGETS = parabfit bin2asc
LIBS = -lm -lgsl -lgslcblas

all: $(TARGETS)

parabfit: pf.c
	gcc -o $@ $< $(LIBS)

bin2asc: bin2asc.c
	gcc -o $@ $<

clean:
	-rm $(TARGETS)
