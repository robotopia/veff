TARGETS = pf bin2asc
LIBS = -lm -lgsl -lgslcblas

all: $(TARGETS)

%: %.c
	gcc -o $@ $< $(LIBS)

clean:
	-rm $(TARGETS)
