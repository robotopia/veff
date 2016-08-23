TARGETS = pf
LIBS = -lm -lgsl -lgslcblas

all: $(TARGETS)

%: %.c
	gcc -o pf pf.c -lm -lgsl -lgslcblas

clean:
	-rm $(TARGETS)
