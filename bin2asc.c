#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Convert float64 to float32


int main ( int argc, char *argv[] ) {
    char *filename = strdup(argv[1]);
    int rows = atoi(argv[2]);
    int cols = atoi(argv[3]);
    char *outfile = "converted.out";

    if (strcmp(argv[1], "-h") == 0 || argc != 4) {
        fprintf(stderr, "usage: %s [filename] [# of x pixels] [# of y pixels]\n", argv[0]);
        exit(0);
    }

    FILE *fi = fopen(filename, "rb");
    FILE *fo = fopen(outfile,  "w");
    
    double f;
    
    int counter = 0;
    int r, c;

    while (1) {
        fread(&f, sizeof(double), 1, fi);
        if (feof(fi))
            break;
        r = counter % cols;
        c = counter / cols;
        fprintf(fo, "%d %d %lf\n", c, r, f);
        counter++;
    }

    fclose(fo);
    fclose(fi);
}
