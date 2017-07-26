#include <stdlib.h>
#include <stdio.h>
#include "par.h"

FILE *open_par( char *filename ) {

    // Simply open the file (for reading) and do basic error checking
    FILE *f = fopen( filename, "r" );
    if (f == NULL) {
        fprintf( stderr, "error: Couldn't open par file '%s' for reading...",
                 filename );
        exit(EXIT_FAILURE);
    }

    // Return the file handle
    return f;
}
