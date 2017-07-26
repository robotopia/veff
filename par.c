#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "par.h"

FILE *open_par( char *filename )
{
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

int get_par_double( FILE *f, char *param, double *val, double *err )
/* Get a specific value/error pair from the par file.
 * INPUTS:
 *   f     = par file handle
 *   param = name of desired parameter
 * OUTPUTS:
 *   val   = where the value will be put
 *   err   = where the error will be put
 */
{

    int MAX_STR_LEN = 4096;
    int nscan;              // The number of items scanned by sscanf
    char line[MAX_STR_LEN]; // Each line of the file gets put into here
    char word[MAX_STR_LEN]; // The parameter names get put into here

    // Read in the file line by line
    rewind(f);
    while (fgets( line, MAX_STR_LEN, f) != NULL)
    {
        nscan = sscanf( line, "%s %lf %lf", word, val, err );
        if (strcmp(word, param) == 0) // = Found a match
            return nscan-1; // Return number of (val,err) found
    }

    *val = NAN;
    *err = NAN;
    return -1; // -1 = Didn't find param in the par file
}
