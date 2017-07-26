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

    int nscan;              // The number of items scanned by sscanf
    char line[MAXSTRLEN];   // Each line of the file gets put into here
    char word[MAXSTRLEN];   // The parameter names get put into here

    // Read in the file line by line
    rewind(f);
    while (fgets( line, MAXSTRLEN, f) != NULL)
    {
        nscan = sscanf( line, "%s %lf %lf", word, val, err );
        if (strcmp(word, param) == 0) // = Found a match
            return nscan-1; // Return number of (val,err) found
    }

    *val = NAN;
    *err = NAN;
    return -1; // -1 = Didn't find param in the par file
}


void read_par( char *filename, struct pardata *pd )
{
    // Make sure neither input argument is a NULL pointer
    if (filename == NULL)
    {
        fprintf( stderr, "error: read_par: filename cannot be NULL\n" );
        exit(EXIT_FAILURE);
    }

    if (pd == NULL)
    {
        fprintf( stderr, "error: read_par: pardata cannot be NULL\n" );
        exit(EXIT_FAILURE);
    }

    // Open the par file
    FILE *fpar = open_par( filename );

    // Read through the file line by line and input the values into the
    // corresponding variables in the pardata struct
    int nscan;              // The number of items scanned by sscanf
    char line[MAXSTRLEN];   // Each line of the file gets put into here
    char param[MAXSTRLEN];  // The parameter names get put into here
    char val[MAXSTRLEN];    // The parameter names get put into here
    char err[MAXSTRLEN];    // The parameter names get put into here
    while (fgets( line, MAXSTRLEN, fpar) != NULL)
    {
        // Read next line
        nscan = sscanf( line, "%s %s %s", param, val, err );

        if (nscan <= 1) continue;  // No values in line, so ignore line

        // Check the param and update the value in the pardata struct
        if (strcmp(param, "PSRJ") == 0)  { pd->psrf = strdup(val);  continue; };
        //... lots more to implement here...
    }

    // Close the par file
    fclose( fpar );
}
