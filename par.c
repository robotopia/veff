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
/* Read the whole contents of a pulsar ephemeris (par) file, and populate
 * a pardata struct with it.
 */
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

    // Firstly, initialise all pointers to NULL
    pd->psrj    = NULL;
    pd->name    = NULL;
    pd->raj     = NULL;
    pd->decj    = NULL;
    pd->survey  = NULL;
    pd->binary  = NULL;
    pd->assoc   = NULL;
    pd->bincomp = NULL;
    pd->type    = NULL;
    pd->osurvey = NULL;
    pd->units   = NULL;

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
        // Do the "error" column first:
        if (nscan == 3) // i.e. if there is an "error" value
        {
            if (strcmp(param, "RAJ"    ) == 0) pd->decj_err   = atof(err);
            else if (strcmp(param, "DECJ"   ) == 0) pd->decj_err   = atof(err);
            else if (strcmp(param, "DM"     ) == 0) pd->dm_err     = atof(err);
            else if (strcmp(param, "F0"     ) == 0) pd->f0_err     = atof(err);
            else if (strcmp(param, "F1"     ) == 0) pd->f1_err     = atof(err);
            else if (strcmp(param, "P0"     ) == 0) pd->p0_err     = atof(err);
            else if (strcmp(param, "P1"     ) == 0) pd->p1_err     = atof(err);
            else if (strcmp(param, "PMRA"   ) == 0) pd->pmra_err   = atof(err);
            else if (strcmp(param, "PMDEC"  ) == 0) pd->pmdec_err  = atof(err);
            else if (strcmp(param, "S150"   ) == 0) pd->s150_err   = atof(err);
            else if (strcmp(param, "S400"   ) == 0) pd->s400_err   = atof(err);
            else if (strcmp(param, "S600"   ) == 0) pd->s600_err   = atof(err);
            else if (strcmp(param, "S1400"  ) == 0) pd->s1400_err  = atof(err);
            else if (strcmp(param, "S5000"  ) == 0) pd->s5000_err  = atof(err);
            else if (strcmp(param, "W50"    ) == 0) pd->w50_err    = atof(err);
            else if (strcmp(param, "W10"    ) == 0) pd->w10_err    = atof(err);
            else if (strcmp(param, "PB"     ) == 0) pd->pb_err     = atof(err);
            else if (strcmp(param, "ECC"    ) == 0) pd->ecc_err    = atof(err);
            else if (strcmp(param, "A1"     ) == 0) pd->a1_err     = atof(err);
            else if (strcmp(param, "T0"     ) == 0) pd->t0_err     = atof(err);
            else if (strcmp(param, "OM"     ) == 0) pd->om_err     = atof(err);
            else if (strcmp(param, "OMDOT"  ) == 0) pd->omdot_err  = atof(err);
            else if (strcmp(param, "BIGOM"  ) == 0) pd->bigom_err  = atof(err);
            else if (strcmp(param, "DIST_A" ) == 0) pd->dist_a_err = atof(err);
            else if (strcmp(param, "PBDOT"  ) == 0) pd->pbdot_err  = atof(err);
            else if (strcmp(param, "SPINDX" ) == 0) pd->spindx_err = atof(err);
            else if (strcmp(param, "RM"     ) == 0) pd->rm_err     = atof(err);
            else if (strcmp(param, "S700"   ) == 0) pd->s700_err   = atof(err);
            else if (strcmp(param, "S3000"  ) == 0) pd->s3000_err  = atof(err);
            else if (strcmp(param, "PX"     ) == 0) pd->px_err     = atof(err);
            else if (strcmp(param, "S1600"  ) == 0) pd->s1600_err  = atof(err);
            else if (strcmp(param, "FB0"    ) == 0) pd->fb0_err    = atof(err);
            else if (strcmp(param, "FB1"    ) == 0) pd->fb1_err    = atof(err);
            else if (strcmp(param, "PMTOT"  ) == 0) pd->pmtot_err  = atof(err);
            else if (strcmp(param, "PML"    ) == 0) pd->pml_err    = atof(err);
            else if (strcmp(param, "PMB"    ) == 0) pd->pmb_err    = atof(err);
            else if (strcmp(param, "KOM"    ) == 0) pd->kom_err    = atof(err);
            else if (strcmp(param, "KIN"    ) == 0) pd->kin_err    = atof(err);
            else if (strcmp(param, "M2"     ) == 0) pd->m2_err     = atof(err);
            else if (strcmp(param, "SINI"   ) == 0) pd->sini_err   = atof(err);
        }

        // Now get the ordinary values

        // Here are the "char *"s:
        if (strcmp(param, "PSRJ"   ) == 0)  pd->psrj    = strdup(val);
        else if (strcmp(param, "NAME"   ) == 0)  pd->name    = strdup(val);
        else if (strcmp(param, "RAJ"    ) == 0)  pd->raj     = strdup(val);
        else if (strcmp(param, "DECJ"   ) == 0)  pd->decj    = strdup(val);
        else if (strcmp(param, "SURVEY" ) == 0)  pd->survey  = strdup(val);
        else if (strcmp(param, "BINARY" ) == 0)  pd->binary  = strdup(val);
        else if (strcmp(param, "ASSOC"  ) == 0)  pd->assoc   = strdup(val);
        else if (strcmp(param, "BINCOMP") == 0)  pd->bincomp = strdup(val);
        else if (strcmp(param, "TYPE"   ) == 0)  pd->type    = strdup(val);
        else if (strcmp(param, "OSURVEY") == 0)  pd->osurvey = strdup(val);
        else if (strcmp(param, "UNITS"  ) == 0)  pd->units   = strdup(val);
        // Here are the doubles:
        else if (strcmp(param, "ELONG"   ) == 0)  pd->elong    = atof(val);
        else if (strcmp(param, "ELAT"    ) == 0)  pd->elat     = atof(val);
        else if (strcmp(param, "DM"      ) == 0)  pd->dm       = atof(val);
        else if (strcmp(param, "PEPOCH"  ) == 0)  pd->pepoch   = atof(val);
        else if (strcmp(param, "F0"      ) == 0)  pd->f0       = atof(val);
        else if (strcmp(param, "F1"      ) == 0)  pd->f1       = atof(val);
        else if (strcmp(param, "P0"      ) == 0)  pd->p0       = atof(val);
        else if (strcmp(param, "P1"      ) == 0)  pd->p1       = atof(val);
        else if (strcmp(param, "DIST_DM" ) == 0)  pd->dist_dm  = atof(val);
        else if (strcmp(param, "DIST_DM1") == 0)  pd->dist_dm1 = atof(val);
        else if (strcmp(param, "PMRA"    ) == 0)  pd->pmra     = atof(val);
        else if (strcmp(param, "PMDEC"   ) == 0)  pd->pmdec    = atof(val);
        else if (strcmp(param, "S150"    ) == 0)  pd->s150     = atof(val);
        else if (strcmp(param, "S400"    ) == 0)  pd->s400     = atof(val);
        else if (strcmp(param, "S600"    ) == 0)  pd->s600     = atof(val);
        else if (strcmp(param, "S1400"   ) == 0)  pd->s1400    = atof(val);
        else if (strcmp(param, "S5000"   ) == 0)  pd->s5000    = atof(val);
        else if (strcmp(param, "W50"     ) == 0)  pd->w50      = atof(val);
        else if (strcmp(param, "W10"     ) == 0)  pd->w10      = atof(val);
        else if (strcmp(param, "POSEPOCH") == 0)  pd->posepoch = atof(val);
        else if (strcmp(param, "DMEPOCH" ) == 0)  pd->dmepoch  = atof(val);
        else if (strcmp(param, "PB"      ) == 0)  pd->pb       = atof(val);
        else if (strcmp(param, "ECC"     ) == 0)  pd->ecc      = atof(val);
        else if (strcmp(param, "A1"      ) == 0)  pd->a1       = atof(val);
        else if (strcmp(param, "T0"      ) == 0)  pd->t0       = atof(val);
        else if (strcmp(param, "OM"      ) == 0)  pd->om       = atof(val);
        else if (strcmp(param, "OMDOT"   ) == 0)  pd->omdot    = atof(val);
        else if (strcmp(param, "BIGOM"   ) == 0)  pd->bigom    = atof(err);
        else if (strcmp(param, "DIST_A"  ) == 0)  pd->dist_a   = atof(val);
        else if (strcmp(param, "AGE"     ) == 0)  pd->age      = atof(val);
        else if (strcmp(param, "PBDOT"   ) == 0)  pd->pbdot    = atof(val);
        else if (strcmp(param, "SPINDX"  ) == 0)  pd->spindx   = atof(val);
        else if (strcmp(param, "RM"      ) == 0)  pd->rm       = atof(val);
        else if (strcmp(param, "S700"    ) == 0)  pd->s700     = atof(val);
        else if (strcmp(param, "S3000"   ) == 0)  pd->s3000    = atof(val);
        else if (strcmp(param, "PX"      ) == 0)  pd->px       = atof(val);
        else if (strcmp(param, "S1600"   ) == 0)  pd->s1600    = atof(val);
        else if (strcmp(param, "FB0"     ) == 0)  pd->fb0      = atof(val);
        else if (strcmp(param, "FB1"     ) == 0)  pd->fb1      = atof(val);
        else if (strcmp(param, "DIST"    ) == 0)  pd->dist     = atof(val);
        else if (strcmp(param, "DIST1"   ) == 0)  pd->dist1    = atof(val);
        else if (strcmp(param, "PMERR_PA") == 0)  pd->pmerr_pa = atof(val);
        else if (strcmp(param, "P1_I"    ) == 0)  pd->p1_i     = atof(val);
        else if (strcmp(param, "AGE_I"   ) == 0)  pd->age_i    = atof(val);
        else if (strcmp(param, "BSURF_I" ) == 0)  pd->bsurf_i  = atof(val);
        else if (strcmp(param, "EDOT_I"  ) == 0)  pd->edot_i   = atof(val);
        else if (strcmp(param, "EDOTD2"  ) == 0)  pd->edotd2   = atof(val);
        else if (strcmp(param, "R_LUM"   ) == 0)  pd->r_lum    = atof(val);
        else if (strcmp(param, "R_LUM14" ) == 0)  pd->r_lum14  = atof(val);
        else if (strcmp(param, "PMTOT"   ) == 0)  pd->pmtot    = atof(val);
        else if (strcmp(param, "VTRANS"  ) == 0)  pd->vtrans   = atof(val);
        else if (strcmp(param, "BSURF"   ) == 0)  pd->bsurf    = atof(val);
        else if (strcmp(param, "B_LC"    ) == 0)  pd->b_lc     = atof(val);
        else if (strcmp(param, "SI414"   ) == 0)  pd->si414    = atof(val);
        else if (strcmp(param, "EDOT"    ) == 0)  pd->edot     = atof(val);
        else if (strcmp(param, "RAJD"    ) == 0)  pd->rajd     = atof(val);
        else if (strcmp(param, "DECJD"   ) == 0)  pd->decjd    = atof(val);
        else if (strcmp(param, "MASSFN"  ) == 0)  pd->massfn   = atof(val);
        else if (strcmp(param, "MINMASS" ) == 0)  pd->minmass  = atof(val);
        else if (strcmp(param, "MEDMASS" ) == 0)  pd->medmass  = atof(val);
        else if (strcmp(param, "DMSINB"  ) == 0)  pd->dmsinb   = atof(val);
        else if (strcmp(param, "GL"      ) == 0)  pd->gl       = atof(val);
        else if (strcmp(param, "GB"      ) == 0)  pd->gb       = atof(val);
        else if (strcmp(param, "XX"      ) == 0)  pd->xx       = atof(val);
        else if (strcmp(param, "YY"      ) == 0)  pd->yy       = atof(val);
        else if (strcmp(param, "ZZ"      ) == 0)  pd->zz       = atof(val);
        else if (strcmp(param, "PML"     ) == 0)  pd->pml      = atof(val);
        else if (strcmp(param, "PMB"     ) == 0)  pd->pmb      = atof(val);
        else if (strcmp(param, "KOM"     ) == 0)  pd->kom      = atof(val);
        else if (strcmp(param, "KIN"     ) == 0)  pd->kin      = atof(val);
        else if (strcmp(param, "M2"      ) == 0)  pd->m2       = atof(val);
        else if (strcmp(param, "UPRMASS" ) == 0)  pd->uprmass  = atof(val);
        else if (strcmp(param, "MINOMDOT") == 0)  pd->minomdot = atof(val);
        else if (strcmp(param, "SINI"    ) == 0)  pd->sini     = atof(val);
        // Here are the ints:
        else if (strcmp(param, "DATE"    ) == 0)  pd->date     = atoi(val);
    }

    // Close the par file
    fclose( fpar );
}


void free_par( struct pardata *pd )
/* Free memory allocated in pardata struct (but not pardata itself) */
{
    free( pd->psrj );
    free( pd->name );
    free( pd->raj );
    free( pd->decj );
    free( pd->survey );
    free( pd->binary );
    free( pd->assoc );
    free( pd->bincomp );
    free( pd->type );
    free( pd->osurvey );
    free( pd->units );
}

