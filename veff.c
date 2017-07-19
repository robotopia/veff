#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "psrcat.h"

#define PI          (4*atan(1.0))
#define DEG2RAD(x)  ((x)*PI/180.0)
#define RAD2DEG(x)  ((x)*180.0/PI)

typedef struct vec_t
{
    double x;
    double y;
    double z;
} vec;

void usage()
{
    printf("usage: veff [OPTIONS]\n");
    printf("  OPTIONS:\n");
    printf("    -e MJD     MJD epoch to calculate [required]\n");
    printf("    -h         Display this help and exit\n");
    printf("    -p PULSAR  Select pulsar PULSAR from psrcat catalogue [required]\n");
    printf("    -v         Turn verbose on (=show working)\n");
    printf("\n");
}

double get_psrcat_value( char *psr, char *param )
{
    double val, err;
    char ref[100];

    int ret;

    ret = callPsrcat_val("public", psr, param, &val, &err, ref);

    switch (ret)
    {
        case 1:
            fprintf(stderr, "error: could not open psrcat catalogue\n");
            exit(EXIT_FAILURE);
            break;
        case 2:
            fprintf(stderr, "error: unrecognised pulsar '%s'\n", psr);
            exit(EXIT_FAILURE);
            break;
        case 3:
            fprintf(stderr, "error: unrecognised parameter\n");
            exit(EXIT_FAILURE);
            break;
    }

    return val;
}

void calc_vearth(vec *vearth)
{
    vearth->x = 0.0;
    vearth->y = 0.0;
    vearth->z = 0.0;
}

int main( int argc, char *argv[] )
{

    char *psr = NULL;
    double epoch = 0.0;
    int verbose = 0;

    // Parse command line options
    if (argc <= 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int co;
    while ((co = getopt(argc, argv, "e:hp:v")) != -1)
    {
        switch (co)
        {
            case 'e':
                epoch = atof(optarg);
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'p':
                psr = strdup(optarg);
                break;
            case 'v':
                verbose = 1;
                break;
            default:
                fprintf(stderr, "error: unknown option '%c'\n", co);
                exit(EXIT_FAILURE);
                break;
        }
    }

    // Error checking on options
    if (psr == NULL)
    {
        fprintf(stderr, "error: option -p is required\n");
        usage();
        exit(EXIT_FAILURE);
    }

    if (epoch <= 0.0)
    {
        fprintf(stderr, "error: invalid epoch %f\n", epoch);
        fprintf(stderr, "Run '%s -h' for options\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Collect the needed values from psrcat
    double sini  = get_psrcat_value( psr, "sini"   );
    double ecc   = get_psrcat_value( psr, "ecc"    );
    double pb    = get_psrcat_value( psr, "pb"     );
    double a1    = get_psrcat_value( psr, "a1"     );
    double om    = get_psrcat_value( psr, "om"     );
    double kom   = get_psrcat_value( psr, "kom"    );
    double dist  = get_psrcat_value( psr, "dist_a" );
    double pmra  = get_psrcat_value( psr, "pmra"   );
    double pmdec = get_psrcat_value( psr, "pmdec"  );

    if (verbose) {
        printf("\nValues from PSRCAT:\n");
        printf("  sini  = %.12f\n", sini);
        printf("  ecc   = %.12f\n", ecc);
        printf("  pb    = %.12f days\n", pb);
        printf("  a1    = %.12f lt sec\n", a1);
        printf("  om    = %.12f deg\n", om);
        printf("  kom   = %.10f deg\n", kom);
        printf("  dist  = %.12f kpc\n", dist);
        printf("  pmra  = %.10f mas/yr\n", pmra);
        printf("  pmdec = %.10f mas/yr\n", pmdec);
    }

    // Derive other values
    double cosi = sqrt(1-sini*sini);
    double komr = DEG2RAD(kom);

    // Other constants
    double c = 2.99792458e8;

    if (verbose) {
        printf("\nOther needed/derived values:\n");
        printf("  c     = %.4f m/s\n", c);
        printf("  cosi  = %.12f\n", cosi);
        printf("  komr  = %.12f rad\n", komr);
    }

    return 0;
}
