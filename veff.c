#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "SpiceUsr.h"
#include "par.h"
#include "vec.h"

#define PI          (4*atan(1.0))
#define DEG2RAD(x)  ((x)*PI/180.0)
#define RAD2DEG(x)  ((x)*180.0/PI)

void usage()
{
    printf("usage: veff [OPTIONS]\n");
    printf("  OPTIONS:\n");
    printf("    -e MJD     MJD epoch to calculate [required]\n");
    printf("    -h         Display this help and exit\n");
    printf("    -p PULSAR  Select pulsar PULSAR from psrcat catalogue [required]\n");
    printf("    -s SPK     SPK = NASA planetary ephemeris file [required] (e.g.\n");
    printf("               http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)\n");
    printf("    -v         Turn verbose on (=show working)\n");
    printf("\n");
}

int main( int argc, char *argv[] )
{

    char   *par       = NULL;  // The name of the pulsar ephemeris par file
    char   *psr       = NULL;  // The name of the pulsar
    char   *ephemfile = NULL;  // The name of the NASA planetary ephemeris file
    double  epoch     = 0.0;   // The epoch in question
    int     verbose   = 0;     // 0 = verbose output off, 1 = verbose output on

    // Parse command line options
    if (argc <= 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int co;
    while ((co = getopt(argc, argv, "e:hp:s:v")) != -1)
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
                par = strdup(optarg);
                break;
            case 's':
                ephemfile = strdup(optarg);
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
    if (par == NULL)
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

    if (ephemfile == NULL)
    {
        fprintf(stderr, "error: option -s is required\n");
        usage();
        exit(EXIT_FAILURE);
    }

    if (verbose)
    {
        printf( "Inputs:\n" );
        printf( "  Ephemeris = %s\n",  ephemfile );
        printf( "  Epoch     = %lf\n", epoch );
        printf( "  Par file  = %s\n",  par );
    }

    // Open ephemeris file
    double jd = epoch + 2400000.5;
    SpiceDouble et = (jd - j2000_c() ) * spd_c();
    SpiceDouble state[6];
    SpiceDouble lt;
    furnsh_c( ephemfile );
    spkezr_c( "earth", et, "j2000", "NONE", "solar system barycenter", state, &lt );
    unload_c( ephemfile );

    // Normalise the velocity vector
    vec v_earth, vn_earth;
    v_earth.x = state[3];
    v_earth.y = state[4];
    v_earth.z = state[5];
    normalise( &v_earth, &vn_earth );

    // Print out position and velocity of Earth
    if (verbose)
    {
        double v_earth_mag = magnitude( &v_earth );

        printf("\nEarth pos (AU):\n  [%lf, %lf, %lf]\n", KM2AU(state[0]), KM2AU(state[1]), KM2AU(state[2]));
        printf("Earth vel (km/s):\n  [%lf, %lf, %lf]\n", v_earth.x, v_earth.y, v_earth.z );
        printf("  Total: %lf km/s\n", v_earth_mag );
        printf("Earth vel (normalised):\n  [%lf, %lf, %lf]\n", vn_earth.x, vn_earth.y, vn_earth.z );
    }

    // Collect the needed values from par file
    FILE *fpar = open_par( par );

    double sini, ecc, pb, a1, om, kom, dist, pmra, pmdec, rajd, decjd;
    double sini_err, ecc_err, pb_err, a1_err, om_err, kom_err, dist_err,
           pmra_err, pmdec_err, rajd_err, decjd_err;

    int a = get_par_double( fpar, "SINI",   &sini,  &sini_err  );
    printf("get_par_double returned %d\n", a);
    get_par_double( fpar, "ECC",    &ecc,   &ecc_err   );
    get_par_double( fpar, "PB",     &pb,    &pb_err    );
    get_par_double( fpar, "A1",     &a1,    &a1_err    );
    get_par_double( fpar, "OM",     &om,    &om_err    );
    get_par_double( fpar, "KOM",    &kom,   &kom_err   );
    get_par_double( fpar, "DIST_A", &dist,  &dist_err  );
    get_par_double( fpar, "PMRA",   &pmra,  &pmra_err  );
    get_par_double( fpar, "PMDEC",  &pmdec, &pmdec_err );
    get_par_double( fpar, "RAJD",   &rajd,  &rajd_err  );
    get_par_double( fpar, "DECJD",  &decjd, &decjd_err );

    if (verbose) {
        printf("\nValues from %s:\n", par);
        printf("  sini  = %.12f\n", sini);
        printf("  ecc   = %.12f\n", ecc);
        printf("  pb    = %.12f days\n", pb);
        printf("  a1    = %.12f lt sec\n", a1);
        printf("  om    = %.12f deg\n", om);
        printf("  kom   = %.10f deg\n", kom);
        printf("  dist  = %.12f kpc\n", dist);
        printf("  pmra  = %.10f mas/yr\n", pmra);
        printf("  pmdec = %.10f mas/yr\n", pmdec);
        printf("  rajd  = %.10f deg\n", rajd);
        printf("  decjd = %.10f deg\n", decjd);
    }

    // Derive other values
    double cosi  = sqrt(1-sini*sini);
    double komr  = DEG2RAD(kom);
    double rajr  = DEG2RAD(rajd);
    double decjr = DEG2RAD(decjd);
    double tvra  = PM2TV(pmra,dist);
    double tvdec = PM2TV(pmdec,dist);

    // Other constants
    double c = 2.99792458e8;

    if (verbose) {
        printf("\nOther needed/derived values:\n");
        printf("  c     = %.4f m/s\n", c);
        printf("  cosi  = %.12f\n", cosi);
        printf("  komr  = %.12f rad\n", komr);
        printf("  rajr  = %.12f rad\n", rajr);
        printf("  decjr = %.12f rad\n", decjr);
        printf("  tvra  = %.12f km/s\n", tvra);
        printf("  tvdec = %.12f km/s\n", tvdec);
    }

    // Set "up" and target reference frame (unit) vectors.
    // "up" = z-axis in J2000 frame
    // I = in direction of positive RA,  in the plane of the sky
    // J = in direction of positive Dec, in the plane of the sky
    // K = the line of sight towards the pulsar (orthogonal to X and Y)
    vec up;
    vec I, J, K;

    up.x = 0.0;
    up.y = 0.0;
    up.z = 1.0;

    // Calculate K
    K.x = cos(decjr) * cos(rajr);
    K.y = cos(decjr) * sin(rajr);
    K.z = sin(decjr);

    cross_norm( &up, &K, &I ); // Calculate I
    cross( &K, &I, &J );       // Calculate J

    // Get the Earth's velocity projected onto I and J (i.e. the plane of the sky)
    double vI = proj_length( &v_earth, &I );
    double vJ = proj_length( &v_earth, &J );

    // Rotate to line up with the pulsar orbit's line of nodes
    double skomr = sin(komr);
    double ckomr = cos(komr);
    double vx = ckomr*vI + skomr*vJ;
    double vy = skomr*vI - ckomr*vJ;

    printf( "\nVelocity of the Earth in 'x,y' coordinates:\n" );
    printf( "  [%lf, %lf]\n", vx, vy );

    // Free up memory
    free( par );
    free( ephemfile );

    fclose( fpar );

    return 0;
}
