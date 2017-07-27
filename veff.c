#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "SpiceUsr.h"
#include "par.h"
#include "vec.h"
#include "ss.h"

#define PI          (4*atan(1.0))
#define DEG2RAD(x)  ((x)*PI/180.0)
#define RAD2DEG(x)  ((x)*180.0/PI)

void usage()
{
    printf("usage: veff [OPTIONS]\n");
    printf("  OPTIONS:\n");
    printf("    -e MJD   MJD epoch to calculate [required]\n");
    printf("    -f FREQ  Observing frequency (MHz) [required]\n");
    printf("    -g       Turn on Geogebra-style output\n");
    printf("    -h       Display this help and exit\n");
    printf("    -p PAR   Supply (PSRCAT-style) ephemeris file [required]\n");
    printf("    -s SPK   NASA planetary ephemeris file [required] (e.g.");
    printf(" from\n             http://naif.jpl.nasa.gov/pub/naif/");
    printf("generic_kernels/spk/planets/de430.bsp)\n");
    printf("    -v      Turn verbose on (=show working)\n");
    printf("\n");
}

int main( int argc, char *argv[] )
{

    char   *par       = NULL;  // The name of the pulsar ephemeris par file
    char   *ephemfile = NULL;  // The name of the NASA planetary ephemeris file
    double  epoch     = 0.0;   // The epoch in question
    double  freq      = 0.0;   // Observing frequency
    int     verbose   = 0;     // 0 = verbose output off, 1 = verbose output on
    int     geogebra  = 0;     // 0 = verbose output off, 1 = verbose output on

    struct pardata pd;         // A container to hold the par file information
    struct sec_spect ss;       // The secondary spectrum information

    /* Just dummy data for now */
    ss.dx = 0.0;
    ss.dy = 0.0;
    ss.min_x = 0.0;
    ss.max_x = 0.0;

    // Parse command line options
    if (argc <= 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int co;
    while ((co = getopt(argc, argv, "e:f:ghp:s:v")) != -1)
    {
        switch (co)
        {
            case 'e':
                epoch = atof(optarg);
                break;
            case 'f':
                freq = atof(optarg);
                break;
            case 'g':
                geogebra = 1;
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

    if (freq <= 0.0)
    {
        fprintf(stderr, "error: invalid frequency %f\n", freq);
        fprintf(stderr, "Run '%s -h' for options\n", argv[0]);
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
    spkezr_c( "earth", et, "j2000", "NONE", "solar system barycenter",
              state, &lt );
    unload_c( ephemfile );

    // Normalise the position and velocity vectors
    vec r_earth, rn_earth; // position vector (kms)
    vec v_earth, vn_earth; // velocity vector (km/s)

    r_earth.x = state[0];   v_earth.x = state[3];
    r_earth.y = state[1];   v_earth.y = state[4];
    r_earth.z = state[2];   v_earth.z = state[5];

    normalise( &v_earth, &vn_earth );
    normalise( &r_earth, &rn_earth );

    // Print out position and velocity of Earth
    if (verbose)
    {
        printf("Earth vel (km/s):\n  [%lf, %lf, %lf]\n",
                v_earth.x, v_earth.y, v_earth.z );
    }

    // Collect the needed values from par file
    read_par( par, &pd );

    if (verbose) {
        printf("\nValues from %s:\n", par);
        printf("  sini  = %.12f\n", pd.sini);
        printf("  ecc   = %.12f\n", pd.ecc);
        printf("  pb    = %.12f days\n", pd.pb);
        printf("  a1    = %.12f lt sec\n", pd.a1);
        printf("  om    = %.12f deg\n", pd.om);
        printf("  kom   = %.10f deg\n", pd.kom);
        printf("  dist  = %.12f kpc\n", pd.dist);
        printf("  pmra  = %.10f mas/yr\n", pd.pmra);
        printf("  pmdec = %.10f mas/yr\n", pd.pmdec);
        printf("  rajd  = %.10f deg\n", pd.rajd);
        printf("  decjd = %.10f deg\n", pd.decjd);
    }

    // Derive other values
    double cosi  = sqrt(1 - pd.sini*pd.sini);
    double komr  = DEG2RAD(pd.kom);
    double rajr  = DEG2RAD(pd.rajd);
    double decjr = DEG2RAD(pd.decjd);
    double tvra  = PM2TV(pd.pmra, pd.dist);
    double tvdec = PM2TV(pd.pmdec, pd.dist);

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

    if (geogebra)
    {
        printf("\nGeogebra values:\n");
        printf("  B_a          = %e\n", pd.a1);
        printf("  D            = %e\n", pd.dist);
        printf("  Dec_{Pulsar} = %e\n", decjr);
        printf("  E_{now}      = (%e, %e, %e)\n",
                  rn_earth.x, rn_earth.y, rn_earth.z );
        printf("  P_b          = %e\n", pd.pb);
        printf("  RA_{Pulsar}  = %e\n", rajr);
        printf("  T_0          = %e\n", pd.t0);
        printf("  ecc          = %e\n", pd.ecc);
        printf("  epoch        = %e\n", epoch);
        printf("  freq         = %e\n", freq);
        printf("  img_{dx}     = %e\n", ss.dx);
        printf("  img_{ratio}  = %e\n", ss.dy/ss.dx);
        printf("  img_w        = %.0lf\n", (ss.max_x-ss.min_x)/ss.dx);
        printf("  sini         = %e\n", pd.sini);
        printf("  v_{E,x}      = %e\n", v_earth.x);
        printf("  v_{E,y}      = %e\n", v_earth.y);
        printf("  v_{E,z}      = %e\n", v_earth.z);
        printf("  v_{μ⟂,Dec}   = %e\n", tvdec);
        printf("  v_{μ⟂,RA}    = %e\n", tvra);
        printf("  Ω            = %e\n", pd.bigom);
        printf("  ω            = %e\n", pd.om);
    }

    // Free up memory
    free( par );
    free( ephemfile );
    free_par( &pd );

    return 0;
}
