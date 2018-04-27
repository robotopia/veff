#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "SpiceUsr.h"
#include "par.h"
#include "vec.h"
#include "ss.h"
#include "hough.h"

#define PI          (4*atan(1.0))
#define DEG2RAD(x)  ((x)*PI/180.0)
#define RAD2DEG(x)  ((x)*180.0/PI)

void usage()
{
    printf("usage: veff [OPTIONS]\n");
    printf("  OPTIONS:\n");
    printf("    -a ETA   Arc curvature (s^3) [required]\n" );
    printf("    -e MJD   MJD epoch to calculate [required]\n");
    printf("    -f FREQ  Observing frequency (MHz) [required]\n");
    printf("    -h       Display this help and exit\n");
    printf("    -p PAR   Supply (PSRCAT-style) ephemeris file [required]\n");
    printf("    -s SPK   NASA planetary ephemeris file [required] (e.g. "  
                        "from\n" );
    printf("             http://naif.jpl.nasa.gov/pub/naif/generic_kernels/"
                        "spk/planets/de430.bsp)\n");
    printf("    -v       Turn verbose on (=show working)\n");
    printf("\n");
}

void get_earth_state( double epoch, char *ephemfile, vec *r_earth, vec *v_earth )
/* Extracts the Earth's position (km) and velocity (km/s) in J2000 coordinates
 * from the supplied planetary ephemeris.
 */
{
    // Open ephemeris file
    double jd = epoch + 2400000.5;
    SpiceDouble et = (jd - j2000_c() ) * spd_c();
    SpiceDouble state[6];
    SpiceDouble lt;
    furnsh_c( ephemfile );
    spkezr_c( "earth", et, "j2000", "NONE", "solar system barycenter",
              state, &lt );
    unload_c( ephemfile );

    // Put answer into "output" variables
    r_earth->x = state[0];   v_earth->x = state[3];
    r_earth->y = state[1];   v_earth->y = state[4];
    r_earth->z = state[2];   v_earth->z = state[5];
}


//void parseoptions( int argc, char *argv[] )
//{
//}


int main( int argc, char *argv[] )
{

    char   *par       = NULL;  // The name of the pulsar ephemeris par file
    char   *ephemfile = NULL;  // The name of the NASA planetary ephemeris file
    double  epoch     = 0.0;   // The epoch in question
    double  freq      = 0.0;   // Observing frequency
    int     verbose   = 0;     // 0 = verbose output off, 1 = verbose output on
    double  arccurve  = 0.0;   // The arc curvature, in s^3

    struct pardata pd;         // A container to hold the par file information

    // Parse command line options
    if (argc <= 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    int co;
    while ((co = getopt(argc, argv, "e:f:hp:s:v")) != -1)
    {
        switch (co)
        {
            case 'e':
                epoch = atof(optarg);
                break;
            case 'f':
                freq = atof(optarg);
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

    // Get Earth's position and velocity
    vec r_earth; // position vector (kms)
    vec v_earth; // velocity vector (km/s)
    get_earth_state( epoch, ephemfile, &r_earth, &v_earth );

    // Normalise the position and velocity vectors
    vec rn_earth, vn_earth;
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
        printf("  kom   = %.12f deg\n", pd.kom);
        printf("  dist  = %.12f kpc\n", pd.dist);
        printf("  pmra  = %.12f mas/yr\n", pd.pmra);
        printf("  pmdec = %.12f mas/yr\n", pd.pmdec);
        printf("  rajd  = %.12f deg\n", pd.rajd);
        printf("  decjd = %.12f deg\n", pd.decjd);
    }

    // Derive other values
    double cosi   = sqrt(1 - pd.sini*pd.sini);
    double komr   = DEG2RAD(pd.kom);
    double rajr   = DEG2RAD(pd.rajd);
    double decjr  = DEG2RAD(pd.decjd);
    double tvra   = PM2TV(pd.pmra, pd.dist);
    double tvdec  = PM2TV(pd.pmdec, pd.dist);

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
    double vearthx = ckomr*vI + skomr*vJ;
    double vearthy = skomr*vI - ckomr*vJ;

    printf( "\nVelocity of the Earth in 'x,y' coordinates (km/s):\n" );
    printf( "  [%lf, %lf]\n", vearthx, vearthy );





    // Calculate the orbital velocity vector
    vec ascnode;
    ascnode.x = ckomr;
    ascnode.y = skomr;
    ascnode.z = 0.0;

    vec ascnode1; // = ascnode cross [0,0,1]
    ascnode1.x =  ascnode.y;
    ascnode1.y = -ascnode.x;
    ascnode1.z = 0.0;

    vec n; // = cosi*[0,0,1] + sini*ascnode1
    n.x = pd.sini * ascnode1.x;
    n.y = pd.sini * ascnode1.y;
    n.z = cosi;

    double a     = pd.a1 * c / pd.sini;
    double theta = 2*PI*fmod( epoch - pd.t0, pd.pb ) / pd.pb;
    double r_a   = (1.0 - pd.ecc*pd.ecc) / (1.0 + pd.ecc*cos(theta));

    vec ascnode2; // = n cross ascnode1
    ascnode2.x = n.y * ascnode.z  -  n.z * ascnode.y;
    ascnode2.y = n.z * ascnode.x  -  n.x * ascnode.z;
    ascnode2.z = n.x * ascnode.y  -  n.y * ascnode.x;

    vec rhat;
    rhat.x = cos(pd.om + theta)*ascnode.x + sin(pd.om + theta)*ascnode2.x;
    rhat.y = cos(pd.om + theta)*ascnode.y + sin(pd.om + theta)*ascnode2.y;
    rhat.z = cos(pd.om + theta)*ascnode.z + sin(pd.om + theta)*ascnode2.z;

    vec vbinhat; // = n cross rhat
    vbinhat.x = n.y*rhat.z - n.z*rhat.y;
    vbinhat.y = n.z*rhat.x - n.x*rhat.z;
    vbinhat.z = n.x*rhat.y - n.y*rhat.x;

    double vbin = 2*PI*a / (pd.pb * 24.0 * 60.0 * 60.0 ) * sqrt(2.0/r_a - 1);

    double vbinra  = vbin * vbinhat.x;
    double vbindec = vbin * vbinhat.y;

    double vbinx = ckomr*vbinra + skomr*vbindec;
    double vbiny = skomr*vbinra - ckomr*vbindec;

    printf( "\n(Orbital) velocity of the pulsar in 'x,y' coordinates (km/s):\n" );
    printf( "  [%lf, %lf]\n", vbinx*1e-3, vbiny*1e-3 );



    // Calculate the proper motion (convert from mas/yr to km/s)
    double tvx = ckomr*tvra + skomr*tvdec;
    double tvy = skomr*tvra - ckomr*tvdec;

    printf( "\nProper motion of the pulsar in 'x,y' coordinates (km/s):\n" );
    printf( "  [%lf, %lf]\n", tvx, tvy );



    // Calculate what s must be


    // Free up memory
    free( par );
    free( ephemfile );
    free_par( &pd );

    return 0;
}
