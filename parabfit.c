#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "ss.h"
#include "hough.h"

// Usage function
void usage()
{
    printf( "usage: parabfit [OPTIONS] DATAFILE\n" );
    printf( "\n" );
    printf( "  DATAFILE\n" );
    printf( "    The file containing the secondary " );
    printf( "spectrum data. Should be a text file with three\n" );
    printf( "    columns of number: \"xpixel ypixel value\"\n\n" );
    printf( "  OPTIONS\n" );
    printf( "    --orig=X,Y      The X,Y values (in units matching " );
    printf( "DATAFILE) corresponding to the origin [default = 0,0]\n" );
    printf( "    --units=X,Y     Units of the secondary spectrum's " );
    printf( "x/y-axes   [default = 0,0]\n" );
    printf( "    --res=X,Y       The resolution of the x/y axis " );
    printf( "[default = (1.0, 1.0)]\n" );
    printf( "    --xrange=LO,HI  The min/max X to consider (mHz)  " );
    printf( "[default = (-20.0,20.0)]\n" );
    printf( "    --yrange=LO,HI  The min/max Y to consider (us)   " );
    printf( "[default = (0.0,20.0)] (values < - will default to 0)\n" );
    printf( "    --mask=X,Y      Data points <= this distance from the " );
    printf( "x/y-axis will be masked   " );
    printf( "[default = (1.5, 1.5)]\n" );
    printf( "    --omask=OX,OY   Data points <= this ellipse at the " );
    printf( "origin will be masked   " );
    printf( "[default = (4.0, 4.0)]\n" );
    printf( "    --cbrange=MIN,MAX  Set dynamic range of output\n" );
    printf( "    --pdist=DX,DY   Parabola \"thickness\" to consider " );
    printf( "[default = 1.0,1.0]\n" );
    printf( "    --dB            Input values are in dB units " );
    printf( "[default = off]\n" );
    printf( "    --logspace      Sample curvatures evenly in log space" );
    printf( "[default = off]\n" );
    printf( "    --curves=START:END:N  Sample N curvatures in the Hough " );
    printf( "Transform, from START to END  " );
    printf( "[default = 0.1:10.0:1000]\n" );
    printf( "    --q1, --q2, --q3, --q4\n" );
    printf( "                    Sum over pixels in quadrants 1,2,3,4. " );
    printf( "If any are supplied, the defaults are not used  " );
    printf( "[default = --q1 and --q2]\n" );
    printf( "    --out=FILE      File to write the transform results to  " );
    printf( "    --hggpi=FILE    File to write a gnuplot script for " );
    printf( "viewing the Hough transform to    " );
    printf( "    --crop=FILE     File to write the cropped solution to   " );
    printf( "    --ssgpi=FILE    File to write a gnuplot script for " );
    printf( "viewing the (cropped) secondary spectrum to    " );
    printf( "    -h, --help      Display this help and exit\n" );
    printf( "\n" );
}

/********
 * MAIN *
 *******/

int main( int argc, char *argv[] )
{
    // Check if option -h was given
    if (argc == 2 && strcmp(argv[1],"-h") == 0)
    {
        usage();
        exit(EXIT_FAILURE);
    }

    // Generic counters
    int i,j;

    // Get values from input file
    struct sec_spect ss;
    struct hough hg;
    char ss_filename[256];
    char *crop  = NULL;
    char *out   = NULL;
    char *ssgpi = NULL;
    char *hggpi = NULL;

    // Set defaults input values
    sprintf(ss_filename, "ss.out");
    ss.x_orig = 0;
    ss.y_orig = 0;
    ss.dx = 1.0;
    ss.dy = 1.0;
    sprintf(ss.xunits, "mHz");
    sprintf(ss.yunits, "μs");
    ss.is_dB = 0;
    ss.cbmin = nan;
    ss.cbmax = nan;

    double min_x = -20.0;
    double max_x = -20.0;
    double min_y = 0.0;
    double max_y = 20.0;

    hg.quadrant = 0; // i.e. no quadrants; empty
    hg.amin     = 0.1;
    hg.amax     = 10.0;
    hg.x0mask   = 1.5;
    hg.y0mask   = 1.5;
    hg.xmask    = 4.0;
    hg.ymask    = 4.0;
    hg.pxdist   = 1.0;
    hg.pydist   = 1.0;
    hg.logspace = 0;
    hg.size     = 1000;

    // Parse options

    int nscan; // Number of items scanned in with sscanf
    int c;
    while (1)
    {

        static struct option long_options[] = {
            {"ssdata",   required_argument, NULL, 's'},
            {"orig",     required_argument, NULL, 'o'},
            {"res",      required_argument, NULL, 'r'},
            {"xrange",   required_argument, NULL, 'x'},
            {"yrange",   required_argument, NULL, 'y'},
            {"mask",     required_argument, NULL, 'm'},
            {"omask",    required_argument, NULL, 'M'},
            {"cbrange",  required_argument, NULL, 'c'},
            {"dB",       no_argument,       NULL, 'l'},
            {"logspace", no_argument,       NULL, 'L'},
            {"curves",   required_argument, NULL, 'n'},
            {"q1",       no_argument,       NULL, 'p'},
            {"q2",       no_argument,       NULL, 'P'},
            {"q3",       no_argument,       NULL, 'q'},
            {"q4",       no_argument,       NULL, 'Q'},
            {"help",     no_argument,       NULL, 'h'},
            {"crop",     required_argument, NULL, 'A'},
            {"out",      required_argument, NULL, 'B'},
            {"ssgpi",    required_argument, NULL, 'C'},
            {"hggpi",    required_argument, NULL, 'D'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "A:B:c:C:D:hlLm:n:o:pPqQr:s:x:y:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
            case 'A':
                crop = strdup(optarg);
                break;
            case 'B':
                out = strdup(optarg);
                break;
            case 'c':
                nscan = sscanf(optarg, "%lf,%lf", &hg.cbmin, &hg.cbmax);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --cbrange=%s as MIN,MAX\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'C':
                ssgpi = strdup(optarg);
                break;
            case 'D':
                hggpi = strdup(optarg);
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'l':
                ss.is_dB = 1;
                break;
            case 'L':
                hg.logspace = 1;
                break;
            case 'm':
                nscan = sscanf(optarg, "%lf,%lf", &hg.mask_x, &hg.mask_y);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --mask=%s as X,Y\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'M':
                nscan = sscanf(optarg, "%lf,%lf", &hg.x0mask, &hg.y0mask);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --omask=%s as DX,DY\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'n':
                nscan = sscanf(optarg, "%lf:%lf:%lf", &hg.amin, &hg.amax, &hg.size);
                if (nscan != 3)
                {
                    fprintf(stderr, "error: couldn't parse --curves=%s as START:END:N\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                nscan = sscanf(optarg, "%lf,%lf", &ss.x_orig, &ss.y_orig);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --orig=%s as X,Y\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'p':
                hg.quadrant |= HG_Q1;
                break;
            case 'P':
                hg.quadrant |= HG_Q2;
                break;
            case 'q':
                hg.quadrant |= HG_Q3;
                break;
            case 'Q':
                hg.quadrant |= HG_Q4;
                break;
            case 'r':
                nscan = sscanf(optarg, "%lf,%lf", &ss.dx, &ss.dy);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --res=%s as X,Y\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 's':
                sprintf(ss.dat_filename, optarg);
            case 'x':
                nscan = sscanf(optarg, "%lf,%lf", &min_x, &max_x);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --xrange=%s as LO,HI\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'y':
                nscan = sscanf(optarg, "%lf,%lf", &min_y, &max_y);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --yrange=%s as LO,HI\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
                break;
            default:
                fprintf(stderr, "error: unrecognised option\n");
                exit(EXIT_FAILURE);
                break;
        }
    }

    // Check if any quadrants were supplied
    if (hg.quadrant == 0)
        hg.quadrant = HG_Q1 | HG_Q2;

    // Check that any gnuplot requests have an accompanying
    // data-out request
    if (crop == NULL && ssgpi != NULL)
    {
        fprintf( stderr, "error: you must supply a --crop output in order " );
        fprintf( stderr, "to be able to write a --ssgpi gnuplot script\n" );
        exit(EXIT_FAILURE);
    }
    if (out == NULL && hggpi != NULL)
    {
        fprintf( stderr, "error: you must supply an --out output in order " );
        fprintf( stderr, "to be able to write a --hggpi gnuplot script\n" );
        exit(EXIT_FAILURE);
    }

    // Read in secondary spectrum
    if (optind >= argc)
    {
        fprintf( stderr, "error: DATAFILE required\n" );
        usage();
        exit(EXIT_FAILURE);
    }

    sprintf(ss_filename, "%s", argv[optind]);
    FILE *f = fopen( ss_filename, "r" );
    ss_read( f, &ss, SS_ASCII );
    fclose(f);

    // Crop it to the desired dimensions
    struct sec_spect ss_cropped;
    ss_crop( &ss, &ss_cropped, min_x, max_x, min_y, max_y );

    // Associate the hough transform with the cropped ss
    hg.ss = &ss_cropped;

    // Hough Transform
    hg_calc_transform( &hg );

    // Write out the results
    if (crop)
    {
        f = fopen( crop, "w" );
        ss_write( f, &ss_cropped, SS_BINARY );
        fclose(f);
    }
    if (out)
    {
        f = fopen( out, "w" );
        hg_write( f, &hg, HG_ASCII );
        fclose(f);
    }
    if (ssgpi)
    {
        f = fopen( ssgpi, "w" );
        ss_write_gnuplot( f, &ss_cropped, crop, 1 );
        fclose(f);
    }
    if (hggpi)
    {
        f = fopen( hggpi, "w" );
        ss_write( f, &hg, out );
        fclose(f);
    }

    // Free memory
    ss_free( &ss );
    ss_free( &ss_cropped );
}
