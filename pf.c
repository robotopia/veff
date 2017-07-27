#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_poly.h>
#include "ss.h"

typedef struct output_parameters_t {
    char   dat_small_filename[110];
    char   png_filename[110];
    char   hough_filename[110];
    char   hough_png_filename[110];
    char   gpi_filename[110];
    int    display_dB;
    double best_a;
} output_parameters;

// Usage function
void usage()
{
    printf( "usage: parabfit [OPTIONS] DATAFILE\n" );
    printf( "\n" );
    printf( "  OPTIONS\n" );
    printf( "    --ssdata=PATH   The file containing the secondary spectrum data. Should be a text file with three\n" );
    printf( "                    columns of number: \"xpixel ypixel value [default = ss.out]\"\n" );
    printf( "    --orig=X,Y      The X,Y values (in units matching DATAFILE) corresponding to the origin [default = 0,0]\n" );
    printf( "    --res=X,Y       The resolution of the x/y axis in mHz/us [default = (1.0, 1.0)]\n" );
    printf( "    --xrange=LO,HI  The min/max X to consider (mHz)  [default = (-20.0,20.0)]\n" );
    printf( "    --yrange=LO,HI  The min/max Y to consider (us)   [default = (0.0,20.0)] (values < - will default to 0)\n" );
    printf( "    --mask=O,X,Y    Data points <= this many pixels from the origin/x-axis/y-axis will be masked   [default = (20.0, 8.0, 8.0)]\n" );
    printf( "    --tmax=N        The largest parabola thickness to consider (integer) [default = 5]\n" );
    printf( "    --dB            Input values are in dB units [default = off]\n" );
    printf( "    -h, --help      Display this help and exit\n" );
    printf( "\n" );
}

void distparab(double x, double y, double a, double *dx, double *dy)
{
  x = -fabs(x);
  double A = 1.0/(2.0*a*a) - y/a;
  double B =  -x/(2.0*a*a);

  double r2, r3; // Placeholders for the roots of the cubic equation
  double xp, yp;

  gsl_poly_solve_cubic( 0.0, A, B, &xp, &r2, &r3 );
  yp = a*xp*xp;

  *dx = xp - x;
  *dy = yp - y;
}

void write_ss_gnuplot_script( FILE *f, struct sec_spect *ss, output_parameters *op )
{
    // Write gnuplot script
    fprintf(f, "print \"Running gnuplot script...\"\n");
    fprintf(f, "set terminal pngcairo size 800,600\n\n");

    fprintf(f, "set out '%s'\n", op->hough_png_filename);
    fprintf(f, "set xlabel \"'a' parameter\"\n");
    if (op->display_dB)
        fprintf(f, "set ylabel \"Mean dB\"\n");
    else
        fprintf(f, "set ylabel \"Mean power\"\n");
    fprintf(f, "set cblabel \"Thickness of parabola (pixels)\"\n");
    fprintf(f, "set logscale x\n");
    fprintf(f, "plot './%s' u 1:3:2 w l palette notitle \n\n", op->hough_filename);

    fprintf(f, "set out '%s'\n", op->png_filename);
    fprintf(f, "unset logscale x\n");
    fprintf(f, "set dgrid3d\n");
    fprintf(f, "set samples 10000\n");
    fprintf(f, "unset xlabel\n");
    fprintf(f, "unset ylabel\n");
    fprintf(f, "set xrange [%lf:%lf]\n", ss->min_x, ss->max_x);
    fprintf(f, "set yrange [%lf:%lf]\n", ss->min_y, ss->max_y);
    fprintf(f, "set cbrange [%f:%f]\n", -60.0, -38.0);
    fprintf(f, "set xlabel \"Doppler frequency (mHz)\"\n");
    fprintf(f, "set ylabel \"Delay (us)\"\n");
    fprintf(f, "set cblabel \"dB\"\n");
    fprintf(f, "set parametric\n");
    fprintf(f, "set dummy t\n");
    fprintf(f, "set trange [-pi:pi]\n");
    fprintf(f, "myx(t) = (t+pi)/(2*pi) * %lf + %lf\n", ss->max_x - ss->min_x, ss->min_x);
    fprintf(f, "myy(t) = (t+pi)/(2*pi) * %lf + %lf\n", ss->max_y - ss->min_y, ss->min_y);
    fprintf(f, "best_a = %lf\n", op->best_a);
    fprintf(f, "plot './%s' using (($1-%d)*%f):(($2-%d)*%f):3 with image notitle, \\\n",
                ss->dat_filename, ss->x_orig, ss->dx, ss->y_orig, ss->dy);
    fprintf(f, "     myx(t),best_a*(myx(t))**2 w l notitle lc rgb 'white', \\\n");
    fprintf(f, "     %lf*sin(t),%lf*cos(t) w l notitle lc rgb 'green', \\\n", ss->mask_ox, ss->mask_oy);
    fprintf(f, "     myx(t), %lf w l notitle lc rgb 'green', \\\n", ss->mask_y);
    fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green', \\\n", -ss->mask_x);
    fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green'\n\n", ss->mask_x);
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
    output_parameters op;

    // Set defaults for input parameters
    sprintf(ss.dat_filename, "ss.out");
    ss.x_orig = 0;
    ss.y_orig = 0;
    ss.dx = 1.0;
    ss.dy = 1.0;
    ss.min_x = -20.0;
    ss.max_x = -20.0;
    ss.min_y = 0.0;
    ss.max_y = 20.0;
    ss.mask_o = 20.0;
    ss.mask_x = 8.0;
    ss.mask_y = 8.0;
    ss.max_t = 5;
    ss.is_dB = 0;

    // Parse options

    int nscan; // Number of items scanned in with sscanf
    int c;
    while (1)
    {

        static struct option long_options[] = {
            {"ssdata",  required_argument, NULL, 's'},
            {"orig",    required_argument, NULL, 'o'},
            {"res",     required_argument, NULL, 'r'},
            {"xrange",  required_argument, NULL, 'x'},
            {"yrange",  required_argument, NULL, 'y'},
            {"mask",    required_argument, NULL, 'm'},
            {"tmax",    required_argument, NULL, 't'},
            {"dB",      no_argument,       NULL, 'l'},
            {"help",    no_argument,       NULL, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "hlm:o:r:s:t:x:y:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
                usage();
                exit(0);
                break;
            case 'l':
                ss.is_dB = 1;
                break;
            case 'm':
                nscan = sscanf(optarg, "%lf,%lf,%lf", &ss.mask_o, &ss.mask_x, &ss.mask_y);
                if (nscan != 3)
                {
                    fprintf(stderr, "error: couldn't parse --mask=%s as O,X,Y\n", optarg);
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
            case 't':
                ss.max_t = atof(optarg);
                break;
            case 'x':
                nscan = sscanf(optarg, "%lf,%lf", &ss.min_x, &ss.max_x);
                if (nscan != 2)
                {
                    fprintf(stderr, "error: couldn't parse --xrange=%s as LO,HI\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            case 'y':
                nscan = sscanf(optarg, "%lf,%lf", &ss.min_y, &ss.max_y);
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

    // Check: min_y >= 0
    if (ss.min_y < 0.0)  ss.min_y = 0.0;

    // Check that max's are more than min's
    if ((ss.max_x < ss.min_x) || (ss.max_y < ss.min_y))
    {
        fprintf(stderr,"error: minimum values cannot be greater than maximum values\n");
        exit(1);
    }

    // Convert mask_o to an ellipse with correct units
    ss.mask_ox = ss.mask_o * ss.dx;
    ss.mask_oy = ss.mask_o * ss.dy;

    // Round mask to pixel boundary
    ss.mask_x = floor(ss.mask_x) + 0.5;
    ss.mask_y = floor(ss.mask_y) + 0.5;

    // First, write out file containing only those data to be considered
    // Open the original file for reading (fr), and a file for writing (fw)
    sprintf(op.dat_small_filename, "%s.small", ss.dat_filename);
    FILE *fr = fopen(ss.dat_filename, "r");
    if (!fr)
    {
        fprintf(stderr,"error: Could not open file '%s'\n", ss.dat_filename);
        exit(EXIT_FAILURE);
    }
    FILE *fw = fopen(op.dat_small_filename, "w");
    if (!fw)
    {
        fprintf(stderr,"error: Could not open file '%s'\n", op.dat_small_filename);
        exit(EXIT_FAILURE);
    }

    // Fields are x,y,val
    double x, y;
    double temp_max_x, temp_min_x;
    double temp_max_y, temp_min_y;
    double val;
    int is_first_time = 1;
    int n_pixels = 0;

    while (!feof(fr))
    {
        // Read in values
        fscanf(fr, "%lf %lf %lf", &x, &y, &val);

        // Convert pixel numbers to correct units
        x = (x - (double)ss.x_orig) * ss.dx;
        y = (y - (double)ss.y_orig) * ss.dy;

        if (is_first_time)
        {
            temp_max_x = temp_min_x = x;
            temp_max_y = temp_min_y = x;
            is_first_time = 0;
        }

        // Keep track of the max and min values of x and y read in so far
        if (x > temp_max_x)  temp_max_x = x;
        if (x < temp_min_x)  temp_min_x = x;
        if (y > temp_max_y)  temp_max_y = y;
        if (y < temp_min_y)  temp_min_y = y;

        // If values are within range, write out to "reduced" file
        if ((x >= ss.min_x) && (x <= ss.max_x) &&
                (y >= ss.min_y) && (y <= ss.max_y))
        {
            fprintf(fw, "%lf %lf %le\n", x, y, val);
            n_pixels++;
        }
    }

    // Close files
    fclose(fr);
    fclose(fw);

    // If specified x and y limits are outside of actual data, trim them to fit actual data
    if (ss.max_x > temp_max_x)  ss.max_x = temp_max_x;
    if (ss.min_x < temp_min_x)  ss.min_x = temp_min_x;
    if (ss.max_y > temp_max_y)  ss.max_y = temp_max_y;
    if (ss.min_y < temp_min_y)  ss.min_y = temp_min_y;

    // Check that mask values are all positive
    if ((ss.mask_o < 0.0) || (ss.mask_x < 0.0) || (ss.mask_y < 0.0))
    {
        fprintf(stderr, "error: mask values must be > 0\n");
        exit(EXIT_FAILURE);
    }

    // Check that max_t is >= 1
    if (ss.max_t < 1)
    {
        fprintf(stderr, "error: max_t value must be >= 1\n");
        exit(EXIT_FAILURE);
    }

    // Convert mask_x and mask_y to correct units;
    ss.mask_x *= ss.dx;
    ss.mask_y *= ss.dy;

    // Find the minimum and maximum unmasked pixels in both x and y directions
    double inner_x, outer_x;
    double inner_y, outer_y;

    // For the x's...
    if (ss.min_x > 0)
    {
        inner_x = (ss.min_x > ss.mask_x ? ss.min_x : ss.mask_x);
        outer_x = ss.max_x;
    }
    else if (ss.max_x < 0)
    {
        inner_x = (fabs(ss.max_x) > ss.mask_x ? fabs(ss.max_x) : ss.mask_x);
        outer_x = fabs(ss.min_x);
    }
    else
    {
        inner_x = ss.mask_x;
        outer_x = fabs(ss.min_x) > ss.max_x ? fabs(ss.min_x) : ss.max_x;
    }

    // ...and for the y's
    if (ss.min_y > 0)
    {
        inner_y = (ss.min_y > ss.mask_y ? ss.min_y : ss.mask_y);
        outer_y = ss.max_y;
    }
    else if (ss.max_y < 0)
    {
        inner_y = (fabs(ss.max_y) > ss.mask_y ? fabs(ss.max_y) : ss.mask_y);
        outer_y = fabs(ss.min_y);
    }
    else
    {
        inner_y = ss.mask_y;
        outer_y = fabs(ss.min_y) > ss.max_y ? fabs(ss.min_y) : ss.max_y;
    }

    // Add a fudge factor, to avoid parabolas with only a minimal number of pixels
    inner_y *= 2.0;
    inner_x *= sqrt(2.0);

    // Set up array of trial 'a' values
    double step_x = ss.dx / 10.0;
    double step_y = ss.dy / 10.0;
    int n_as_x = (int)((outer_x - inner_x) / step_x) + 1;
    int n_as_y = (int)((outer_y - inner_y) / step_y) + 1;
    int n_as = n_as_x + n_as_y + 2; // +2 just a safety buffer in case something goes wrong with the int casting above
    int n = 0;
    double trial_as[n_as];
    for (y = inner_y; y <= outer_y; y += step_y) // Run up pixels bordering on right
    {
        trial_as[n] = y / (outer_x*outer_x);
        n++;
    }

    for (x = outer_x; x >= inner_x; x -= step_x)
    {
        trial_as[n] = outer_y / (x*x);
        n++;
        if (n > n_as)
        {
            fprintf(stderr, "error: more trial 'a's written (%d) than space allocated (%d), ", n, n_as);
            fprintf(stderr, "inner_x = %f, x = %f, step_x = %f\n", inner_x, x, step_x);
            exit(EXIT_FAILURE);
        }
    }

    // Set up array of trial 't' (thickness) values
    int n_ts = ss.max_t;
    double trial_ts[n_ts];
    for (i = 0; i < n_ts; i++)
        trial_ts[i] = (double)(i+1);

    // Create the "Hough" transform with parameters "a" and "t"
    double a, t; // 'a' parameter & thickness 't'

    double parab_dist_x; // Shortest distance from point to parabola (x component)
    double parab_dist_y; // Shortest distance from point to parabola (y component)
    double dist_pixels; // The same distance as above in pixel units
    double orig_dist; // Distance from pixel to origin
    int pixelcount[n_as][n_ts]; // Count how many pixels are in this sum
    double hough[n_as][n_ts]; // The actual result of the transform
    int old_percent_done = 0, percent_done = 0;

    // Initialse hough and pixelcount to zero
    for (i = 0; i < n_as; i++)
        for (j = 0; j < n_ts; j++)
        {
            hough[i][j] = 0.0;
            pixelcount[i][j] = 0;
        }

    // Load the reduced file
    FILE *f = fopen(op.dat_small_filename, "r");
    if (!f)
    {
        fprintf(stderr,"error: Could not open file '%s'\n", op.dat_small_filename);
        exit(EXIT_FAILURE);
    }

    printf("Calculating Hough transform..."); fflush(stdout);
    printf("%4d%%", percent_done);       fflush(stdout);

    // Go through the data and sum up the eligible points
    n = 0; // Just for counting how many pixels have been processed, to report percentage done
    while (!feof(f))
    {
        // Output percentage done
        n++;
        percent_done = n * 100 / n_pixels;
        if (percent_done > old_percent_done)
        {
            printf("\b\b\b\b\b%4d%%", percent_done); fflush(stdout);
            old_percent_done = percent_done;
        }

        // Get next lot of values
        fscanf(f, "%lf %lf %lf\n", &x, &y, &val);

        // Convert values to real, not-log powers, if necessary
        if (ss.is_dB)  val = pow(10.0,val/10.0);

        if (fabs(x) <= ss.mask_x) // don't count pixels too close to y-axis
            continue;

        if (fabs(y) <= ss.mask_y) // don't count pixels too close to x-axis
            continue;

        orig_dist = hypot(x/ss.mask_ox, y/ss.mask_oy);
        if (orig_dist <= 1) // don't count pixels too close to origin
            continue;

        for (i = 0; i < n_as; i++)
        {
            a = trial_as[i];

            distparab(x, y, a, &parab_dist_x, &parab_dist_y); // Calculate distance away from parabola
            dist_pixels = hypot(parab_dist_x/ss.dx, parab_dist_y/ss.dy);

            for (j = 0; j < ss.max_t; j++)
            {
                t = trial_ts[j];

                if (dist_pixels <= t && x < 0.0)
                {
                    pixelcount[i][j]++;
                    hough[i][j] += val;
                }
            } // end for thickness 't'
        } // end for parameter 'a'
    }

    fclose(f);

    // Just print a new line
    printf("\n");

    // Use the mean value rather than the sum (i.e. normalise w.r.t number of pixels summed)
    printf("Finding best 'a' parameter... \n");
    op.best_a = trial_as[0];     // Initialise best 'a' parameter to the first trial 'a'
    double best_hough = hough[0][0] / (double)pixelcount[0][0]; // Initialise best hough value to the first hough value
    for (i = 0; i < n_as; i++)
        for (j = 0; j < n_ts; j++)
        {
            hough[i][j] = hough[i][j] / (double)pixelcount[i][j];
            if (hough[i][j] > best_hough)
            {
                best_hough = hough[i][j];
                op.best_a  = trial_as[i];
            }
        }
    printf("  Found a = %f\n", op.best_a);

    ////////////////////
    // OUTPUT RESULTS //
    ////////////////////

    sprintf(op.png_filename,       "%s.png",       ss.dat_filename);
    sprintf(op.hough_filename,     "%s.hough.dat", ss.dat_filename);
    sprintf(op.hough_png_filename, "%s.hough.png", ss.dat_filename);
    sprintf(op.gpi_filename,       "%s.gpi",       ss.dat_filename);

    // Save Hough to file
    f = fopen(op.hough_filename, "w");
    printf("Writing hough file...\n");
    for (j = 0; j < n_ts; j++)
        for (i = 0; i < n_as; i++) // reverse order to make plots look right
        {
            if (i == 0)
                fprintf(f,"\n");
            if (ss.is_dB)
                hough[i][j] = 10.0*log10(hough[i][j]);
            if (!isnan(hough[i][j]))
                fprintf(f, "%e %f %e %d\n", trial_as[i], trial_ts[j], hough[i][j], pixelcount[i][j]);
        }
    fclose(f);

    // Write out gnuplot script
    f = fopen(op.gpi_filename, "w");
    printf("Writing gnuplot script...\n");
    write_ss_gnuplot_script(f, &ss, &op);
    fclose(f);
}
