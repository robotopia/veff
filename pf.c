#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_poly.h>

typedef struct output_parameters_t {
    char   dat_small_filename[110];
    char   png_filename[110];
    char   hough_filename[110];
    char   hough_png_filename[110];
    char   gpi_filename[110];
    int    display_dB;
    double best_a;
} output_parameters;

typedef struct input_parameters_t {
    char   dat_filename[100];
    int    x_orig, y_orig;
    double dx;
    double dy;
    double min_x, max_x;
    double min_y, max_y;
    double mask_o, mask_x, mask_y;
    double mask_ox, mask_oy;
    int    max_t;
    int    is_dB;
} input_parameters;

// Usage function
void usage()
{
    printf( "usage: parabfit [OPTIONS] DATAFILE\n" );
    printf( "\n" );
    printf( "  OPTIONS\n" );
    printf( "    --ssdata=PATH The file containing the secondary spectrum data. Should be a text file with three\n" );
    printf( "                  columns of number: \"xpixel ypixel value [default = ss.out]\"\n" );
    printf( "    --xorig=N     The x value (in units matching DATAFILE) corresponding to the origin [default = 0]\n" );
    printf( "    --yorig=N     The y value (in units matching DATAFILE) corresponding to the origin [default = 0]\n" );
    printf( "    --dx=N        The resolution of the x axis in mHz [default = 1.0]\n" );
    printf( "    --dy=N        The resolution of the y axis in us  [default = 1.0]\n" );
    printf( "    --xmin=N      The minimum x to consider (mHz)  [default = -20.0]\n" );
    printf( "    --xmax=N      The maximum x to consider (mHz)  [default =  20.0]\n" );
    printf( "    --ymin=N      The minimum y to consider (us)   [default =   0.0] (values < 0 will default to 0)\n" );
    printf( "    --ymax=N      The maximum y to consider (us)   [default =  20.0]\n" );
    printf( "    --omask=N     Data points <= this many pixels from the origin will be masked   [default =  20.0]\n" );
    printf( "    --xmask=N     Data points <= this many pixels from the x-axis will be masked   [default =  8.0]\n" );
    printf( "    --ymask=N     Data points <= this many pixels from the y-axis will be masked   [default =  8.0]\n" );
    printf( "    --tmax=N      The largest parabola thickness to consider (integer) [default = 5]\n" );
    printf( "    --dB          Input values are in dB units [default = off]\n" );
    printf( "    -h, --help    Display this help and exit\n" );
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

void write_gnuplot_script( FILE *f, input_parameters *ip, output_parameters *op )
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
    fprintf(f, "set xrange [%lf:%lf]\n", ip->min_x, ip->max_x);
    fprintf(f, "set yrange [%lf:%lf]\n", ip->min_y, ip->max_y);
    fprintf(f, "set cbrange [%f:%f]\n", -60.0, -38.0);
    fprintf(f, "set xlabel \"Doppler frequency (mHz)\"\n");
    fprintf(f, "set ylabel \"Delay (us)\"\n");
    fprintf(f, "set cblabel \"dB\"\n");
    fprintf(f, "set parametric\n");
    fprintf(f, "set dummy t\n");
    fprintf(f, "set trange [-pi:pi]\n");
    fprintf(f, "myx(t) = (t+pi)/(2*pi) * %lf + %lf\n", ip->max_x - ip->min_x, ip->min_x);
    fprintf(f, "myy(t) = (t+pi)/(2*pi) * %lf + %lf\n", ip->max_y - ip->min_y, ip->min_y);
    fprintf(f, "best_a = %lf\n", op->best_a);
    fprintf(f, "plot './%s' using (($1-%d)*%f):(($2-%d)*%f):3 with image notitle, \\\n",
                ip->dat_filename, ip->x_orig, ip->dx, ip->y_orig, ip->dy);
    fprintf(f, "     myx(t),best_a*(myx(t))**2 w l notitle lc rgb 'white', \\\n");
    fprintf(f, "     %lf*sin(t),%lf*cos(t) w l notitle lc rgb 'green', \\\n", ip->mask_ox, ip->mask_oy);
    fprintf(f, "     myx(t), %lf w l notitle lc rgb 'green', \\\n", ip->mask_y);
    fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green', \\\n", -ip->mask_x);
    fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green'\n\n", ip->mask_x);
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
    input_parameters ip;
    output_parameters op;

    // Set defaults for input parameters
    sprintf(ip.dat_filename, "ss.out");
    ip.x_orig = 0;
    ip.y_orig = 0;
    ip.dx = 1.0;
    ip.dy = 1.0;
    ip.min_x = -20.0;
    ip.max_x = -20.0;
    ip.min_y = 0.0;
    ip.max_y = 20.0;
    ip.mask_o = 20.0;
    ip.mask_x = 8.0;
    ip.mask_y = 8.0;
    ip.max_t = 5;
    ip.is_dB = 0;

    // Parse options

    int c;
    while (1)
    {

        static struct option long_options[] = {
            {"ssdata",  required_argument, NULL, 's'},
            {"xorig",   required_argument, NULL, 'o'},
            {"yorig",   required_argument, NULL, 'O'},
            {"dx",      required_argument, NULL, 'd'},
            {"dy",      required_argument, NULL, 'D'},
            {"xmin",    required_argument, NULL, 'x'},
            {"xmax",    required_argument, NULL, 'X'},
            {"ymin",    required_argument, NULL, 'y'},
            {"ymax",    required_argument, NULL, 'Y'},
            {"omask",   required_argument, NULL, 'k'},
            {"xmask",   required_argument, NULL, 'm'},
            {"ymask",   required_argument, NULL, 'M'},
            {"tmax",    required_argument, NULL, 'T'},
            {"dB",      no_argument,       NULL, 'L'},
            {"help",    no_argument,       NULL, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long(argc, argv, "d:D:hk:Lm:M:o:O:s:T:x:X:y:Y:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
            case 'd':
                ip.dx = atof(optarg);
                break;
            case 'D':
                ip.dy = atof(optarg);
                break;
            case 'h':
                usage();
                exit(0);
                break;
            case 'k':
                ip.mask_o = atof(optarg);
                break;
            case 'L':
                ip.is_dB = 1;
                break;
            case 'm':
                ip.mask_x = atof(optarg);
                break;
            case 'M':
                ip.mask_y = atof(optarg);
                break;
            case 'o':
                ip.x_orig = atof(optarg);
                break;
            case 'O':
                ip.y_orig = atof(optarg);
                break;
            case 's':
                sprintf(ip.dat_filename, optarg);
            case 'T':
                ip.max_t = atof(optarg);
                break;
            case 'x':
                ip.min_x = atof(optarg);
                break;
            case 'X':
                ip.max_x = atof(optarg);
                break;
            case 'y':
                ip.min_y = atof(optarg);
                break;
            case 'Y':
                ip.max_y = atof(optarg);
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
    if (ip.min_y < 0.0)  ip.min_y = 0.0;

    // Check that max's are more than min's
    if ((ip.max_x < ip.min_x) || (ip.max_y < ip.min_y))
    {
        fprintf(stderr,"error: minimum values cannot be greater than maximum values\n");
        exit(1);
    }

    // Convert mask_o to an ellipse with correct units
    ip.mask_ox = ip.mask_o * ip.dx;
    ip.mask_oy = ip.mask_o * ip.dy;

    // Round mask to pixel boundary
    ip.mask_x = floor(ip.mask_x) + 0.5;
    ip.mask_y = floor(ip.mask_y) + 0.5;

    // First, write out file containing only those data to be considered
    // Open the original file for reading (fr), and a file for writing (fw)
    sprintf(op.dat_small_filename, "%s.small", ip.dat_filename);
    FILE *fr = fopen(ip.dat_filename, "r");
    if (!fr)
    {
        fprintf(stderr,"error: Could not open file '%s'\n", ip.dat_filename);
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
        x = (x - (double)ip.x_orig) * ip.dx;
        y = (y - (double)ip.y_orig) * ip.dy;

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
        if ((x >= ip.min_x) && (x <= ip.max_x) &&
                (y >= ip.min_y) && (y <= ip.max_y))
        {
            fprintf(fw, "%lf %lf %le\n", x, y, val);
            n_pixels++;
        }
    }

    // Close files
    fclose(fr);
    fclose(fw);

    // If specified x and y limits are outside of actual data, trim them to fit actual data
    if (ip.max_x > temp_max_x)  ip.max_x = temp_max_x;
    if (ip.min_x < temp_min_x)  ip.min_x = temp_min_x;
    if (ip.max_y > temp_max_y)  ip.max_y = temp_max_y;
    if (ip.min_y < temp_min_y)  ip.min_y = temp_min_y;

    // Check that mask values are all positive
    if ((ip.mask_o < 0.0) || (ip.mask_x < 0.0) || (ip.mask_y < 0.0))
    {
        fprintf(stderr, "error: mask values must be > 0\n");
        exit(EXIT_FAILURE);
    }

    // Check that max_t is >= 1
    if (ip.max_t < 1)
    {
        fprintf(stderr, "error: max_t value must be >= 1\n");
        exit(EXIT_FAILURE);
    }

    // Convert mask_x and mask_y to correct units;
    ip.mask_x *= ip.dx;
    ip.mask_y *= ip.dy;

    // Find the minimum and maximum unmasked pixels in both x and y directions
    double inner_x, outer_x;
    double inner_y, outer_y;

    // For the x's...
    if (ip.min_x > 0)
    {
        inner_x = (ip.min_x > ip.mask_x ? ip.min_x : ip.mask_x);
        outer_x = ip.max_x;
    }
    else if (ip.max_x < 0)
    {
        inner_x = (fabs(ip.max_x) > ip.mask_x ? fabs(ip.max_x) : ip.mask_x);
        outer_x = fabs(ip.min_x);
    }
    else
    {
        inner_x = ip.mask_x;
        outer_x = fabs(ip.min_x) > ip.max_x ? fabs(ip.min_x) : ip.max_x;
    }

    // ...and for the y's
    if (ip.min_y > 0)
    {
        inner_y = (ip.min_y > ip.mask_y ? ip.min_y : ip.mask_y);
        outer_y = ip.max_y;
    }
    else if (ip.max_y < 0)
    {
        inner_y = (fabs(ip.max_y) > ip.mask_y ? fabs(ip.max_y) : ip.mask_y);
        outer_y = fabs(ip.min_y);
    }
    else
    {
        inner_y = ip.mask_y;
        outer_y = fabs(ip.min_y) > ip.max_y ? fabs(ip.min_y) : ip.max_y;
    }

    // Add a fudge factor, to avoid parabolas with only a minimal number of pixels
    inner_y *= 2.0;
    inner_x *= sqrt(2.0);

    // Set up array of trial 'a' values
    double step_x = ip.dx / 10.0;
    double step_y = ip.dy / 10.0;
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
    int n_ts = ip.max_t;
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
        if (ip.is_dB)  val = pow(10.0,val/10.0);

        if (fabs(x) <= ip.mask_x) // don't count pixels too close to y-axis
            continue;

        if (fabs(y) <= ip.mask_y) // don't count pixels too close to x-axis
            continue;

        orig_dist = hypot(x/ip.mask_ox, y/ip.mask_oy);
        if (orig_dist <= 1) // don't count pixels too close to origin
            continue;

        for (i = 0; i < n_as; i++)
        {
            a = trial_as[i];

            distparab(x, y, a, &parab_dist_x, &parab_dist_y); // Calculate distance away from parabola
            dist_pixels = hypot(parab_dist_x/ip.dx, parab_dist_y/ip.dy);

            for (j = 0; j < ip.max_t; j++)
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

    sprintf(op.png_filename,       "%s.png",       ip.dat_filename);
    sprintf(op.hough_filename,     "%s.hough.dat", ip.dat_filename);
    sprintf(op.hough_png_filename, "%s.hough.png", ip.dat_filename);
    sprintf(op.gpi_filename,       "%s.gpi",       ip.dat_filename);

    // Save Hough to file
    f = fopen(op.hough_filename, "w");
    printf("Writing hough file...\n");
    for (j = 0; j < n_ts; j++)
        for (i = 0; i < n_as; i++) // reverse order to make plots look right
        {
            if (i == 0)
                fprintf(f,"\n");
            if (ip.is_dB)
                hough[i][j] = 10.0*log10(hough[i][j]);
            if (!isnan(hough[i][j]))
                fprintf(f, "%e %f %e %d\n", trial_as[i], trial_ts[j], hough[i][j], pixelcount[i][j]);
        }
    fclose(f);

    // Write out gnuplot script
    f = fopen(op.gpi_filename, "w");
    printf("Writing gnuplot script...\n");
    write_gnuplot_script(f, &ip, &op);
    fclose(f);
}
