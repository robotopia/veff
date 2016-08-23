#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_poly.h>

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

int main( int argc, char *argv[] )
{
  // Check if option -h was given
  if (argc == 2 && strcmp(argv[1],"-h") == 0)
  {
    fprintf(stderr,"usage: %s\n", argv[0]);
    fprintf(stderr,"  The file 'pf.model' must be present and it must contain the following, each on separate lines, in this order:\n\n");
    fprintf(stderr,"    inputfile      - the ASCII file containing columns pixelx, pixely, value\n");
    fprintf(stderr,"    x_orig y_orig  - the x,y values (in units matching the inputfile) corresponding to the origin\n");
    fprintf(stderr,"    dx dy          - the resolution of the x,y axes in mHz, us\n");
    fprintf(stderr,"    min_x max_x    - the range of x's to consider (mHz)\n");
    fprintf(stderr,"    min_y max_y    - the range of y's to consider (us); values < 0 will default to 0\n");
    fprintf(stderr,"    mask_o mask_x mask_y - data points <= than this many pixels from origin/y_axis/x_axis will be masked\n");
    fprintf(stderr,"    max_t          - the largest parabola thickness to consider (integer >= 1)\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  // Generic counters
  int i,j;

  // Get values from input file
  char dat_filename[100];
  int x_orig, y_orig;
  double dx;
  double dy;
  double min_x, max_x;
  double min_y, max_y;
  double mask_o, mask_x, mask_y;
  int max_t;

  double mask_ox, mask_oy;

  const char *model_filename = "pf.model";
  FILE *f = fopen(model_filename, "r");
  if (!f)
  {
    fprintf(stderr,"error: Could not open file '%s'\n", model_filename);
    exit(1);
  }

  fscanf(f, "%s", dat_filename);
  fscanf(f, "%d %d", &x_orig, &y_orig);
  fscanf(f, "%lf %lf", &dx, &dy);
  fscanf(f, "%lf %lf", &min_x, &max_x);
  fscanf(f, "%lf %lf", &min_y, &max_y);
  fscanf(f, "%lf %lf %lf", &mask_o, &mask_x, &mask_y);
  fscanf(f, "%d", &max_t);

  fclose(f);

  // Check: min_y >= 0
  if (min_y < 0.0)  min_y = 0.0;

  // Check that max's are more than min's
  if ((max_x < min_x) || (max_y < min_y))
  {
    fprintf(stderr,"error: minimum values cannot be greater than maximum values\n");
    exit(1);
  }

  // Convert mask_o to an ellipse with correct units
  mask_ox = mask_o * dx;
  mask_oy = mask_o * dy;

  // Round mask to pixel boundary
  mask_x = floor(mask_x) + 0.5;
  mask_y = floor(mask_y) + 0.5;

  // First, write out file containing only those data to be considered
  // Open the original file for reading (fr), and a file for writing (fw)
  char dat_small_filename[110];
  sprintf(dat_small_filename, "%s.small", dat_filename);
  FILE *fr = fopen(dat_filename, "r");
  if (!fr)
  {
    fprintf(stderr,"error: Could not open file '%s'\n", dat_filename);
    exit(1);
  }
  FILE *fw = fopen(dat_small_filename, "w");
  if (!fw)
  {
    fprintf(stderr,"error: Could not open file '%s'\n", dat_small_filename);
    exit(1);
  }

  // Fields are x,y,val
  double x, y;
  double temp_max_x, temp_min_x;
  double temp_max_y, temp_min_y;
  double val;
  int is_dB = 1; // indicates power values are log (dB) values
  int is_first_time = 1;
  int n_pixels = 0;

  while (!feof(fr))
  {
    // Read in values
    fscanf(fr, "%lf %lf %lf", &x, &y, &val);

    // Convert pixel numbers to correct units
    x = (x - (double)x_orig) * dx;
    y = (y - (double)y_orig) * dy;

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
    if ((x >= min_x) && (x <= max_x) && (y >= min_y) && (y <= max_y))
    {
      fprintf(fw, "%lf %lf %le\n", x, y, val);
      n_pixels++;
    }
  }

  // Close files
  fclose(fr);
  fclose(fw);

  // If specified x and y limits are outside of actual data, trim them to fit actual data
  if (max_x > temp_max_x)  max_x = temp_max_x;
  if (min_x < temp_min_x)  min_x = temp_min_x;
  if (max_y > temp_max_y)  max_y = temp_max_y;
  if (min_y < temp_min_y)  min_y = temp_min_y;

  // Check that mask values are all positive
  if ((mask_o < 0.0) || (mask_x < 0.0) || (mask_y < 0.0))
  {
    fprintf(stderr, "error: mask values must be > 0\n");
    exit(1);
  }

  // Check that max_t is >= 1
  if (max_t < 1)
  {
    fprintf(stderr, "error: max_t value must be >= 1\n");
    exit(1);
  }

  // Convert mask_x and mask_y to correct units;
  mask_x *= dx;
  mask_y *= dy;

  // Find the minimum and maximum unmasked pixels in both x and y directions
  double inner_x, outer_x;
  double inner_y, outer_y;

  // For the x's...
  if (min_x > 0)
  {
    inner_x = (min_x > mask_x ? min_x : mask_x);
    outer_x = max_x;
  }
  else if (max_x < 0)
  {
    inner_x = (fabs(max_x) > mask_x ? fabs(max_x) : mask_x);
    outer_x = fabs(min_x);
  }
  else
  {
    inner_x = mask_x;
    outer_x = fabs(min_x) > max_x ? fabs(min_x) : max_x;
  }

  // ...and for the y's
  if (min_y > 0)
  {
    inner_y = (min_y > mask_y ? min_y : mask_y);
    outer_y = max_y;
  }
  else if (max_y < 0)
  {
    inner_y = (fabs(max_y) > mask_y ? fabs(max_y) : mask_y);
    outer_y = fabs(min_y);
  }
  else
  {
    inner_y = mask_y;
    outer_y = fabs(min_y) > max_y ? fabs(min_y) : max_y;
  }

  // Add a fudge factor, to avoid parabolas with only a minimal number of pixels
  inner_y *= 2.0;
  inner_x *= sqrt(2.0);

  // Set up array of trial 'a' values
  double step_x = dx / 10.0;
  double step_y = dy / 10.0;
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

int m = 0;
  for (x = outer_x; x >= inner_x; x -= step_x)
  {
    trial_as[n] = outer_y / (x*x);
    n++;
    if (n > n_as)
    {
      fprintf(stderr,"error: more trial 'a's written (%d) than space allocated (%d), inner_x = %f, x = %f, step_x = %f\n", n, n_as, inner_x, x, step_x);
      exit(1);
    }
m++;
  }

  // Set up array of trial 't' (thickness) values
  int n_ts = max_t;
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
  f = fopen(dat_small_filename, "r");
  if (!f)
  {
    fprintf(stderr,"error: Could not open file '%s'\n", dat_filename);
    exit(1);
  }

  printf("Calculating Hough transform..."); fflush(stdout);
  printf("%4d%% done", percent_done);       fflush(stdout);

  // Go through the data and sum up the eligible points
  n = 0; // Just for counting how many pixels have been processed, to report percentage done
  while (!feof(f))
  {
    // Output percentage done
    n++;
    percent_done = n * 100 / n_pixels;
    if (percent_done > old_percent_done)
    {
      printf("\b\b\b\b\b\b\b\b\b\b%4d%% done", percent_done); fflush(stdout);
      old_percent_done = percent_done;
    }

    // Get next lot of values
    fscanf(f, "%lf %lf %lf\n", &x, &y, &val);

    // Convert values to real, not-log powers, if necessary
    if (is_dB)  val = pow(10.0,val/10.0);

    if (fabs(x) <= mask_x) // don't count pixels too close to y-axis
      continue;

    if (fabs(y) <= mask_y) // don't count pixels too close to x-axis
      continue;

    orig_dist = hypot(x/mask_ox, y/mask_oy);
    if (orig_dist <= 1) // don't count pixels too close to origin
      continue;

    for (i = 0; i < n_as; i++)
    {
      a = trial_as[i];

      distparab(x, y, a, &parab_dist_x, &parab_dist_y); // Calculate distance away from parabola
      dist_pixels = hypot(parab_dist_x/dx, parab_dist_y/dy);

      for (j = 0; j < max_t; j++)
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
  double best_a = trial_as[0];     // Initialise best 'a' parameter to the first trial 'a'
  double best_hough = hough[0][0] / (double)pixelcount[0][0]; // Initialise best hough value to the first hough value
  for (i = 0; i < n_as; i++)
  for (j = 0; j < n_ts; j++)
  {
    hough[i][j] = hough[i][j] / (double)pixelcount[i][j];
    if (hough[i][j] > best_hough)
    {
      best_hough = hough[i][j];
      best_a     = trial_as[i];
    }
  }
  printf("  Found a = %f\n", best_a);

  ////////////////////
  // OUTPUT RESULTS //
  ////////////////////

  char hough_filename[110];
  char gpi_filename[110];
  char png_filename[110];
  char hough_png_filename[110];
  sprintf(hough_filename, "%s.hough.dat", dat_filename);
  sprintf(gpi_filename, "%s.gpi", dat_filename);
  sprintf(png_filename, "%s.png", dat_filename);
  sprintf(hough_png_filename, "%s.hough.png", dat_filename);

  // Save Hough to file
  f = fopen(hough_filename, "w");
  printf("Writing hough file...\n");
  for (j = 0; j < n_ts; j++)
  for (i = 0; i < n_as; i++) // reverse order to make plots look right
  {
    if (i == 0)
      fprintf(f,"\n");
    if (is_dB)
      hough[i][j] = 10.0*log10(hough[i][j]);
    if (!isnan(hough[i][j]))
      fprintf(f, "%e %f %e %d\n", trial_as[i], trial_ts[j], hough[i][j], pixelcount[i][j]);
  }
  fclose(f);

  // Write gnuplot script
  f = fopen(gpi_filename, "w");
  fprintf(f, "print \"Running gnuplot script...\"\n");
  fprintf(f, "set terminal pngcairo size 800,600\n\n");

  fprintf(f, "set out '%s'\n", hough_png_filename);
  fprintf(f, "set xlabel \"'a' parameter\"\n");
  if (is_dB)
    fprintf(f, "set ylabel \"Mean dB\"\n");
  else
    fprintf(f, "set ylabel \"Mean power\"\n");
  fprintf(f, "set cblabel \"Thickness of parabola (pixels)\"\n");
  fprintf(f, "set logscale x\n");
  fprintf(f, "plot './%s' u 1:3:2 w l palette notitle \n\n", hough_filename);

  fprintf(f, "set out '%s'\n", png_filename);
  fprintf(f, "unset logscale x\n");
  fprintf(f, "set dgrid3d\n");
  fprintf(f, "set samples 10000\n");
  fprintf(f, "unset xlabel\n");
  fprintf(f, "unset ylabel\n");
  fprintf(f, "set xrange [%lf:%lf]\n", min_x, max_x);
  fprintf(f, "set yrange [%lf:%lf]\n", min_y, max_y);
  fprintf(f, "set cbrange [%f:%f]\n", -60.0, -38.0);
  fprintf(f, "set xlabel \"Doppler frequency (mHz)\"\n");
  fprintf(f, "set ylabel \"Delay (us)\"\n");
  fprintf(f, "set cblabel \"dB\"\n");
  fprintf(f, "set parametric\n");
  fprintf(f, "set dummy t\n");
  fprintf(f, "set trange [-pi:pi]\n");
  fprintf(f, "myx(t) = (t+pi)/(2*pi) * %lf + %lf\n", max_x - min_x, min_x);
  fprintf(f, "myy(t) = (t+pi)/(2*pi) * %lf + %lf\n", max_y - min_y, min_y);
  fprintf(f, "plot './%s' using (($1-%d)*%f):(($2-%d)*%f):3 with image notitle, \\\n", dat_filename, x_orig, dx, y_orig, dy);
  fprintf(f, "     myx(t),%lf*(myx(t))**2 w l notitle lc rgb 'white', \\\n", best_a);
  fprintf(f, "     %lf*sin(t),%lf*cos(t) w l notitle lc rgb 'green', \\\n", mask_ox, mask_oy);
  fprintf(f, "     myx(t), %lf w l notitle lc rgb 'green', \\\n", mask_y);
  fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green', \\\n", -mask_x);
  fprintf(f, "     %lf, myy(t) w l notitle lc rgb 'green'\n\n", mask_x);
  fclose(f);
}
