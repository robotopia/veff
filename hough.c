#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ss.h"
#include "hough.h"

void hg_malloc( struct hough *hg, int size )
/* Allocate memory for hough struct */
{
    hg->size = size;
    hg->npixels   =    (int *)malloc( size * sizeof(int)    );
    hg->transform = (double *)malloc( size * sizeof(double) );
}


void hg_free( struct hough *hg )
/* Free memory for hough struct */
{
    free( hg->npixels   );
    free( hg->transform );
}



double hg_idx_to_a( struct hough *hg, int idx )
/* Convert from "index" units to curvature units */
{
    if (hg->logspace) // logarithmic interpolation
        return hg->amin*pow(hg->amax / hg->amin,
                            (double)idx/(double)(hg->size-1));
    else // linear interpolation
        return hg->amin + (double)idx*(hg->amax - hg->amin) /
                                      (double)(hg->size-1);
}


double hg_a_to_idx( struct hough *hg, double a )
/* Convert from curvature units to "index" units */
{
    if (hg->logspace) // logarithmic interpolation
        return (double)(hg->size - 1) * log(a / hg->amin) /
                                        log(hg->amax / hg->amin);
    else // linear interpolation
        return (a - hg->amin) / ((hg->amax - hg->amin)/(hg->size - 1));
}



void hg_calc_transform( struct hough *hg )
/* Calculate the hough transform of parabolas in the associated secondary
 * spectrum.
 * This function assumes that memory has already been allocated in the
 * hough struct.
 */
{
    struct sec_spect *ss = hg->ss; // A shorthand

    // Check that we have a valid pointer to sec spect data
    if (ss == NULL)
    {
        fprintf( stderr, "error: cannot calculate Hough Transform without " );
        fprintf( stderr, "a secondary spectrum\n" );
        exit(EXIT_FAILURE);
    }

    // Check to make sure we have finite ranges for "proximity to parabola"
    if (hg->pxdist <= 0.0 && hg->pydist <= 0.0)
    {
        fprintf( stderr, "error: need at least one valid \"parabola " );
        fprintf( stderr, "proximity\" specification:\n" );
        fprintf( stderr, "  pxdist = %lf\n", hg->pxdist );
        fprintf( stderr, "  pydist = %lf\n", hg->pydist );
        exit(EXIT_FAILURE);
    }

    // Check that the arange is sensible
    if (hg->amin >= hg->amax)
    {
        fprintf( stderr, "error: bad curvature limits (%lf,%lf)\n",
                hg->amin, hg->amax );
        exit(EXIT_FAILURE);
    }

    // Reset values to zero
    int i;
    for (i = 0; i < hg->size; i++)
    {
        hg->npixels[i]   = 0;
        hg->transform[i] = 0.0;
    }

    // Loop over pixels, first in x-direction
    int xidx, yidx, aidx;
    double x, y;    // x coord; y coord
    double xn, yn;  // scaled distance from axes (for use with x/y0mask)
    for (xidx = 0; xidx < ss->xsize; xidx++)
    {
        x  = ss_xunits( ss, xidx );
        xn = x / hg->x0mask;

        // Ignore pixels closer than "xmask" to the y-axis
        if (fabs(x) < hg->xmask)
            continue;

        // Loop over pixels in y-direction
        for (yidx = 0; yidx < ss->ysize; yidx++)
        {
            y = ss_yunits( ss, yidx );

            // Ignore pixels closer than "ymask" to the x-axis
            if (fabs(y) < hg->ymask)
                continue;

            // Ignore pixels in the wrong quadrants
            if (!(hg->quadrant & HG_Q1) && (x > 0.0) && (y > 0.0))
                continue;
            if (!(hg->quadrant & HG_Q2) && (x < 0.0) && (y > 0.0))
                continue;
            if (!(hg->quadrant & HG_Q3) && (x < 0.0) && (y < 0.0))
                continue;
            if (!(hg->quadrant & HG_Q4) && (x > 0.0) && (y < 0.0))
                continue;

            // Ignore pixels too close to the origin
            yn = y / hg->y0mask;
            if (hypot(xn,yn) <= 1.0)
                continue;

            // Find all sufficiently close parabolas
            double xleft   = x - hg->pxdist;
            double xright  = x + hg->pxdist;
            double ybottom = y - hg->pydist;
            double ytop    = y + hg->pydist;
            double axmin = y / (xleft*xleft);
            double axmax = y / (xright*xright);
            double aymin = ybottom / (x*x);
            double aymax = ytop    / (x*x);
            // Choose the most limiting of the two mins and maxs
            double amin = (x/y < 0.0 ? (axmin > aymin ? axmin : aymin) : 
                                       (axmin < aymin ? axmin : aymin));
            double amax = (x/y < 0.0 ? (axmax < aymax ? axmax : aymax) :
                                       (axmax > aymax ? axmax : aymax));
            // Make sure they don't exceed the allowed limits
            int amin_idx = (int) hg_a_to_idx(hg, amin < hg->amin ? hg->amin : amin);
            int amax_idx = (int)(hg_a_to_idx(hg, amax > hg->amax ? hg->amax : amax)+1.0);
            for (aidx = amin_idx; aidx <= amax_idx; aidx++)
            {
                hg->transform[aidx] += ss->data[xidx][yidx];
                hg->npixels[aidx]++;
            }
        }
    }
}


double hg_best_a( struct hough *hg )
/* Returns the curvature "a" with the highest value in the transform */
{
    double h;
    double best_h = hg->transform[0] / (double)hg->npixels[0];
    int aidx, best_aidx = 0;
    for (aidx = 0; aidx < hg->size; aidx++)
    {
        h = hg->transform[aidx] / (double)hg->npixels[aidx];
        if (h > best_h)
        {
            best_aidx = aidx;
            best_h    = h;
        }
    }

    return hg_idx_to_a( hg, best_aidx );
}


void hg_read( FILE *f, struct hough *hg )
/* Read in hough transform data from file.
 * f = handle of file containing hough data
 * hg = pointer to hough struct to be written to
 */
{
    int size;
    fread( &hg->quadrant,  sizeof(int), 1, f );
    fread( &hg->amin,      sizeof(double), 1, f );
    fread( &hg->amax,      sizeof(double), 1, f );
    fread( &hg->x0mask,    sizeof(double), 1, f );
    fread( &hg->y0mask,    sizeof(double), 1, f );
    fread( &hg->xmask,     sizeof(double), 1, f );
    fread( &hg->ymask,     sizeof(double), 1, f );
    fread( &hg->pxdist,    sizeof(double), 1, f );
    fread( &hg->pydist,    sizeof(double), 1, f );
    fread( &hg->logspace,  sizeof(int),    1, f );
    fread( &size,          sizeof(int),    1, f );
    hg_malloc( hg, size );
    fread( hg->npixels,    sizeof(int),    hg->size, f );
    fread( hg->transform,  sizeof(double), hg->size, f );
}


void hg_write( FILE *f, struct hough *hg, int filetype )
/* Dump contents of hough struct to file.
 * f = file handle
 * hg = pointer to hough struct to be written out
 */
{
    int i; // Generic loop counter

    switch (filetype)
    {
        case HG_BINARY:
            fwrite( &hg->quadrant,  sizeof(int),    1, f );
            fwrite( &hg->amin,      sizeof(double), 1, f );
            fwrite( &hg->amax,      sizeof(double), 1, f );
            fwrite( &hg->x0mask,    sizeof(double), 1, f );
            fwrite( &hg->y0mask,    sizeof(double), 1, f );
            fwrite( &hg->xmask,     sizeof(double), 1, f );
            fwrite( &hg->ymask,     sizeof(double), 1, f );
            fwrite( &hg->pxdist,    sizeof(double), 1, f );
            fwrite( &hg->pydist,    sizeof(double), 1, f );
            fwrite( &hg->logspace,  sizeof(int),    1, f );
            fwrite( &hg->size,      sizeof(int),    1, f );
            fwrite( hg->npixels,    sizeof(int),    hg->size, f );
            fwrite( hg->transform,  sizeof(double), hg->size, f );
            break;
        case HG_ASCII:
            for (i = 0; i < hg->size; i++)
                fprintf( f, "%lf %lf %d\n",
                            hg_idx_to_a( hg, (double)i ),
                            hg->transform[i],
                            hg->npixels[i] );
            break;
        default:
            fprintf( stderr, "error: unrecognised hough output format\n" );
            exit(EXIT_FAILURE);
            break;
    }
}


void hg_write_gnuplot( FILE *f, struct hough *hg, char *filename )
/* Writes a gnuplot script for viewing the hough transform result */
{
    if (hg->logspace)
        fprintf( f, "set logscale x\n\n" );

    fprintf( f, "set xlabel \"Curvature, η\"\n" );
    if (hg->ss->is_dB)
        fprintf( f, "set ylabel \"Mean dB\"\n" );
    else
        fprintf( f, "set ylabel \"Mean amplitude\"\n" );
    fprintf( f, "set xrange [*:*] noextend\n\n" );

    fprintf( f, "plot '%s' ", filename );
    fprintf( f, "using 1:($2/$3) with lines " );
    fprintf( f, "title \"Hough Transform with parabola distance " );
    fprintf( f, "(%.2f %s, %.2f %s)\"", hg->pxdist, hg->ss->xunits,
                                        hg->pydist, hg->ss->yunits );
}


void hg_write_ssmarkup_gnuplot( FILE *f, struct hough *hg, char *filename,
        int border )
/* Writes a gnuplot script for viewing the secondary spectrum, including
 * mark-up lines showing the best fit parabola and the masks.
 */
{
    struct sec_spect *ss = hg->ss; // A short-hand
    double xmin = ss_xunits( ss, 0 );
    double xmax = ss_xunits( ss, ss->xsize-1 );
    double ymin = ss_yunits( ss, 0 );
    double ymax = ss_yunits( ss, ss->ysize-1 );
    double best_a = hg_best_a( hg );

    // First, write the secondary spectrum gnuplot script like normal
    ss_write_gnuplot( f, ss, filename, border );
    fprintf( f, ", \\\n" ); // (continue on the next line)

    // Now append the lines for the masks
    fprintf( f, "%lf*x**2 lc rgb 'green' notitle, \\\n", best_a );
    fprintf( f, "%lf*sqrt(1-(x/%lf)**2) lc rgb 'green' notitle, \\\n",
            hg->y0mask, hg->x0mask );
    fprintf( f, "-%lf*sqrt(1-(x/%lf)**2) lc rgb 'green' notitle, \\\n",
            hg->y0mask, hg->x0mask );
    fprintf( f, "'-' w l lc rgb 'green' notitle, \\\n" );
    fprintf( f, "'-' w l lc rgb 'green' notitle, \\\n" );
    fprintf( f, "'-' w l lc rgb 'green' notitle, \\\n" );
    fprintf( f, "'-' w l lc rgb 'green' notitle\n" );
    fprintf( f, "%lf %lf\n%lf %lf\ne\n", -hg->xmask, ymin, -hg->xmask, ymax );
    fprintf( f, "%lf %lf\n%lf %lf\ne\n",  hg->xmask, ymin,  hg->xmask, ymax );
    fprintf( f, "%lf %lf\n%lf %lf\ne\n", xmin, -hg->ymask, xmax, -hg->ymask );
    fprintf( f, "%lf %lf\n%lf %lf\ne\n", xmin,  hg->ymask, xmax,  hg->ymask );
}
