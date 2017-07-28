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
    if (logspace) // logarithmic interpolation
        return hg->amin*pow(hg->amax / hg->amin,
                            (double)idx/(double)(hg->size-1));
    else // linear interpolation
        return hg->amin + (double)idx*(hg->amax - hg->amin) /
                                      (double)(hg->size-1);
}


double hg_a_to_idx( struct hough *hg, double a )
/* Convert from curvature units to "index" units */
{
    if (logspace) // logarithmic interpolation
        return (double)(hg->size - 1) * log(a / hg->amin) /
                                        log(hg->amax / hg->amin);
    else // linear interpolation
        return (a - hg->amin) / ((hg->amax - hg->amin)/hg->size);
}



void hg_calc_transform( struct hough *hg )
/* Calculate the hough transform of parabolas in the associated secondary
 * spectrum.
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
    if (pxdist <= 0.0 && pydist <= 0.0)
    {
        fprintf( stderr, "error: need at least one valid \"parabola " );
        fprintf( stderr, "proximity\" specification:\n" );
        fprintf( stderr, "  pxdist = %lf\n", ss->pxdist );
        fprintf( stderr, "  pydist = %lf\n", ss->pydist );
        exit(EXIT_FAILURE);
    }

    // Loop over pixels, first in x-direction
    int xidx, yidx, aidx;
    double x, y, a;
    double xdist, ydist;
    for (xidx = 0; xidx < ss->xsize; xidx++)
    {
        x = ss_xunits( ss, xidx );

        // Ignore pixels closer than "xmask" to the y-axis
        if (fabs(x) < ss->xmask)
            continue;

        // Loop over pixels in y-direction
        for (yidx = 0; yidx < ss->ysize; yidx++)
        {
            y = ss_yunits( ss, yidx );

            // Ignore pixels closer than "ymask" to the x-axis
            if (fabs(y) < ss->ymask)
                continue;

            // Ignore pixels in the wrong quadrants
            if (!(hg->quadrant & HG_Q1) && (x > 0.0) && (y > 0.0))  continue;
            if (!(hg->quadrant & HG_Q2) && (x < 0.0) && (y > 0.0))  continue;
            if (!(hg->quadrant & HG_Q3) && (x < 0.0) && (y < 0.0))  continue;
            if (!(hg->quadrant & HG_Q4) && (x > 0.0) && (y < 0.0))  continue;

            // Ignore pixels too close to the origin
            // YET TO IMPLEMENT!

            // Test closeness to parabola
            for (aidx = 0; aidx < hg->size; aidx++)
            {
                a = hg_idx_to_a( hg, aidx );

                // Check the horizontal distance to the parabola
                if (hg->pxdist > 0.0)
                {
                    xdist = x - sqrt(y/a);
                    if (fabs(xdist) > pxdist)
                        continue;
                }

                // Check the vertical distance to the parabola
                if (hg->pydist > 0.0)
                {
                    ydist = y - a*x*x;
                    if (fabs(ydist) > pydist)
                        continue;

                }

                // We've survived all the checks, so include this pixel
                // in the Hough Transform
                // YET TO IMPLEMENT
            }

        }
    }
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
    fread( &hg->npixels,   sizeof(int),    hg->size, f );
    fread( &hg->transform, sizeof(double), hg->size, f );
}


void hg_write( FILE *f, struct hough *hg )
/* Dump contents of hough struct to file.
 * f = file handle
 * hg = pointer to hough struct to be written out
 */
{
    fwrite( &hg->quadrant,  sizeof(int), 1, f );
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
    fwrite( &hg->npixels,   sizeof(int),    hg->size, f );
    fwrite( &hg->transform, sizeof(double), hg->size, f );
}


void hg_write_gnuplot( FILE *f, struct hough *hg, char *filename )
{
}
