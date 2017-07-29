#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ss.h"

double ss_xidx( struct sec_spect *ss, double x )
/* Translates from physical units to "index" units
 */
{
    return x / ss->dx + ss->x_orig;
}

double ss_yidx( struct sec_spect *ss, double y )
/* Translates from physical units to "index" units
 */
{
    return y / ss->dy + ss->y_orig;
}

double ss_xunits( struct sec_spect *ss, double xidx )
/* Translates from "index" units to physical units
 */
{
    return (xidx - ss->x_orig) * ss->dx;
}

double ss_yunits( struct sec_spect *ss, double yidx )
/* Translates from "index" units to physical units
 */
{
    return (yidx - ss->y_orig) * ss->dy;
}


void ss_read( FILE *f, struct sec_spect *ss, int filetype )
/* Read in data from file stream and store it in a secondary spectrum
 *   f = file handle
 *   ss = pointer to struct to read data into
 *   filetype = code for data format:
 *       0 = binary dump of sec_spect struct
 *       1 = three-column ASCII [x, y, data] (x & y assumed 1-offset)
 *
 * Note: if filetype == 1, then no meta-information about the spectrum
 * is written to ss (origin, resolution, etc).
 *
 * This function allocates memory (in ss). Free with
 *   ss_free()
 */
{
    int xsize, ysize; // The amount of data to be read in
    int xidx, yidx;   // x and y "coordinates" for read-in data
    double val;       // A dummy variable for reading in unneeded numbers
    int nscan;        // The number of items read in by fscanf
    int i;            // A generic counter

    switch (filetype)
    {
        case SS_BINARY: // binary dump of sec_spect struct
            fread( &ss->x_orig, sizeof(double), 1, f );
            fread( &ss->y_orig, sizeof(double), 1, f );
            fread( &ss->dx,     sizeof(double), 1, f );
            fread( &ss->dy,     sizeof(double), 1, f );
            fread( ss->xunits,  sizeof(char),  32, f );
            fread( ss->yunits,  sizeof(char),  32, f );
            fread( &ss->is_dB,  sizeof(int),    1, f );
            fread( &ss->cbmin,  sizeof(double), 1, f );
            fread( &ss->cbmax,  sizeof(double), 1, f );
            fread( &xsize,      sizeof(int),    1, f );
            fread( &ysize,      sizeof(int),    1, f );
            ss_malloc( ss, xsize, ysize );
            for (i = 0; i < xsize; i++)
                fread( &ss->data[i], sizeof(double), ysize, f );
            break;
        case SS_ASCII: // three-column ASCII [x, y, data]
            // Assumes x and y start counting at 1
            // First, read the whole file to get the last three numbers
            // (= [xsize ysize last_value])
            while ((nscan = fscanf(f, "%d %d %lf", &xsize, &ysize, &val)) != EOF)
            // while loop will exit when exactly 0 items are read in
            {
                if (nscan < 3)
                {
                    fprintf( stderr, "error: error reading in secondary " );
                    fprintf( stderr, "spectrumfrom ASCII. Check file.\n" );
                    exit(EXIT_FAILURE);
                }
            }

            // Allocate memory
            ss_malloc( ss, xsize, ysize );

            // Now rewind and read them all in properly
            rewind( f );
            while (fscanf(f, "%d %d %lf", &xidx, &yidx, &val) != EOF)
                ss->data[xidx-1][yidx-1] = val;
            break;
        default:
            fprintf( stderr, "error: filetype must be either\n" );
            fprintf( stderr, "  0 = binary dump of sec_spect struct\n" );
            fprintf( stderr, "  1 = three-column ASCII (x, y, data)\n" );
            exit(EXIT_FAILURE);
            break;
    }
}


void ss_write( FILE *f, struct sec_spect *ss, int filetype )
/* Write secondary spectrum data to file stream
 *   f = file handle
 *   ss = pointer to struct to be written out
 *   filetype = code for data format:
 *       0 = binary dump of sec_spect struct
 *       1 = three-column ASCII [x, y, data] (x & y are 1-offset)
 *
 * Note: if filetype == 1, then no meta-information about the spectrum
 * is written (origin, resolution, etc).
 *
 */
{
    int i;     // A generic counter
    int x, y;  // x and y idx counters

    switch (filetype)
    {
        case SS_BINARY: // binary dump of sec_spect struct
            fwrite( &ss->x_orig, sizeof(double), 1, f );
            fwrite( &ss->y_orig, sizeof(double), 1, f );
            fwrite( &ss->dx,     sizeof(double), 1, f );
            fwrite( &ss->dy,     sizeof(double), 1, f );
            fwrite( ss->xunits,  sizeof(char),  32, f );
            fwrite( ss->yunits,  sizeof(char),  32, f );
            fwrite( &ss->is_dB,  sizeof(int),    1, f );
            fwrite( &ss->cbmin,  sizeof(double), 1, f );
            fwrite( &ss->cbmax,  sizeof(double), 1, f );
            fwrite( &ss->xsize,  sizeof(int),    1, f );
            fwrite( &ss->ysize,  sizeof(int),    1, f );
            for (i = 0; i < ss->xsize; i++)
                fwrite( ss->data[i], sizeof(double), ss->ysize, f );
            break;
        case SS_ASCII: // three-column ASCII [x, y, data]
            // Start x and y counting at 1
            for (x = 0; x < ss->xsize ; x++)
            for (y = 0; y < ss->ysize ; y++)
                fprintf(f, "%d %d %lf\n", x+1, y+1, ss->data[x][y]);
            break;
        default:
            fprintf( stderr, "error: filetype must be either\n" );
            fprintf( stderr, "  0 = binary dump of sec_spect struct\n" );
            fprintf( stderr, "  1 = three-column ASCII (x, y, data)\n" );
            exit(EXIT_FAILURE);
            break;
    }
}


void ss_crop( struct sec_spect *old_ss, struct sec_spect *new_ss,
        double xmin, double xmax, double ymin, double ymax )
/* Crops a secondary spectrum to the desired new limits.
 * INPUTS:
 *   old_ss    = pointer to secondary spectrum to be cropped
 *   xmin/xmax = new left/right border (in physical units)
 *   ymin/ymay = new top/bottom border (in physical units)
 * OUTPUTS:
 *   new_ss    = pointer to secondary spectrum with unalloc-
 *               ated "data"
 *
 * This function allocates memory (in new_ss). Free with
 *   ss_free()
 */
{
    int xmin_idx = floor(ss_xidx( old_ss, xmin ));
    int xmax_idx =  ceil(ss_xidx( old_ss, xmax ));
    int ymin_idx = floor(ss_yidx( old_ss, ymin ));
    int ymax_idx =  ceil(ss_yidx( old_ss, ymax ));

    // Make sure we're within limits
    if (xmin_idx <  0            )  xmin_idx = 0;
    if (ymin_idx <  0            )  ymin_idx = 0;
    if (xmax_idx >= old_ss->xsize)  xmax_idx = old_ss->xsize - 1;
    if (ymax_idx >= old_ss->ysize)  ymax_idx = old_ss->ysize - 1;

    // Calculate the new meta values
    new_ss->x_orig = old_ss->x_orig - xmin_idx;
    new_ss->y_orig = old_ss->y_orig - ymin_idx;
    new_ss->dx     = old_ss->dx;
    new_ss->dy     = old_ss->dy;
    strcpy( new_ss->xunits, old_ss->xunits );
    strcpy( new_ss->yunits, old_ss->yunits );
    new_ss->is_dB  = old_ss->is_dB;
    new_ss->cbmin  = old_ss->cbmin;
    new_ss->cbmax  = old_ss->cbmax;

    // Allocate memory
    ss_malloc( new_ss, xmax_idx - xmin_idx + 1,
                       ymax_idx - ymin_idx + 1 );

    // Copy across values
    int xidx, yidx;
    for (xidx = xmin_idx; xidx <= xmax_idx; xidx++)
    for (yidx = ymin_idx; yidx <= ymax_idx; yidx++)
        new_ss->data[xidx-xmin_idx][yidx-ymin_idx] = old_ss->data[xidx][yidx];
}

void ss_malloc( struct sec_spect *ss, int xsize, int ysize )
/* Allocates secondary spectrum memory according to the supplied sizes
 */
{
    ss->xsize = xsize;
    ss->ysize = ysize;
    ss->data = (double **)malloc( xsize * sizeof(double *) );
    int i;
    for (i = 0; i < xsize; i++)
        ss->data[i] = (double *)malloc( ysize * sizeof(double) );
}


void ss_free( struct sec_spect *ss )
{
    int i;
    for (i = 0; i < ss->xsize; i++)
        free( ss->data[i] );
    free( ss->data );
}


void ss_write_gnuplot( FILE *f, struct sec_spect *ss,
        char *filename, int border )
{
    // Get the limits of the plot
    double xmin = ss_xunits( ss, 0 );
    double xmax = ss_xunits( ss, ss->xsize-1 );
    double ymin = ss_yunits( ss, 0 );
    double ymax = ss_yunits( ss, ss->ysize-1 );

    int cbrange = 0;
    char cbmin[16], cbmax[16];
    if (isnan(ss->cbmin))
        sprintf( cbmin, "*" );
    else
    {
        cbrange = 1;
        fprintf( f, "cbmin = %lf\n", ss->cbmin );
        sprintf( cbmin, "cbmin" );
    }

    if (isnan(ss->cbmax))
        sprintf( cbmax, "*" );
    else
    {
        cbrange = 1;
        fprintf( f, "cbmax = %lf\n", ss->cbmax );
        sprintf( cbmax, "cbmax" );
    }

    if (cbrange)
        fprintf( f, "set cbrange [%s:%s] noextend\n\n", cbmin, cbmax );

    fprintf( f, "set xrange [%lf:%lf] noextend\n",   xmin, xmax );
    fprintf( f, "set yrange [%lf:%lf] noextend\n\n", ymin, ymax );

    if (border)
    {
        fprintf( f, "set xlabel \"Doppler frequency (%s)\"\n", ss->xunits);
        fprintf( f, "set ylabel \"Delay (%s)\"\n", ss->yunits);

        if (ss->is_dB)
            fprintf( f, "set cblabel \"dB\"\n\n");
        else
            fprintf( f, "set cblabel \"Amplitude\"\n\n");
    }
    else
    {
        fprintf( f, "set lmargin at screen 0\n" );
        fprintf( f, "set rmargin at screen 1\n" );
        fprintf( f, "set bmargin at screen 0\n" );
        fprintf( f, "set tmargin at screen 1\n" );
        fprintf( f, "unset colorbox\n\n" );
    }

    fprintf( f, "plot '%s' binary \\\n", filename );
    fprintf( f, "    skip=%d \\\n",  6*sizeof(double) +
                                    64*sizeof(char)   +
                                     3*sizeof(int)      );
    fprintf( f, "    array=(%d,%d) \\\n", ss->xsize, ss->ysize );
    fprintf( f, "    scan=yx \\\n" );
    fprintf( f, "    origin=(%lf,%lf) \\\n", xmin, ymin );
    fprintf( f, "    format=\"%%double\" \\\n" );
    fprintf( f, "    dx=%lf dy=%lf \\\n", ss->dx, ss->dy );
    fprintf( f, "    with image notitle\n" );
}
