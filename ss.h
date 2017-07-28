#ifndef SS_H
#define SS_H

#include <stdio.h>

// Secondary spectrum
struct sec_spect {
    double x_orig, y_orig;            // The x/y coordinates of the origin (in "index" units)
    double dx, dy;                    // The x/y axis resolution
    char   xunits[32], yunits[32];    // The names of the x,y-units
    int    is_dB;                     // The data are in logarithmic units?
    double cbmin, cbmax;              // The best dynamic range for viewing the data
    int    xsize, ysize;              // The number of pixels in each direction
    double **data;                    // 2D array: the secondary spectrum itself
};

double ss_xidx( struct sec_spect *ss, double x );
double ss_yidx( struct sec_spect *ss, double y );

double ss_xunits( struct sec_spect *ss, double xidx );
double ss_yunits( struct sec_spect *ss, double yidx );

void ss_read( FILE *f, struct sec_spect *ss, int filetype );
void ss_write( FILE *f, struct sec_spect *ss, int filetype );
void ss_crop( struct sec_spect *old_ss, struct sec_spect *new_ss,
        double xmin, double xmax, double ymin, double ymax );

void ss_malloc( struct sec_spect *ss, int xsize, int ysize );
void ss_free( struct sec_spect *ss );

void ss_write_gnuplot( FILE *f, struct sec_spect *ss,
       char *filename, int border );

#endif
