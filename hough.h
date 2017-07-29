#ifndef HOUGH_H
#define HOUGH_H

#include "ss.h"

// Quadrant flags
#define HG_Q1  1<<0
#define HG_Q2  1<<1
#define HG_Q3  1<<2
#define HG_Q4  1<<3

#define HG_BINARY  0
#define HG_ASCII   1

struct hough {
    struct sec_spect *ss;      // A pointer to the secondary spectrum on which the Hough Transform is to be operated
    int    quadrant;           // Which quadrant(s) to include: HG_Q1 | HG_Q2 | HG_Q3 | HG_Q4
                               //    Q2 | Q1
                               //    ---+---
                               //    Q3 | Q4
    double amin;               // First parabolic curvature to try
    double amax;               // Last  parabolic curvature to try
    double x0mask;             // Ignore pixels within ellipse defined by x0mask and y0mask ...
    double y0mask;             //   centered on origin
    double xmask;              // Ignore pixels within this range of y-axis
    double ymask;              // Ignore pixels within this range of x-axis
    double pxdist;             // Include pixels within this (x-)range of parabola (neg = ignore)
    double pydist;             // Include pixels within this (y-)range of parabola (neg = ignore)
    int    logspace;           // (boolean) Space out points [amin:amax] logarithmically
    int    size;               // Number of points in the Hough Transform
    int    *npixels;           // The number of ss pixels that went into each hough pixel
    double *transform;         // The result of the Hough Transform
};

void hg_malloc( struct hough *hg, int size );
void hg_free( struct hough *hg );

double hg_idx_to_a( struct hough *hg, int idx );
double hg_a_to_idx( struct hough *hg, double a );

void hg_calc_transform( struct hough *hg );

void hg_read( FILE *f, struct hough *hg );
void hg_write( FILE *f, struct hough *hg, int filetype );

void hg_write_gnuplot( FILE *f, struct hough *hg, char *filename );

#endif
