#ifndef SS_H
#define SS_H

// Secondary spectrum
struct sec_spect {
    double **data;                    // 2D array: the secondary spectrum itself
    char   dat_filename[100];         // The name of the file whence the data came
    int    filetype;                  // For future support of different data formats
    int    x_orig, y_orig;            // The x/y-index of the origin
    double dx, dy;                    // The x/y axis resolution
    double min_x, min_y;              // The minimum x/y index to be used in
                                      // the parabola-finding algorithm
    double max_x, max_y;              // The maximum x/y index to be used in
                                      // the parabola-finding algorithm
    double mask_o, mask_x, mask_y;    // If a pixel is this many pixels away from the
                                      // origin/x-axis/y-axis, it will be ignored in
                                      // the parabola-finding algorithm
    double mask_ox, mask_oy;          // The x/y components of "mask_o"
    int    max_t;                     // The maximum distance away (in pixels) from
                                      // the parabola still considered part of the
                                      // parabola
    int    is_dB;                     // The data are in logarithmic units?
    double cbmin, cbmax;              // The best dynamic range for viewing the data
};

#endif
