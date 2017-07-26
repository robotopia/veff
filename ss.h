#ifndef SS_H
#define SS_H

// Secondary spectrum
struct ss {
    double **data;
    int    x_orig;
    int    y_orig;
    double dx, dy;
    double min_x, min_y;
    double max_x, max_y;
    double mask_o, mask_x, mask_y;
    double mask_ox, mask_oy;
    int    max_t;
    int    is_dB;
};

#endif
