#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ss.h"
#include "../hough.h"

int main()
{
    // Set up the cropped secondary spectrum
    FILE *f = fopen("ss.out", "r");

    struct sec_spect ss;
    struct sec_spect ss_cropped;
    ss_read( f, &ss, 1 );

    fclose(f);

    ss.x_orig = ss.y_orig = 257;
    ss.dx     = ss.dy     = 0.195;
    sprintf( ss.xunits, "mHz" );
    sprintf( ss.yunits, "μs" );

    ss.is_dB = 1;

    ss.cbmin = -52.0;
    ss.cbmax = -38.0;

    ss_crop( &ss, &ss_cropped, -20.0, 20.0, 0.0, 20.0 );

    // Hough Transform
    struct hough hg;
    hg.ss        = &ss_cropped;
    hg.quadrant  = HG_Q2;
    hg.amin      = 0.1;
    hg.amax      = 2.0;
    hg.x0mask    = 4.0;
    hg.y0mask    = 4.0;
    hg.xmask     = 1.6;
    hg.ymask     = 1.4;
    hg.pxdist    = 0.2;
    hg.pydist    = 0.2;
    hg.logspace  = 1;
    hg_malloc( &hg, 200 );

    hg_calc_transform( &hg );

    // Write out results
    f = fopen("mwa_cropped.ss", "w");
    ss_write( f, &ss_cropped, SS_BINARY );
    fclose(f);

    f = fopen("mwa_cropped.ss.gpi", "w");
    ss_write_gnuplot( f, &ss_cropped, "mwa_cropped.ss", 1 );
    fclose(f);

    f = fopen("mwa_cropped_Q2.hg", "w");
    hg_write( f, &hg, HG_BINARY );
    fclose(f);

    f = fopen("mwa_cropped_Q2.hg.txt", "w");
    hg_write( f, &hg, HG_ASCII );
    fclose(f);

    f = fopen("mwa_cropped_Q2.hg.gpi", "w");
    hg_write_gnuplot( f, &hg, "mwa_cropped_Q2.hg" );
    fclose(f);

    ss_free( &ss );
    ss_free( &ss_cropped );
    hg_free( &hg );

    return 0;
}
