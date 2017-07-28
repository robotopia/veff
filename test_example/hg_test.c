#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../ss.h"
#include "../hough.h"

int main()
{
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

    f = fopen("mwa_cropped.ss", "w");
    ss_write( f, &ss_cropped, SS_BINARY );
    fclose(f);

    f = fopen("mwa_cropped.ss.gpi", "w");
    ss_write_gnuplot( f, &ss_cropped, "mwa_cropped.ss", 1 );
    fclose(f);

    ss_free( &ss );
    ss_free( &ss_cropped );

    return 0;
}
