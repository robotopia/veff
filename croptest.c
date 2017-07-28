#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ss.h"

int main()
{
    FILE *f = fopen("ss.out", "r");

    struct sec_spect ss;
    ss_read( f, &ss, 1 );
    fclose(f);

    ss.x_orig = ss.y_orig = 257;
    ss.dx     = ss.dy     = 0.195;
    sprintf( ss.xunits, "mHz" );
    sprintf( ss.yunits, "μs" );

    ss.is_dB = 1;

    ss.cbmin = -52.0;
    ss.cbmax = -38.0;

    f = fopen("mwa.ss", "w");
    ss_write( f, &ss, 0 );
    fclose(f);

    f = fopen("mwa.ss.gpi", "w");
    ss_write_gnuplot(f, &ss, "mwa.ss", 0);
    fclose(f);

    ss_free( &ss );

    return 0;
}
