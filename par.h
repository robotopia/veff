#ifndef PAR_H
#define PAR_H

#include <stdio.h>

#define AU  1.49597871e8 // km
#define KM2AU(x)    ((x)/AU)
#define KPC2KM(x)   ((x)*3.0857e16)
#define YRS2SEC(x)  ((x)*365.25*8.64e4)
#define MAS2RAD(x)  (DEG2RAD((x)/3.6e6))
#define PM2TV(kpc,masyr)  (KPC2KM(kpc)*MAS2RAD(masyr)/YRS2SEC(1.0))

FILE *open_par( char *filename );

#endif
