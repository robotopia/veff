#ifndef PAR_H
#define PAR_H

#include <stdio.h>

#define MAXSTRLEN 4096

#define AU  1.49597871e8 // km
#define KM2AU(x)    ((x)/AU)
#define KPC2KM(x)   ((x)*3.0857e16)
#define YRS2SEC(x)  ((x)*365.25*8.64e4)
#define MAS2RAD(x)  (DEG2RAD((x)/3.6e6))
#define PM2TV(kpc,masyr)  (KPC2KM(kpc)*MAS2RAD(masyr)/YRS2SEC(1.0))

struct pardata {
    char   *psrj;
    char   *name;
    char   *raj;
    double raj_err;
    char   *decj;
    double decj_err;
    double elong;
    double elat;
    double dm;
    double dm_err;
    double pepoch;
    double f0;
    double f0_err;
    double f1;
    double f1_err;
    double p0;
    double p0_err;
    double p1;
    double p1_err;
    double dist_dm;
    double dist_dm1;
    char   *survey;
    double pmra;
    double pmra_err;
    double pmdec;
    double pmdec_err;
    double s150;
    double s150_err;
    double s400;
    double s400_err;
    double s600;
    double s600_err;
    double s1400;
    double s1400_err;
    double s5000;
    double s5000_err;
    double w50;
    double w50_err;
    double w10;
    double w10_err;
    double posepoch;
    double dmepoch;
    char   *binary;
    double pb;
    double pb_err;
    double ecc;
    double ecc_err;
    double a1;
    double a1_err;
    double t0;
    double t0_err;
    double om;
    double om_err;
    double omdot;
    double omdot_err;
    double dist_a;
    double dist_a_err;
    char   *assoc;
    char   *bincomp;
    double age;
    double pbdot;
    double pbdot_err;
    double spindx;
    double spindx_err;
    double rm;
    double rm_err;
    double s700;
    double s700_err;
    double s3000;
    double s3000_err;
    double px;
    double px_err;
    char   *type;
    double s1600;
    double s1600_err;
    double fb0;
    double fb0_err;
    double fb1;
    double fb1_err;
    int    date;
    double dist;
    double dist1;
    double pmerr_pa;
    double p1_i;
    double age_i;
    double bsurf_i;
    double edot_i;
    double edotd2;
    double r_lum;
    double r_lum14;
    double pmtot;
    double pmtot_err;
    double vtrans;
    double bsurf;
    double b_lc;
    double si414;
    double edot;
    double rajd;
    double decjd;
    char   *osurvey;
    double massfn;
    double minmass;
    double medmass;
    double dmsinb;
    double gl;
    double gb;
    double xx;
    double yy;
    double zz;
    double pml;
    double pml_err;
    double pmb;
    double pmb_err;
    double kom;
    double kom_err;
    double kin;
    double kin_err;
    double m2;
    double m2_err;
    double uprmass;
    double minomdot;
    char   *units;
    double sini;
    double sini_err;
};

FILE *open_par( char *filename );
int get_par_double( FILE *f, char *param, double *val, double *err );
void read_par( char *filename, struct pardata *pd );
void free_par( struct pardata *pd );

#endif
