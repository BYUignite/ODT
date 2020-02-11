/* ckwyp.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
*/

#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Table of constant values */

static double c_b8 = 10.;

#include <cmath>

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/*     Automatically generated code */

/*     Tianfeng Lu */
/*     University of Connecticut */
/*     191 Auditorium Road U-3139 */
/*     Storrs, CT 06269, USA */
/*     Email: tlu@engr.uconn.edu */

/*     August 03, 2014 */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
// /* Subroutine */ int ckwyp_(double *p, double *t, double *y,
// 	long int *ickwrk, double *rckwrk, double *wdot)
// {



void getProblemSpecificRR(double rho, double t, double p,
        double *y, double *wdot) {

    // rho input not used but fits the generic interface: units kg/m3
    // t input units are K
    // p input units are Pa
    // y input mass fractions
    // wdot is output: units are kmol/m^3*s


    p *= 10.0;      // convert Pa to dynes/cm^2



    static double c__[86], rb[625], rf[625], xq[13];
    extern /* Subroutine */ int qssa_(double *, double *, double *
            ), rdot_(double *, double *, double *), ratt_(
                double *, double *, double *, double *), ratx_(
                    double *, double *, double *, double *,
                    double *), ytcp_(double *, double *, double *,
                        double *);
    static double rklow[43];


    /* Parameter adjustments */
    --wdot;
    //--rckwrk;
    //--ickwrk;
    --y;

    /* Function Body */
    ytcp_(&p, &t, &y[1], c__);
    ratt_(&t, rf, rb, rklow);
    ratx_(&t, c__, rf, rb, rklow);
    qssa_(rf, rb, xq);
    rdot_(rf, rb, &wdot[1]);

    // convert units from mol/cm3*s to kmol/m3*s

    for(int k=1; k<=86; k++)
        wdot[k] *= 1000.0;

    //    return 0;
} /* ckwyp_ */

///*                                                                      C */
///* ----------------------------------------------------------------------C */
///*                                                                      C */
///* Subroutine */ int ckwyr_(double *rho, double *t, double *y,
//	long int *ickwrk, double *rckwrk, double *wdot)
//{
//    static double c__[86], rb[625], rf[625], xq[13];
//    extern /* Subroutine */ int qssa_(double *, double *, double *
//	    ), rdot_(double *, double *, double *), ratt_(
//	    double *, double *, double *, double *), ratx_(
//	    double *, double *, double *, double *,
//	    double *), ytcr_(double *, double *, double *,
//	    double *);
//    static double rklow[43];
//
//
//    /* Parameter adjustments */
//    --wdot;
//    --rckwrk;
//    --ickwrk;
//    --y;
//
//    /* Function Body */
//    ytcr_(rho, t, &y[1], c__);
//    ratt_(t, rf, rb, rklow);
//    ratx_(t, c__, rf, rb, rklow);
//    qssa_(rf, rb, xq);
//    rdot_(rf, rb, &wdot[1]);
//
//    return 0;
//} /* ckwyr_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int ytcp_(double *p, double *t, double *y,
        double *c__)
{
    static long int k;
    static double sum;


    /* Parameter adjustments */
    --c__;
    --y;

    /* Function Body */
    c__[1] = y[1] / 1.007969975471497;
    c__[2] = y[2] / 2.015939950942993;
    c__[3] = y[3] / 15.99940013885498;
    c__[4] = y[4] / 31.99880027770996;
    c__[5] = y[5] / 17.00737011432648;
    c__[6] = y[6] / 18.01534008979797;
    c__[7] = y[7] / 33.00677025318146;
    c__[8] = y[8] / 34.01474022865295;
    c__[9] = y[9] / 28.0105504989624;
    c__[10] = y[10] / 44.00995063781738;
    c__[11] = y[11] / 30.0264904499054;
    c__[12] = y[12] / 29.0185204744339;
    c__[13] = y[13] / 32.04243040084839;
    c__[14] = y[14] / 16.04303026199341;
    c__[15] = y[15] / 15.03506028652191;
    c__[16] = y[16] / 14.02709031105042;
    c__[17] = y[17] / 30.07012057304382;
    c__[18] = y[18] / 29.06215059757233;
    c__[19] = y[19] / 28.05418062210083;
    c__[20] = y[20] / 27.04621064662933;
    c__[21] = y[21] / 26.03824067115784;
    c__[22] = y[22] / 44.05358076095581;
    c__[23] = y[23] / 44.05358076095581;
    c__[24] = y[24] / 43.04561078548431;
    c__[25] = y[25] / 42.03764081001282;
    c__[26] = y[26] / 41.02967083454132;
    c__[27] = y[27] / 44.05358076095581;
    c__[28] = y[28] / 56.06473112106323;
    c__[29] = y[29] / 44.09721088409424;
    c__[30] = y[30] / 42.08127093315125;
    c__[31] = y[31] / 41.07330095767975;
    c__[32] = y[32] / 41.07330095767975;
    c__[33] = y[33] / 41.07330095767975;
    c__[34] = y[34] / 40.06533098220825;
    c__[35] = y[35] / 39.05736100673676;
    c__[36] = y[36] / 54.09242129325867;
    c__[37] = y[37] / 52.07648134231567;
    c__[38] = y[38] / 51.06851136684418;
    c__[39] = y[39] / 51.06851136684418;
    c__[40] = y[40] / 50.06054139137268;
    c__[41] = y[41] / 53.08445131778717;
    c__[42] = y[42] / 53.08445131778717;
    c__[43] = y[43] / 54.09242129325867;
    c__[44] = y[44] / 66.05994153022766;
    c__[45] = y[45] / 78.11472201347351;
    c__[46] = y[46] / 66.10357165336609;
    c__[47] = y[47] / 77.10675203800201;
    c__[48] = y[48] / 93.106152176857;
    c__[49] = y[49] / 65.09560167789459;
    c__[50] = y[50] / 80.08703184127808;
    c__[51] = y[51] / 91.13384234905243;
    c__[52] = y[52] / 92.14181232452393;
    c__[53] = y[53] / 102.1370227336884;
    c__[54] = y[54] / 106.1252725124359;
    c__[55] = y[55] / 101.1290527582169;
    c__[56] = y[56] / 104.1529626846314;
    c__[57] = y[57] / 127.1672934293747;
    c__[58] = y[58] / 128.1752634048462;
    c__[59] = y[59] / 152.197564125061;
    c__[60] = y[60] / 177.2278348207474;
    c__[61] = y[61] / 151.1895941495895;
    c__[62] = y[62] / 202.2581055164337;
    c__[63] = y[63] / 201.2501355409622;
    c__[64] = y[64] / 201.2501355409622;
    c__[65] = y[65] / 201.2501355409622;
    c__[66] = y[66] / 106.1689026355743;
    c__[67] = y[67] / 115.1561430692673;
    c__[68] = y[68] / 142.2023537158966;
    c__[69] = y[69] / 130.1475732326508;
    c__[70] = y[70] / 116.1641130447388;
    c__[71] = y[71] / 226.2804062366486;
    c__[72] = y[72] / 252.3186469078064;
    c__[73] = y[73] / 252.3186469078064;
    c__[74] = y[74] / 226.2804062366486;
    c__[75] = y[75] / 226.2804062366486;
    c__[76] = y[76] / 226.2804062366486;
    c__[77] = y[77] / 225.2724362611771;
    c__[78] = y[78] / 225.2724362611771;
    c__[79] = y[79] / 251.3106769323349;
    c__[80] = y[80] / 225.2724362611771;
    c__[81] = y[81] / 251.3106769323349;
    c__[82] = y[82] / 276.3409476280212;
    c__[83] = y[83] / 276.3409476280212;
    c__[84] = y[84] / 275.3329776525497;
    c__[85] = y[85] / 300.3632483482361;
    c__[86] = y[86] / 28.01339912414551;

    sum = 0.;
    for (k = 1; k <= 86; ++k) {
        sum += c__[k];
    }
    sum = *p / (sum * *t * 83145100.);

    for (k = 1; k <= 86; ++k) {
        c__[k] *= sum;
    }

    return 0;
} /* ytcp_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int ytcr_(double *rho, double *t, double *y,
        double *c__)
{

    static long int k;

    /* Parameter adjustments */
    --c__;
    --y;

    /* Function Body */

    /*     H */
    c__[1] = y[1] / 1.007969975471497;
    /*     H2 */
    c__[2] = y[2] / 2.015939950942993;
    /*     O */
    c__[3] = y[3] / 15.99940013885498;
    /*     O2 */
    c__[4] = y[4] / 31.99880027770996;
    /*     OH */
    c__[5] = y[5] / 17.00737011432648;
    /*     H2O */
    c__[6] = y[6] / 18.01534008979797;
    /*     HO2 */
    c__[7] = y[7] / 33.00677025318146;
    /*     H2O2 */
    c__[8] = y[8] / 34.01474022865295;
    /*     CO */
    c__[9] = y[9] / 28.0105504989624;
    /*     CO2 */
    c__[10] = y[10] / 44.00995063781738;
    /*     CH2O */
    c__[11] = y[11] / 30.0264904499054;
    /*     HCO */
    c__[12] = y[12] / 29.0185204744339;
    /*     CH3OH */
    c__[13] = y[13] / 32.04243040084839;
    /*     CH4 */
    c__[14] = y[14] / 16.04303026199341;
    /*     CH3 */
    c__[15] = y[15] / 15.03506028652191;
    /*     CH2 */
    c__[16] = y[16] / 14.02709031105042;
    /*     C2H6 */
    c__[17] = y[17] / 30.07012057304382;
    /*     C2H5 */
    c__[18] = y[18] / 29.06215059757233;
    /*     C2H4 */
    c__[19] = y[19] / 28.05418062210083;
    /*     C2H3 */
    c__[20] = y[20] / 27.04621064662933;
    /*     C2H2 */
    c__[21] = y[21] / 26.03824067115784;
    /*     CH3CHO */
    c__[22] = y[22] / 44.05358076095581;
    /*     C2H3OH */
    c__[23] = y[23] / 44.05358076095581;
    /*     CH2CHO */
    c__[24] = y[24] / 43.04561078548431;
    /*     CH2CO */
    c__[25] = y[25] / 42.03764081001282;
    /*     HCCO */
    c__[26] = y[26] / 41.02967083454132;
    /*     C2H4O1-2 */
    c__[27] = y[27] / 44.05358076095581;
    /*     C2H3CHO */
    c__[28] = y[28] / 56.06473112106323;
    /*     C3H8 */
    c__[29] = y[29] / 44.09721088409424;
    /*     C3H6 */
    c__[30] = y[30] / 42.08127093315125;
    /*     C3H5-A */
    c__[31] = y[31] / 41.07330095767975;
    /*     C3H5-S */
    c__[32] = y[32] / 41.07330095767975;
    /*     C3H5-T */
    c__[33] = y[33] / 41.07330095767975;
    /*     C3H4-P */
    c__[34] = y[34] / 40.06533098220825;
    /*     C3H3 */
    c__[35] = y[35] / 39.05736100673676;
    /*     C4H6 */
    c__[36] = y[36] / 54.09242129325867;
    /*     C4H4 */
    c__[37] = y[37] / 52.07648134231567;
    /*     C4H3-I */
    c__[38] = y[38] / 51.06851136684418;
    /*     C4H3-N */
    c__[39] = y[39] / 51.06851136684418;
    /*     C4H2 */
    c__[40] = y[40] / 50.06054139137268;
    /*     C4H5-I */
    c__[41] = y[41] / 53.08445131778717;
    /*     C4H5-2 */
    c__[42] = y[42] / 53.08445131778717;
    /*     C4H6-2 */
    c__[43] = y[43] / 54.09242129325867;
    /*     H2C4O */
    c__[44] = y[44] / 66.05994153022766;
    /*     A1 */
    c__[45] = y[45] / 78.11472201347351;
    /*     c-C5H6 */
    c__[46] = y[46] / 66.10357165336609;
    /*     A1- */
    c__[47] = y[47] / 77.10675203800201;
    /*     C6H5O */
    c__[48] = y[48] / 93.106152176857;
    /*     c-C5H5 */
    c__[49] = y[49] / 65.09560167789459;
    /*     C5H4O */
    c__[50] = y[50] / 80.08703184127808;
    /*     C6H5CH2 */
    c__[51] = y[51] / 91.13384234905243;
    /*     C6H5CH3 */
    c__[52] = y[52] / 92.14181232452393;
    /*     A1C2H */
    c__[53] = y[53] / 102.1370227336884;
    /*     A1CHO */
    c__[54] = y[54] / 106.1252725124359;
    /*     A1C2H* */
    c__[55] = y[55] / 101.1290527582169;
    /*     A1C2H3 */
    c__[56] = y[56] / 104.1529626846314;
    /*     A2-1 */
    c__[57] = y[57] / 127.1672934293747;
    /*     A2 */
    c__[58] = y[58] / 128.1752634048462;
    /*     A2R5 */
    c__[59] = y[59] / 152.197564125061;
    /*     A3-4 */
    c__[60] = y[60] / 177.2278348207474;
    /*     A2R5- */
    c__[61] = y[61] / 151.1895941495895;
    /*     A4 */
    c__[62] = y[62] / 202.2581055164337;
    /*     A4-1 */
    c__[63] = y[63] / 201.2501355409622;
    /*     A4-2 */
    c__[64] = y[64] / 201.2501355409622;
    /*     A4-4 */
    c__[65] = y[65] / 201.2501355409622;
    /*     A1C2H5 */
    c__[66] = y[66] / 106.1689026355743;
    /*     C9H7 */
    c__[67] = y[67] / 115.1561430692673;
    /*     A2CH3 */
    c__[68] = y[68] / 142.2023537158966;
    /*     C9H6O */
    c__[69] = y[69] / 130.1475732326508;
    /*     C9H8 */
    c__[70] = y[70] / 116.1641130447388;
    /*     A4R5 */
    c__[71] = y[71] / 226.2804062366486;
    /*     BAPYR */
    c__[72] = y[72] / 252.3186469078064;
    /*     BEPYREN */
    c__[73] = y[73] / 252.3186469078064;
    /*     PYC2H-1 */
    c__[74] = y[74] / 226.2804062366486;
    /*     PYC2H-2 */
    c__[75] = y[75] / 226.2804062366486;
    /*     PYC2H-4 */
    c__[76] = y[76] / 226.2804062366486;
    /*     PYC2H-1JP */
    c__[77] = y[77] / 225.2724362611771;
    /*     PYC2H-2JS */
    c__[78] = y[78] / 225.2724362611771;
    /*     BAPYRJS */
    c__[79] = y[79] / 251.3106769323349;
    /*     PYC2H-4JS */
    c__[80] = y[80] / 225.2724362611771;
    /*     BEPYRENJS */
    c__[81] = y[81] / 251.3106769323349;
    /*     BGHIPER */
    c__[82] = y[82] / 276.3409476280212;
    /*     ANTHAN */
    c__[83] = y[83] / 276.3409476280212;
    /*     BGHIPEJS1 */
    c__[84] = y[84] / 275.3329776525497;
    /*     CORONEN */
    c__[85] = y[85] / 300.3632483482361;
    /*     N2 */
    c__[86] = y[86] / 28.01339912414551;

    for (k = 1; k <= 86; ++k) {
        c__[k] = *rho * c__[k];
    }

    return 0;
} /* ytcr_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int ratt_(double *t, double *rf, double *rb,
        double *rklow)
{
    /* Builtin functions */
    double log(double), exp(double);

    /* Local variables */
    static double eg[97], ti, ti2, eqk, smh[97], pfac1, pfac2, pfac3,
                  alogt;
    extern /* Subroutine */ int rdsmh_(double *, double *);


    /* Parameter adjustments */
    --rklow;
    --rb;
    --rf;

    /* Function Body */
    alogt = log(*t);
    ti = 1. / *t;
    ti2 = ti * ti;

    rdsmh_(t, smh);
    eg[0] = exp(smh[0]);
    eg[1] = exp(smh[1]);
    eg[2] = exp(smh[2]);
    eg[3] = exp(smh[3]);
    eg[4] = exp(smh[4]);
    eg[5] = exp(smh[5]);
    eg[6] = exp(smh[6]);
    eg[7] = exp(smh[7]);
    eg[8] = exp(smh[8]);
    eg[9] = exp(smh[9]);
    eg[10] = exp(smh[10]);
    eg[11] = exp(smh[11]);
    eg[12] = exp(smh[12]);
    eg[13] = exp(smh[13]);
    eg[14] = exp(smh[14]);
    eg[15] = exp(smh[15]);
    eg[16] = exp(smh[16]);
    eg[17] = exp(smh[17]);
    eg[18] = exp(smh[18]);
    eg[19] = exp(smh[19]);
    eg[20] = exp(smh[20]);
    eg[21] = exp(smh[21]);
    eg[22] = exp(smh[22]);
    eg[23] = exp(smh[23]);
    eg[24] = exp(smh[24]);
    eg[25] = exp(smh[25]);
    eg[26] = exp(smh[26]);
    eg[27] = exp(smh[27]);
    eg[28] = exp(smh[28]);
    eg[29] = exp(smh[29]);
    eg[30] = exp(smh[30]);
    eg[31] = exp(smh[31]);
    eg[32] = exp(smh[32]);
    eg[33] = exp(smh[33]);
    eg[34] = exp(smh[34]);
    eg[35] = exp(smh[35]);
    eg[36] = exp(smh[36]);
    eg[37] = exp(smh[37]);
    eg[38] = exp(smh[38]);
    eg[39] = exp(smh[39]);
    eg[40] = exp(smh[40]);
    eg[41] = exp(smh[41]);
    eg[42] = exp(smh[42]);
    eg[43] = exp(smh[43]);
    eg[44] = exp(smh[44]);
    eg[45] = exp(smh[45]);
    eg[46] = exp(smh[46]);
    eg[47] = exp(smh[47]);
    eg[48] = exp(smh[48]);
    eg[49] = exp(smh[49]);
    eg[50] = exp(smh[50]);
    eg[51] = exp(smh[51]);
    eg[52] = exp(smh[52]);
    eg[53] = exp(smh[53]);
    eg[54] = exp(smh[54]);
    eg[55] = exp(smh[55]);
    eg[56] = exp(smh[56]);
    eg[57] = exp(smh[57]);
    eg[58] = exp(smh[58]);
    eg[59] = exp(smh[59]);
    eg[60] = exp(smh[60]);
    eg[61] = exp(smh[61]);
    eg[62] = exp(smh[62]);
    eg[63] = exp(smh[63]);
    eg[64] = exp(smh[64]);
    eg[66] = exp(smh[66]);
    eg[68] = exp(smh[68]);
    eg[69] = exp(smh[69]);
    eg[71] = exp(smh[71]);
    eg[72] = exp(smh[72]);
    eg[73] = exp(smh[73]);
    eg[74] = exp(smh[74]);
    eg[75] = exp(smh[75]);
    eg[76] = exp(smh[76]);
    eg[77] = exp(smh[77]);
    eg[78] = exp(smh[78]);
    eg[79] = exp(smh[79]);
    eg[80] = exp(smh[80]);
    eg[81] = exp(smh[81]);
    eg[82] = exp(smh[82]);
    eg[83] = exp(smh[83]);
    eg[84] = exp(smh[84]);
    eg[85] = exp(smh[85]);
    eg[86] = exp(smh[86]);
    eg[87] = exp(smh[87]);
    eg[88] = exp(smh[88]);
    eg[89] = exp(smh[89]);
    eg[90] = exp(smh[90]);
    eg[91] = exp(smh[91]);
    eg[92] = exp(smh[92]);
    eg[93] = exp(smh[93]);
    eg[94] = exp(smh[94]);
    eg[96] = exp(smh[96]);
    pfac1 = 1013250. / (*t * 83145100.);
    pfac2 = pfac1 * pfac1;
    pfac3 = pfac2 * pfac1;


    /*     R1: H + O2 = O + OH */
    rf[1] = exp(32.27541201506992 - ti * 7692.169981543812);
    eqk = eg[2] * eg[4] / eg[0] / eg[3];
    rb[1] = rf[1] / max(eqk,1e-200);
    /*     R2: H2 + O = H + OH */
    rf[2] = exp(alogt * 2.67 + 10.83565163356657 - ti * 3166.239272790375);
    eqk = eg[0] * eg[4] / eg[1] / eg[2];
    rb[2] = rf[2] / max(eqk,1e-200);
    /*     R3: H2 + OH = H + H2O */
    rf[3] = exp(31.41065493331095 - ti * 3517.484506803039);
    eqk = eg[0] * eg[5] / eg[1] / eg[4];
    rb[3] = rf[3] / max(eqk,1e-200);
    /*     R4: O + H2O = 2OH */
    rf[4] = exp(alogt * 2.02 + 14.90407251077888 - ti * 6743.103346374924);
    eqk = eg[4] * eg[4] / eg[2] / eg[5];
    rb[4] = rf[4] / max(eqk,1e-200);
    /*     R5: H2 = 2H */
    rf[5] = exp(45.27016052855837 - alogt * 1.4 - ti * 52535.82010160761);
    eqk = eg[0] * eg[0] / eg[1] * pfac1;
    rb[5] = rf[5] / max(eqk,1e-200);
    /*     R6: 2O = O2 */
    rf[6] = exp(36.35766453152699 - alogt * .5);
    eqk = eg[3] / eg[2] / eg[2] / pfac1;
    rb[6] = rf[6] / max(eqk,1e-200);
    /*     R7: H + O = OH */
    rf[7] = exp(42.99706847840676 - alogt * 1.);
    eqk = eg[4] / eg[0] / eg[2] / pfac1;
    rb[7] = rf[7] / max(eqk,1e-200);
    /*     R8: H + OH = H2O */
    rf[8] = exp(51.90963501436438 - alogt * 2.);
    eqk = eg[5] / eg[0] / eg[4] / pfac1;
    rb[8] = rf[8] / max(eqk,1e-200);
    /*     R9: H + O2 = HO2 */
    rf[9] = exp(alogt * .44 + 29.16788833552781);
    eqk = eg[6] / eg[0] / eg[3] / pfac1;
    rb[9] = rf[9] / max(eqk,1e-200);
    /*     R10: H + HO2 = 2OH */
    rf[10] = exp(31.89073886371465 - ti * 148.4489169537763);
    eqk = eg[4] * eg[4] / eg[0] / eg[6];
    rb[10] = rf[10] / max(eqk,1e-200);
    /*     R11: H2 + O2 = H + HO2 */
    rf[11] = exp(alogt * 2.433 + 13.15695802216883 - ti * 26923.09815207098);
    eqk = eg[0] * eg[6] / eg[1] / eg[3];
    rb[11] = rf[11] / max(eqk,1e-200);
    /*     R12: O + HO2 = O2 + OH */
    rf[12] = 3.25e13;
    eqk = eg[3] * eg[4] / eg[2] / eg[6];
    rb[12] = rf[12] / max(eqk,1e-200);
    /*     R13: OH + HO2 = O2 + H2O */
    rf[13] = exp(ti * 250.0986838170401 + 30.83214021920749);
    eqk = eg[3] * eg[5] / eg[4] / eg[6];
    rb[13] = rf[13] / max(eqk,1e-200);
    /*     R14: 2HO2 = O2 + H2O2 */
    rf[14] = exp(ti * 820.243168253069 + 25.59080028740199);
    eqk = eg[3] * eg[7] / eg[6] / eg[6];
    rb[14] = rf[14] / max(eqk,1e-200);
    /*     R15: 2HO2 = O2 + H2O2 */
    rf[15] = exp(33.53310785188531 - ti * 6038.600011679037);
    eqk = eg[3] * eg[7] / eg[6] / eg[6];
    rb[15] = rf[15] / max(eqk,1e-200);
    /*     R16: H2O2 = 2OH */
    rf[16] = exp(alogt * .9 + 28.32416829648849 - ti * 24531.30933077844);
    eqk = eg[4] * eg[4] / eg[7] * pfac1;
    rb[16] = rf[16] / max(eqk,1e-200);
    /*     R17: H2O2 = 2OH */
    rf[17] = exp(alogt * .9 + 28.32416829648849 - ti * 24531.30933077844);
    eqk = eg[4] * eg[4] / eg[7] * pfac1;
    rb[17] = rf[17] / max(eqk,1e-200);
    /*     R18: H + H2O2 = OH + H2O */
    rf[18] = exp(30.81323295642516 - ti * 1997.770170530481);
    eqk = eg[4] * eg[5] / eg[0] / eg[7];
    rb[18] = rf[18] / max(eqk,1e-200);
    /*     R19: H + H2O2 = H2 + HO2 */
    rf[19] = exp(alogt * 1. + 23.79131877208003 - ti * 3019.300005839518);
    eqk = eg[1] * eg[6] / eg[0] / eg[7];
    rb[19] = rf[19] / max(eqk,1e-200);
    /*     R20: O + H2O2 = OH + HO2 */
    rf[20] = exp(alogt * 2. + 16.07205171245691 - ti * 1997.770170530481);
    eqk = eg[4] * eg[6] / eg[2] / eg[7];
    rb[20] = rf[20] / max(eqk,1e-200);
    /*     R21: OH + H2O2 = H2O + HO2 */
    rf[21] = exp(28.18490622915499 - ti * 160.0229003094945);
    eqk = eg[5] * eg[6] / eg[4] / eg[7];
    rb[21] = rf[21] / max(eqk,1e-200);
    /*     R22: OH + H2O2 = H2O + HO2 */
    rf[22] = exp(31.96043780033013 - ti * 3657.881957074576);
    eqk = eg[5] * eg[6] / eg[4] / eg[7];
    rb[22] = rf[22] / max(eqk,1e-200);
    /*     R23: O + CO = CO2 */
    rf[23] = exp(23.33480513766778 - ti * 1199.668535653569);
    eqk = eg[9] / eg[2] / eg[8] / pfac1;
    rb[23] = rf[23] / max(eqk,1e-200);
    /*     R24: O2 + CO = O + CO2 */
    rf[24] = exp(27.74345654525834 - ti * 24003.43504642417);
    eqk = eg[2] * eg[9] / eg[3] / eg[8];
    rb[24] = rf[24] / max(eqk,1e-200);
    /*     R25: OH + CO = H + CO2 */
    rf[25] = exp(alogt * 2.053 + 11.15839108553061 + ti * 178.9941686795194);
    eqk = eg[0] * eg[9] / eg[4] / eg[8];
    rb[25] = rf[25] / max(eqk,1e-200);
    /*     R26: OH + CO = H + CO2 */
    rf[26] = exp(29.38143762162222 - alogt * .664 - ti * 166.9672903229254);
    eqk = eg[0] * eg[9] / eg[4] / eg[8];
    rb[26] = rf[26] / max(eqk,1e-200);
    /*     R27: HO2 + CO = OH + CO2 */
    rf[27] = exp(alogt * 2.18 + 11.96400108433045 - ti * 9027.70701746016);
    eqk = eg[4] * eg[9] / eg[6] / eg[8];
    rb[27] = rf[27] / max(eqk,1e-200);
    /*     R28: HCO = H + CO */
    rf[28] = exp(alogt * .66 + 27.06890219777501 - ti * 7482.831847805605);
    eqk = eg[0] * eg[8] / eg[11] * pfac1;
    rb[28] = rf[28] / max(eqk,1e-200);
    /*     R29: O2 + HCO = HO2 + CO */
    rf[29] = exp(29.65653431558283 - ti * 206.3188337323671);
    eqk = eg[6] * eg[8] / eg[3] / eg[11];
    rb[29] = rf[29] / max(eqk,1e-200);
    /*     R30: H + HCO = H2 + CO */
    rf[30] = 7.34e13;
    eqk = eg[1] * eg[8] / eg[0] / eg[11];
    rb[30] = rf[30] / max(eqk,1e-200);
    /*     R31: O + HCO = OH + CO */
    rf[31] = 3.02e13;
    eqk = eg[4] * eg[8] / eg[2] / eg[11];
    rb[31] = rf[31] / max(eqk,1e-200);
    /*     R32: O + HCO = H + CO2 */
    rf[32] = 3e13;
    eqk = eg[0] * eg[9] / eg[2] / eg[11];
    rb[32] = rf[32] / max(eqk,1e-200);
    /*     R33: OH + HCO = H2O + CO */
    rf[33] = 1.02e14;
    eqk = eg[5] * eg[8] / eg[4] / eg[11];
    rb[33] = rf[33] / max(eqk,1e-200);
    /*     R34: HO2 + HCO = H + OH + CO2 */
    rf[34] = 3e13;
    rb[34] = 0.;
    /*     R35: 2HCO = H2 + 2CO */
    rf[35] = 3e12;
    rb[35] = 0.;
    /*     R36: HCO + CH3 = CO + CH4 */
    rf[36] = 2.65e13;
    eqk = eg[8] * eg[15] / eg[11] / eg[16];
    rb[36] = rf[36] / max(eqk,1e-200);
    /*     R37: O2 + CH2O = HO2 + HCO */
    rf[37] = exp(36.62692987719255 - ti * 26881.83438532451);
    eqk = eg[6] * eg[11] / eg[3] / eg[10];
    rb[37] = rf[37] / max(eqk,1e-200);
    /*     R38: 2HCO = CO + CH2O */
    rf[38] = 1.8e13;
    eqk = eg[8] * eg[10] / eg[11] / eg[11];
    rb[38] = rf[38] / max(eqk,1e-200);
    /*     R39: H + HCO = CH2O */
    rf[39] = exp(alogt * .48 + 27.7171988121696 + ti * 130.8363335863791);
    eqk = eg[10] / eg[0] / eg[11] / pfac1;
    rb[39] = rf[39] / max(eqk,1e-200);
    /*     R40: H2 + CO = CH2O */
    rf[40] = exp(alogt * 1.5 + 17.57671067365784 - ti * 40056.04674413761);
    eqk = eg[10] / eg[1] / eg[8] / pfac1;
    rb[40] = rf[40] / max(eqk,1e-200);
    /*     R41: OH + CH2O = H2O + HCO */
    rf[41] = exp(alogt * 1.63 + 18.17478020551554 + ti * 530.8935843601153);
    eqk = eg[5] * eg[11] / eg[4] / eg[10];
    rb[41] = rf[41] / max(eqk,1e-200);
    /*     R42: H + CH2O = H2 + HCO */
    rf[42] = exp(alogt * 1.9 + 17.8655548612898 - ti * 1378.81366933338);
    eqk = eg[1] * eg[11] / eg[0] / eg[10];
    rb[42] = rf[42] / max(eqk,1e-200);
    /*     R43: O + CH2O = OH + HCO */
    rf[43] = exp(alogt * 1.15 + 22.55744602205842 - ti * 1137.269668866219);
    eqk = eg[4] * eg[11] / eg[2] / eg[10];
    rb[43] = rf[43] / max(eqk,1e-200);
    /*     R44: CH2O + CH3 = HCO + CH4 */
    rf[44] = exp(alogt * 3.36 + 3.6454498961866 - ti * 2169.870270863334);
    eqk = eg[11] * eg[15] / eg[10] / eg[16];
    rb[44] = rf[44] / max(eqk,1e-200);
    /*     R45: HO2 + CH2O = H2O2 + HCO */
    rf[45] = exp(alogt * 2.7 + 9.84161214881804 - ti * 5797.056011211875);
    eqk = eg[7] * eg[11] / eg[6] / eg[10];
    rb[45] = rf[45] / max(eqk,1e-200);
    /*     R46: CH3O = H + CH2O */
    rf[46] = exp(31.85052882110465 - ti * 13169.1801921367);
    eqk = eg[0] * eg[10] / eg[14] * pfac1;
    rb[46] = rf[46] / max(eqk,1e-200);
    /*     R47: O2 + CH3O = HO2 + CH2O */
    rf[47] = exp(alogt * 9.5 - 42.27206804249851 + ti * 2768.194888687198);
    eqk = eg[6] * eg[10] / eg[3] / eg[14];
    rb[47] = rf[47] / max(eqk,1e-200);
    /*     R48: CH2O + CH3O = HCO + CH3OH */
    rf[48] = exp(27.21853139288342 - ti * 1154.379035565976);
    eqk = eg[11] * eg[12] / eg[10] / eg[14];
    rb[48] = rf[48] / max(eqk,1e-200);
    /*     R49: CH3OH + CH3 = CH3O + CH4 */
    rf[49] = exp(alogt * 3.1 + 2.667228206581955 - ti * 3489.807590082843);
    eqk = eg[14] * eg[15] / eg[12] / eg[16];
    rb[49] = rf[49] / max(eqk,1e-200);
    /*     R50: CH3O + CH3 = CH2O + CH4 */
    rf[50] = 1.2e13;
    eqk = eg[10] * eg[15] / eg[14] / eg[16];
    rb[50] = rf[50] / max(eqk,1e-200);
    /*     R51: H + CH3O = H2 + CH2O */
    rf[51] = 2e13;
    eqk = eg[1] * eg[10] / eg[0] / eg[14];
    rb[51] = rf[51] / max(eqk,1e-200);
    /*     R52: HO2 + CH3O = H2O2 + CH2O */
    rf[52] = 3.01e11;
    eqk = eg[7] * eg[10] / eg[6] / eg[14];
    rb[52] = rf[52] / max(eqk,1e-200);
    /*     R53: H + CH2O = CH2OH */
    rf[53] = exp(alogt * .454 + 27.01483497650473 - ti * 1811.580003503711);
    eqk = eg[13] / eg[0] / eg[10] / pfac1;
    rb[53] = rf[53] / max(eqk,1e-200);
    /*     R54: O2 + CH2OH = HO2 + CH2O */
    rf[54] = exp(34.95088604573752 - alogt * 1.);
    eqk = eg[6] * eg[10] / eg[3] / eg[13];
    rb[54] = rf[54] / max(eqk,1e-200);
    /*     R55: O2 + CH2OH = HO2 + CH2O */
    rf[55] = exp(33.1158180494192 - ti * 2524.638021549477);
    eqk = eg[6] * eg[10] / eg[3] / eg[13];
    rb[55] = rf[55] / max(eqk,1e-200);
    /*     R56: H + CH2OH = H2 + CH2O */
    rf[56] = 6e12;
    eqk = eg[1] * eg[10] / eg[0] / eg[13];
    rb[56] = rf[56] / max(eqk,1e-200);
    /*     R57: HO2 + CH2OH = H2O2 + CH2O */
    rf[57] = 1.2e13;
    eqk = eg[7] * eg[10] / eg[6] / eg[13];
    rb[57] = rf[57] / max(eqk,1e-200);
    /*     R58: HCO + CH2OH = 2CH2O */
    rf[58] = 1.8e14;
    eqk = eg[10] * eg[10] / eg[11] / eg[13];
    rb[58] = rf[58] / max(eqk,1e-200);
    /*     R59: CH2OH + CH3O = CH2O + CH3OH */
    rf[59] = 2.4e13;
    eqk = eg[10] * eg[12] / eg[13] / eg[14];
    rb[59] = rf[59] / max(eqk,1e-200);
    /*     R60: HCO + CH3OH = CH2O + CH2OH */
    rf[60] = exp(alogt * 2.9 + 9.172638504792172 - ti * 6597.170512759347);
    eqk = eg[10] * eg[13] / eg[11] / eg[12];
    rb[60] = rf[60] / max(eqk,1e-200);
    /*     R61: OH + CH2OH = H2O + CH2O */
    rf[61] = 2.4e13;
    eqk = eg[5] * eg[10] / eg[4] / eg[13];
    rb[61] = rf[61] / max(eqk,1e-200);
    /*     R62: O + CH2OH = OH + CH2O */
    rf[62] = 4.2e13;
    eqk = eg[4] * eg[10] / eg[2] / eg[13];
    rb[62] = rf[62] / max(eqk,1e-200);
    /*     R63: 2CH2OH = CH2O + CH3OH */
    rf[63] = 3e12;
    eqk = eg[10] * eg[12] / eg[13] / eg[13];
    rb[63] = rf[63] / max(eqk,1e-200);
    /*     R64: CH3OH = OH + CH3 */
    rf[64] = exp(42.18082079778394 - alogt * .615 - ti * 46567.97235339876);
    eqk = eg[4] * eg[16] / eg[12] * pfac1;
    rb[64] = rf[64] / max(eqk,1e-200);
    /*     R65: CH3OH = H2O + CH2(S) */
    rf[65] = exp(42.58468513718147 - alogt * 1.017 - ti * 46151.00702259231);
    eqk = eg[5] * eg[18] / eg[12] * pfac1;
    rb[65] = rf[65] / max(eqk,1e-200);
    /*     R66: CH3OH = H + CH2OH */
    rf[66] = exp(alogt * 5.038 - 4.841398976850956 - ti * 42505.40355220815);
    eqk = eg[0] * eg[13] / eg[12] * pfac1;
    rb[66] = rf[66] / max(eqk,1e-200);
    /*     R67: H + CH3OH = H2 + CH2OH */
    rf[67] = exp(alogt * 2.55 + 12.63460302656933 - ti * 2737.498671961163);
    eqk = eg[1] * eg[13] / eg[0] / eg[12];
    rb[67] = rf[67] / max(eqk,1e-200);
    /*     R68: H + CH3OH = H2 + CH3O */
    rf[68] = exp(alogt * 2.56 + 12.20106010370663 - ti * 5183.131676691172);
    eqk = eg[1] * eg[14] / eg[0] / eg[12];
    rb[68] = rf[68] / max(eqk,1e-200);
    /*     R69: O + CH3OH = OH + CH2OH */
    rf[69] = exp(alogt * 2.5 + 12.86876061860541 - ti * 1549.907336330953);
    eqk = eg[4] * eg[13] / eg[2] / eg[12];
    rb[69] = rf[69] / max(eqk,1e-200);
    /*     R70: OH + CH3OH = H2O + CH2OH */
    rf[70] = exp(alogt * 2.65 + 10.33526996896167 + ti * 405.9448857851232);
    eqk = eg[5] * eg[13] / eg[4] / eg[12];
    rb[70] = rf[70] / max(eqk,1e-200);
    /*     R71: OH + CH3OH = H2O + CH3O */
    rf[71] = exp(alogt * 3.03 + 5.010635294096256 + ti * 383.9543174092587);
    eqk = eg[5] * eg[14] / eg[4] / eg[12];
    rb[71] = rf[71] / max(eqk,1e-200);
    /*     R72: O2 + CH3OH = HO2 + CH2OH */
    rf[72] = exp(30.65144600207291 - ti * 22594.42837703239);
    eqk = eg[6] * eg[13] / eg[3] / eg[12];
    rb[72] = rf[72] / max(eqk,1e-200);
    /*     R73: HO2 + CH3OH = H2O2 + CH2OH */
    rf[73] = exp(alogt * 2.55 + 9.287301413112312 - ti * 5298.871510248354);
    eqk = eg[7] * eg[13] / eg[6] / eg[12];
    rb[73] = rf[73] / max(eqk,1e-200);
    /*     R74: CH3OH + CH3 = CH2OH + CH4 */
    rf[74] = exp(alogt * 3.17 + 3.462606009790799 - ti * 3609.069940313504);
    eqk = eg[13] * eg[15] / eg[12] / eg[16];
    rb[74] = rf[74] / max(eqk,1e-200);
    /*     R75: CH3O = CH2OH */
    rf[75] = exp(26.42704831160261 - ti * 2050.104703965033);
    eqk = eg[13] / eg[14];
    rb[75] = rf[75] / max(eqk,1e-200);
    /*     R76: 2CH3O = CH2O + CH3OH */
    rf[76] = 6.03e13;
    eqk = eg[10] * eg[12] / eg[14] / eg[14];
    rb[76] = rf[76] / max(eqk,1e-200);
    /*     R77: H + CH3 = CH4 */
    rf[77] = exp(37.08037838837523 - alogt * .63 - ti * 192.7319837060892);
    eqk = eg[15] / eg[0] / eg[16] / pfac1;
    rb[77] = rf[77] / max(eqk,1e-200);
    /*     R78: H + CH4 = H2 + CH3 */
    rf[78] = exp(alogt * 2.5 + 13.32775020712928 - ti * 4824.33819266391);
    eqk = eg[1] * eg[16] / eg[0] / eg[15];
    rb[78] = rf[78] / max(eqk,1e-200);
    /*     R79: OH + CH4 = H2O + CH3 */
    rf[79] = exp(alogt * 2.6 + 10.97335737233858 - ti * 1102.044502131424);
    eqk = eg[5] * eg[16] / eg[4] / eg[15];
    rb[79] = rf[79] / max(eqk,1e-200);
    /*     R80: O + CH4 = OH + CH3 */
    rf[80] = exp(alogt * 1.5 + 20.74306846424259 - ti * 4327.663341703309);
    eqk = eg[4] * eg[16] / eg[2] / eg[15];
    rb[80] = rf[80] / max(eqk,1e-200);
    /*     R81: HO2 + CH4 = H2O2 + CH3 */
    rf[81] = exp(alogt * 3.74 + 2.830267833826459 - ti * 10572.58218711471);
    eqk = eg[7] * eg[16] / eg[6] / eg[15];
    rb[81] = rf[81] / max(eqk,1e-200);
    /*     R82: CH4 + CH2 = 2CH3 */
    rf[82] = exp(alogt * 2. + 14.71567190790855 - ti * 4161.601841382136);
    eqk = eg[16] * eg[16] / eg[15] / eg[17];
    rb[82] = rf[82] / max(eqk,1e-200);
    /*     R83: OH + CH3 = H2O + CH2(S) */
    /*     Reaction of PLOG type */
    rf[83] = 1.;
    eqk = eg[5] * eg[18] / eg[4] / eg[16];
    rb[83] = rf[83] / max(eqk,1e-200);
    /*     R84: OH + CH3 = H2 + CH2O */
    /*     Reaction of PLOG type */
    rf[84] = 1.;
    eqk = eg[1] * eg[10] / eg[4] / eg[16];
    rb[84] = rf[84] / max(eqk,1e-200);
    /*     R85: OH + CH3 = H + CH2OH */
    /*     Reaction of PLOG type */
    rf[85] = 1.;
    eqk = eg[0] * eg[13] / eg[4] / eg[16];
    rb[85] = rf[85] / max(eqk,1e-200);
    /*     R86: OH + CH3 = H + CH3O */
    /*     Reaction of PLOG type */
    rf[86] = 1.;
    eqk = eg[0] * eg[14] / eg[4] / eg[16];
    rb[86] = rf[86] / max(eqk,1e-200);
    /*     R87: HO2 + CH3 = OH + CH3O */
    rf[87] = exp(alogt * .269 + 27.63102111592855 + ti * 345.9614590024448);
    eqk = eg[4] * eg[14] / eg[6] / eg[16];
    rb[87] = rf[87] / max(eqk,1e-200);
    /*     R88: HO2 + CH3 = O2 + CH4 */
    rf[88] = exp(alogt * 2.23 + 11.6613454700885 + ti * 1520.720769607837);
    eqk = eg[3] * eg[15] / eg[6] / eg[16];
    rb[88] = rf[88] / max(eqk,1e-200);
    /*     R89: O + CH3 = H + CH2O */
    rf[89] = exp(alogt * .05 + 31.64560070968179 + ti * 68.43746679902908);
    eqk = eg[0] * eg[10] / eg[2] / eg[16];
    rb[89] = rf[89] / max(eqk,1e-200);
    /*     R90: O2 + CH3 = O + CH3O */
    rf[90] = exp(29.65203873747067 - ti * 14251.09602756253);
    eqk = eg[2] * eg[14] / eg[3] / eg[16];
    rb[90] = rf[90] / max(eqk,1e-200);
    /*     R91: O2 + CH3 = OH + CH2O */
    rf[91] = exp(alogt * 3.283 + .9711576333149952 - ti * 4078.571091221549);
    eqk = eg[4] * eg[10] / eg[3] / eg[16];
    rb[91] = rf[91] / max(eqk,1e-200);
    /*     R92: CH2(S) = CH2 */
    rf[92] = exp(30.33907131703076 - ti * 301.9300005839518);
    eqk = eg[17] / eg[18];
    rb[92] = rf[92] / max(eqk,1e-200);
    /*     R93: O + CH2(S) = H2 + CO */
    rf[93] = 1.5e13;
    eqk = eg[1] * eg[8] / eg[2] / eg[18];
    rb[93] = rf[93] / max(eqk,1e-200);
    /*     R94: O + CH2(S) = H + HCO */
    rf[94] = 1.5e13;
    eqk = eg[0] * eg[11] / eg[2] / eg[18];
    rb[94] = rf[94] / max(eqk,1e-200);
    /*     R95: OH + CH2(S) = H + CH2O */
    rf[95] = 3e13;
    eqk = eg[0] * eg[10] / eg[4] / eg[18];
    rb[95] = rf[95] / max(eqk,1e-200);
    /*     R96: H2 + CH2(S) = H + CH3 */
    rf[96] = 7e13;
    eqk = eg[0] * eg[16] / eg[1] / eg[18];
    rb[96] = rf[96] / max(eqk,1e-200);
    /*     R97: O2 + CH2(S) = H + OH + CO */
    rf[97] = 2.8e13;
    rb[97] = 0.;
    /*     R98: O2 + CH2(S) = H2O + CO */
    rf[98] = 1.2e13;
    eqk = eg[5] * eg[8] / eg[3] / eg[18];
    rb[98] = rf[98] / max(eqk,1e-200);
    /*     R99: CH2(S) = CH2 */
    rf[99] = 3e13;
    eqk = eg[17] / eg[18];
    rb[99] = rf[99] / max(eqk,1e-200);
    /*     R100: CH2(S) = CH2 */
    rf[100] = 9e12;
    eqk = eg[17] / eg[18];
    rb[100] = rf[100] / max(eqk,1e-200);
    /*     R101: CH2(S) = CH2 */
    rf[101] = 7e12;
    eqk = eg[17] / eg[18];
    rb[101] = rf[101] / max(eqk,1e-200);
    /*     R102: CO2 + CH2(S) = CO + CH2O */
    rf[102] = 1.4e13;
    eqk = eg[8] * eg[10] / eg[9] / eg[18];
    rb[102] = rf[102] / max(eqk,1e-200);
    /*     R103: H + CH2 = CH3 */
    rf[103] = exp(37.75765221977888 - alogt * .8);
    eqk = eg[16] / eg[0] / eg[17] / pfac1;
    rb[103] = rf[103] / max(eqk,1e-200);
    /*     R104: O2 + CH2 = OH + HCO */
    rf[104] = exp(29.99187511704657 - ti * 754.8250014598796);
    eqk = eg[4] * eg[11] / eg[3] / eg[17];
    rb[104] = rf[104] / max(eqk,1e-200);
    /*     R105: O2 + CH2 = 2H + CO2 */
    rf[105] = exp(28.60180003308677 - ti * 754.8250014598796);
    rb[105] = 0.;
    /*     R106: O + CH2 = 2H + CO */
    rf[106] = 5e13;
    rb[106] = 0.;
    /*     R107: 2CH3 = C2H6 */
    rf[107] = exp(35.36163518199229 - alogt * .69 - ti * 88.01259517022196);
    eqk = eg[19] / eg[16] / eg[16] / pfac1;
    rb[107] = rf[107] / max(eqk,1e-200);
    /*     R108: H + C2H5 = C2H6 */
    rf[108] = exp(40.79452643666405 - alogt * .99 - ti * 795.082334871073);
    eqk = eg[19] / eg[0] / eg[20] / pfac1;
    rb[108] = rf[108] / max(eqk,1e-200);
    /*     R109: H + C2H6 = H2 + C2H5 */
    rf[109] = exp(alogt * 1.9 + 18.56044268632752 - ti * 3789.221507328595);
    eqk = eg[1] * eg[20] / eg[0] / eg[19];
    rb[109] = rf[109] / max(eqk,1e-200);
    /*     R110: O + C2H6 = OH + C2H5 */
    rf[110] = exp(alogt * 2.4 + 15.0824581614516 - ti * 2933.753172340732);
    eqk = eg[4] * eg[20] / eg[2] / eg[19];
    rb[110] = rf[110] / max(eqk,1e-200);
    /*     R111: OH + C2H6 = H2O + C2H5 */
    rf[111] = exp(alogt * 1.9 + 16.51013773873434 - ti * 478.0558342579237);
    eqk = eg[5] * eg[20] / eg[4] / eg[19];
    rb[111] = rf[111] / max(eqk,1e-200);
    /*     R112: O2 + C2H6 = HO2 + C2H5 */
    rf[112] = exp(31.73035321966169 - ti * 26101.84855048264);
    eqk = eg[6] * eg[20] / eg[3] / eg[19];
    rb[112] = rf[112] / max(eqk,1e-200);
    /*     R113: CH3 + C2H6 = CH4 + C2H5 */
    rf[113] = exp(alogt * 4. - .6014799920341214 - ti * 4166.634008058535);
    eqk = eg[15] * eg[20] / eg[16] / eg[19];
    rb[113] = rf[113] / max(eqk,1e-200);
    /*     R114: HO2 + C2H6 = H2O2 + C2H5 */
    rf[114] = exp(alogt * 3.61 + 3.543853682063679 - ti * 8514.426016467442);
    eqk = eg[7] * eg[20] / eg[6] / eg[19];
    rb[114] = rf[114] / max(eqk,1e-200);
    /*     R115: CH3O + C2H6 = CH3OH + C2H5 */
    rf[115] = exp(26.20806277043707 - ti * 3567.80617356703);
    eqk = eg[12] * eg[20] / eg[14] / eg[19];
    rb[115] = rf[115] / max(eqk,1e-200);
    /*     R116: CH2(S) + C2H6 = CH3 + C2H5 */
    rf[116] = 1.2e14;
    eqk = eg[16] * eg[20] / eg[18] / eg[19];
    rb[116] = rf[116] / max(eqk,1e-200);
    /*     R117: H + C2H4 = C2H5 */
    rf[117] = exp(alogt * 1.463 + 20.67920945074949 - ti * 681.8585846520912);
    eqk = eg[20] / eg[0] / eg[21] / pfac1;
    rb[117] = rf[117] / max(eqk,1e-200);
    /*     R118: 2C2H4 = C2H5 + C2H3 */
    rf[118] = exp(33.80896522997915 - ti * 35995.08823628345);
    eqk = eg[20] * eg[22] / eg[21] / eg[21];
    rb[118] = rf[118] / max(eqk,1e-200);
    /*     R119: CH3 + C2H5 = CH4 + C2H4 */
    rf[119] = exp(alogt * 2.45 + 9.375854810453756 + ti * 1469.895886176205);
    eqk = eg[15] * eg[21] / eg[16] / eg[20];
    rb[119] = rf[119] / max(eqk,1e-200);
    /*     R120: 2CH3 = H + C2H5 */
    /*     Reaction of PLOG type */
    rf[120] = 1.;
    eqk = eg[0] * eg[20] / eg[16] / eg[16];
    rb[120] = rf[120] / max(eqk,1e-200);
    /*     R121: H + C2H5 = H2 + C2H4 */
    rf[121] = 2e12;
    eqk = eg[1] * eg[21] / eg[0] / eg[20];
    rb[121] = rf[121] / max(eqk,1e-200);
    /*     R122: O + C2H5 = H + CH3CHO */
    rf[122] = 1.1e14;
    eqk = eg[0] * eg[25] / eg[2] / eg[20];
    rb[122] = rf[122] / max(eqk,1e-200);
    /*     R123: O2 + C2H5 = HO2 + C2H4 */
    /*     Reaction of PLOG type */
    rf[123] = 1.;
    eqk = eg[6] * eg[21] / eg[3] / eg[20];
    rb[123] = rf[123] / max(eqk,1e-200);
    /*     R124: O2 + C2H5 = HO2 + C2H4 */
    rf[124] = exp(alogt * 3.51 + 1.888432356488316 - ti * 7125.548013781263);
    eqk = eg[6] * eg[21] / eg[3] / eg[20];
    rb[124] = rf[124] / max(eqk,1e-200);
    /*     R125: O2 + C2H5 = OH + C2H4O1-2 */
    /*     Reaction of PLOG type */
    rf[125] = 1.;
    eqk = eg[4] * eg[31] / eg[3] / eg[20];
    rb[125] = rf[125] / max(eqk,1e-200);
    /*     R126: O2 + C2H5 = OH + CH3CHO */
    /*     Reaction of PLOG type */
    rf[126] = 1.;
    eqk = eg[4] * eg[25] / eg[3] / eg[20];
    rb[126] = rf[126] / max(eqk,1e-200);
    /*     R127: C2H4O1-2 = HCO + CH3 */
    rf[127] = exp(31.22283885719935 - ti * 28783.9933890034);
    eqk = eg[11] * eg[16] / eg[31] * pfac1;
    rb[127] = rf[127] / max(eqk,1e-200);
    /*     R128: C2H4O1-2 = CH3CHO */
    rf[128] = exp(29.6334466149597 - ti * 27073.05671902768);
    eqk = eg[25] / eg[31];
    rb[128] = rf[128] / max(eqk,1e-200);
    /*     R129: OH + C2H4O1-2 = H2O + C2H3O1-2 */
    rf[129] = exp(30.51021957322659 - ti * 1816.61217018011);
    eqk = eg[5] * eg[32] / eg[4] / eg[31];
    rb[129] = rf[129] / max(eqk,1e-200);
    /*     R130: H + C2H4O1-2 = H2 + C2H3O1-2 */
    rf[130] = exp(32.01304775060243 - ti * 4871.137342754422);
    eqk = eg[1] * eg[32] / eg[0] / eg[31];
    rb[130] = rf[130] / max(eqk,1e-200);
    /*     R131: HO2 + C2H4O1-2 = H2O2 + C2H3O1-2 */
    rf[131] = exp(30.05582384164685 - ti * 15312.88319628276);
    eqk = eg[7] * eg[32] / eg[6] / eg[31];
    rb[131] = rf[131] / max(eqk,1e-200);
    /*     R132: CH3 + C2H4O1-2 = CH4 + C2H3O1-2 */
    rf[132] = exp(27.69867976440236 - ti * 5953.05317818025);
    eqk = eg[15] * eg[32] / eg[16] / eg[31];
    rb[132] = rf[132] / max(eqk,1e-200);
    /*     R133: CH3O + C2H4O1-2 = CH3OH + C2H3O1-2 */
    rf[133] = exp(25.51075757972846 - ti * 3396.712506569458);
    eqk = eg[12] * eg[32] / eg[14] / eg[31];
    rb[133] = rf[133] / max(eqk,1e-200);
    /*     R134: C2H3O1-2 = CH3CO */
    rf[134] = exp(34.37625746541291 - ti * 7045.033346958876);
    eqk = eg[27] / eg[32];
    rb[134] = rf[134] / max(eqk,1e-200);
    /*     R135: C2H3O1-2 = CH2CHO */
    rf[135] = exp(32.23619130191664 - ti * 7045.033346958876);
    eqk = eg[28] / eg[32];
    rb[135] = rf[135] / max(eqk,1e-200);
    /*     R136: CH3CHO = HCO + CH3 */
    rf[136] = exp(51.55296007042564 - alogt * 1.74 - ti * 43455.27533404526);
    eqk = eg[11] * eg[16] / eg[25] * pfac1;
    rb[136] = rf[136] / max(eqk,1e-200);
    /*     R137: CH3CHO = CO + CH4 */
    rf[137] = exp(49.35491883318287 - alogt * 1.74 - ti * 43455.27533404526);
    eqk = eg[8] * eg[15] / eg[25] * pfac1;
    rb[137] = rf[137] / max(eqk,1e-200);
    /*     R138: H + CH3CHO = H2 + CH3CO */
    rf[138] = exp(alogt * 2.58 + 11.78295260218329 - ti * 613.924334520702);
    eqk = eg[1] * eg[27] / eg[0] / eg[25];
    rb[138] = rf[138] / max(eqk,1e-200);
    /*     R139: H + CH3CHO = H2 + CH2CHO */
    rf[139] = exp(alogt * 3.1 + 7.908387159290043 - ti * 2621.758838403981);
    eqk = eg[1] * eg[28] / eg[0] / eg[25];
    rb[139] = rf[139] / max(eqk,1e-200);
    /*     R140: O + CH3CHO = OH + CH3CO */
    rf[140] = exp(29.4127302493031 - ti * 940.0087351513699);
    eqk = eg[4] * eg[27] / eg[2] / eg[25];
    rb[140] = rf[140] / max(eqk,1e-200);
    /*     R141: OH + CH3CHO = H2O + CH3CO */
    rf[141] = exp(ti * 311.4911172691103 + 28.84593386029282);
    eqk = eg[5] * eg[27] / eg[4] / eg[25];
    rb[141] = rf[141] / max(eqk,1e-200);
    /*     R142: O2 + CH3CHO = HO2 + CH3CO */
    rf[142] = exp(31.03554628768338 - ti * 19700.93253810285);
    eqk = eg[6] * eg[27] / eg[3] / eg[25];
    rb[142] = rf[142] / max(eqk,1e-200);
    /*     R143: CH3 + CH3CHO = CH4 + CH3CO */
    rf[143] = exp(alogt * 4.58 - 7.253066464270554 - ti * 989.3239685800821);
    eqk = eg[15] * eg[27] / eg[16] / eg[25];
    rb[143] = rf[143] / max(eqk,1e-200);
    /*     R144: HO2 + CH3CHO = H2O2 + CH3CO */
    rf[144] = exp(28.73296119468933 - ti * 5998.342678267843);
    eqk = eg[7] * eg[27] / eg[6] / eg[25];
    rb[144] = rf[144] / max(eqk,1e-200);
    /*     R145: OH + CH3CHO = H2O + CH2CHO */
    rf[145] = exp(alogt * 2.4 + 12.05524975579559 - ti * 410.1215841265345);
    eqk = eg[5] * eg[28] / eg[4] / eg[25];
    rb[145] = rf[145] / max(eqk,1e-200);
    /*     R146: CH3CO = CO + CH3 */
    rf[146] = exp(alogt * .63 + 27.69867976440236 - ti * 8504.361683114643);
    eqk = eg[8] * eg[16] / eg[27] * pfac1;
    rb[146] = rf[146] / max(eqk,1e-200);
    /*     R147: H + CH3CO = H2 + CH2CO */
    rf[147] = 2e13;
    eqk = eg[1] * eg[29] / eg[0] / eg[27];
    rb[147] = rf[147] / max(eqk,1e-200);
    /*     R148: O + CH3CO = OH + CH2CO */
    rf[148] = 2e13;
    eqk = eg[4] * eg[29] / eg[2] / eg[27];
    rb[148] = rf[148] / max(eqk,1e-200);
    /*     R149: CH3 + CH3CO = CH4 + CH2CO */
    rf[149] = 5e13;
    eqk = eg[15] * eg[29] / eg[16] / eg[27];
    rb[149] = rf[149] / max(eqk,1e-200);
    /*     R150: CH2CHO = H + CH2CO */
    rf[150] = exp(34.8964508391825 - alogt * .15 - ti * 22946.68004438034);
    eqk = eg[0] * eg[29] / eg[28] * pfac1;
    rb[150] = rf[150] / max(eqk,1e-200);
    /*     R151: CH2CHO = CO + CH3 */
    rf[151] = exp(alogt * .29 + 28.70602353895752 - ti * 20279.63170588876);
    eqk = eg[8] * eg[16] / eg[28] * pfac1;
    rb[151] = rf[151] / max(eqk,1e-200);
    /*     R152: O2 + CH2CHO = HO2 + CH2CO */
    /*     Reaction of PLOG type */
    rf[152] = 1.;
    eqk = eg[6] * eg[29] / eg[3] / eg[28];
    rb[152] = rf[152] / max(eqk,1e-200);
    /*     R153: O2 + CH2CHO = OH + CO + CH2O */
    /*     Reaction of PLOG type */
    rf[153] = 1.;
    rb[153] = 0.;
    /*     R154: CO + CH2 = CH2CO */
    rf[154] = 8.1e11;
    eqk = eg[29] / eg[8] / eg[17] / pfac1;
    rb[154] = rf[154] / max(eqk,1e-200);
    /*     R155: CH3CO = H + CH2CO */
    rf[155] = exp(alogt * 1.917 + 18.3601873635234 - ti * 22638.30887045059);
    eqk = eg[0] * eg[29] / eg[27] * pfac1;
    rb[155] = rf[155] / max(eqk,1e-200);
    /*     R156: H + CH2CO = H2 + HCCO */
    rf[156] = exp(34.87596266226556 - alogt * .171 - ti * 4419.852635214943);
    eqk = eg[1] * eg[30] / eg[0] / eg[29];
    rb[156] = rf[156] / max(eqk,1e-200);
    /*     R157: H + CH2CO = CO + CH3 */
    rf[157] = exp(31.97534588341842 - alogt * .171 - ti * 2105.055964071312);
    eqk = eg[8] * eg[16] / eg[0] / eg[29];
    rb[157] = rf[157] / max(eqk,1e-200);
    /*     R158: O + CH2CO = CO2 + CH2 */
    rf[158] = exp(28.19063690386397 - ti * 679.3425013138916);
    eqk = eg[9] * eg[17] / eg[2] / eg[29];
    rb[158] = rf[158] / max(eqk,1e-200);
    /*     R159: O + CH2CO = OH + HCCO */
    rf[159] = exp(29.93360620892259 - ti * 4025.733341119358);
    eqk = eg[4] * eg[30] / eg[2] / eg[29];
    rb[159] = rf[159] / max(eqk,1e-200);
    /*     R160: OH + CH2CO = H2O + HCCO */
    rf[160] = exp(29.93360620892259 - ti * 1006.433335279839);
    eqk = eg[5] * eg[30] / eg[4] / eg[29];
    rb[160] = rf[160] / max(eqk,1e-200);
    /*     R161: OH + CH2CO = CO + CH2OH */
    rf[161] = exp(ti * 508.2488343163189 + 28.32416829648849);
    eqk = eg[8] * eg[13] / eg[4] / eg[29];
    rb[161] = rf[161] / max(eqk,1e-200);
    /*     R162: CH3 + CH2CO = CO + C2H5 */
    rf[162] = exp(alogt * 2.312 + 10.77247701129227 - ti * 4764.455409214759);
    eqk = eg[8] * eg[20] / eg[16] / eg[29];
    rb[162] = rf[162] / max(eqk,1e-200);
    /*     R163: CH2(S) + CH2CO = CO + C2H4 */
    rf[163] = 1.6e14;
    eqk = eg[8] * eg[21] / eg[18] / eg[29];
    rb[163] = rf[163] / max(eqk,1e-200);
    /*     R164: OH + HCCO = H2 + 2CO */
    rf[164] = 1e14;
    rb[164] = 0.;
    /*     R165: O + HCCO = H + 2CO */
    rf[165] = 8e13;
    rb[165] = 0.;
    /*     R166: H + HCCO = CO + CH2(S) */
    rf[166] = 1e14;
    eqk = eg[8] * eg[18] / eg[0] / eg[30];
    rb[166] = rf[166] / max(eqk,1e-200);
    /*     R167: O2 + HCCO = OH + 2CO */
    rf[167] = exp(25.97553926499304 - alogt * .02 - ti * 513.2810009927181);
    rb[167] = 0.;
    /*     R168: O2 + HCCO = H + CO + CO2 */
    rf[168] = exp(29.19546166243191 - alogt * .142 - ti * 578.6991677859077);
    rb[168] = 0.;
    /*     R169: H + C2H3 = C2H4 */
    rf[169] = exp(alogt * .27 + 29.43602581190662 - ti * 140.9006669391775);
    eqk = eg[21] / eg[0] / eg[22] / pfac1;
    rb[169] = rf[169] / max(eqk,1e-200);
    /*     R170: C2H4 = H2 + H2CC */
    rf[170] = exp(alogt * .44 + 29.71046265760839 - ti * 44670.54358639567);
    eqk = eg[1] * eg[42] / eg[21] * pfac1;
    rb[170] = rf[170] / max(eqk,1e-200);
    /*     R171: H + C2H4 = H2 + C2H3 */
    rf[171] = exp(alogt * 1.93 + 17.74143646856141 - ti * 6516.65584593696);
    eqk = eg[1] * eg[22] / eg[0] / eg[21];
    rb[171] = rf[171] / max(eqk,1e-200);
    /*     R172: O + C2H4 = HCO + CH3 */
    rf[172] = exp(alogt * 1.88 + 15.82412719386383 - ti * 92.0886501781053);
    eqk = eg[11] * eg[16] / eg[2] / eg[21];
    rb[172] = rf[172] / max(eqk,1e-200);
    /*     R173: O + C2H4 = H + CH2CHO */
    rf[173] = exp(alogt * 1.88 + 15.62347140653034 - ti * 92.0886501781053);
    eqk = eg[0] * eg[28] / eg[2] / eg[21];
    rb[173] = rf[173] / max(eqk,1e-200);
    /*     R174: OH + C2H4 = H2O + C2H3 */
    rf[174] = exp(alogt * 2.745 + 10.01234195744821 - ti * 1114.876527156242);
    eqk = eg[5] * eg[22] / eg[4] / eg[21];
    rb[174] = rf[174] / max(eqk,1e-200);
    /*     R175: OH + C2H4 = CH2O + CH3 */
    /*     Reaction of PLOG type */
    rf[175] = 1.;
    eqk = eg[10] * eg[16] / eg[4] / eg[21];
    rb[175] = rf[175] / max(eqk,1e-200);
    /*     R176: OH + C2H4 = H + CH3CHO */
    /*     Reaction of PLOG type */
    rf[176] = 1.;
    eqk = eg[0] * eg[25] / eg[4] / eg[21];
    rb[176] = rf[176] / max(eqk,1e-200);
    /*     R177: OH + C2H4 = H + C2H3OH */
    /*     Reaction of PLOG type */
    rf[177] = 1.;
    eqk = eg[0] * eg[26] / eg[4] / eg[21];
    rb[177] = rf[177] / max(eqk,1e-200);
    /*     R178: O2 + C2H3OH = HO2 + CH2CHO */
    rf[178] = exp(alogt * .21 + 26.99802785818835 - ti * 20043.119872098);
    eqk = eg[6] * eg[28] / eg[3] / eg[26];
    rb[178] = rf[178] / max(eqk,1e-200);
    /*     R179: O + C2H3OH = OH + CH2CHO */
    rf[179] = exp(alogt * 1.9 + 14.44411921738665 + ti * 432.7663341703309);
    eqk = eg[4] * eg[28] / eg[2] / eg[26];
    rb[179] = rf[179] / max(eqk,1e-200);
    /*     R180: OH + C2H3OH = H2O + CH2CHO */
    rf[180] = exp(alogt * 1.1 + 21.92623814093876 - ti * 271.9886088593766);
    eqk = eg[5] * eg[28] / eg[4] / eg[26];
    rb[180] = rf[180] / max(eqk,1e-200);
    /*     R181: CH3 + C2H3OH = CH4 + CH2CHO */
    rf[181] = exp(alogt * 5.9 - 17.71264495089867 - ti * 529.3839343571955);
    eqk = eg[15] * eg[28] / eg[16] / eg[26];
    rb[181] = rf[181] / max(eqk,1e-200);
    /*     R182: H + C2H3OH = H2 + CH2CHO */
    rf[182] = exp(alogt * 3.077 + 7.299797366758161 - ti * 3633.22434036022);
    eqk = eg[1] * eg[28] / eg[0] / eg[26];
    rb[182] = rf[182] / max(eqk,1e-200);
    /*     R183: C2H3OH = CH3CHO */
    rf[183] = exp(alogt * 1.67 + 11.9117015849276 - ti * 3426.905506627853);
    eqk = eg[25] / eg[26];
    rb[183] = rf[183] / max(eqk,1e-200);
    /*     R184: C2H3OH = CH3CHO */
    /*     Reaction of PLOG type */
    rf[184] = 1.;
    eqk = eg[25] / eg[26];
    rb[184] = rf[184] / max(eqk,1e-200);
    /*     R185: CH3 + C2H4 = CH4 + C2H3 */
    rf[185] = exp(alogt * 3.7 + 1.890095369948917 - ti * 4780.558342579237);
    eqk = eg[15] * eg[22] / eg[16] / eg[21];
    rb[185] = rf[185] / max(eqk,1e-200);
    /*     R186: O2 + C2H4 = HO2 + C2H3 */
    rf[186] = exp(31.37344133697052 - ti * 28996.90436108186);
    eqk = eg[6] * eg[22] / eg[3] / eg[21];
    rb[186] = rf[186] / max(eqk,1e-200);
    /*     R187: CH3O + C2H4 = CH3OH + C2H3 */
    rf[187] = exp(25.51075757972846 - ti * 3396.712506569458);
    eqk = eg[12] * eg[22] / eg[14] / eg[21];
    rb[187] = rf[187] / max(eqk,1e-200);
    /*     R188: HO2 + C2H4 = OH + C2H4O1-2 */
    rf[188] = exp(27.04672834028068 - ti * 8650.294516730219);
    eqk = eg[4] * eg[31] / eg[6] / eg[21];
    rb[188] = rf[188] / max(eqk,1e-200);
    /*     R189: CH3 + CH2(S) = H + C2H4 */
    rf[189] = 2e13;
    eqk = eg[0] * eg[21] / eg[16] / eg[18];
    rb[189] = rf[189] / max(eqk,1e-200);
    /*     R190: H + C2H2 = C2H3 */
    rf[190] = exp(alogt * 1.266 + 23.56234430045502 - ti * 1363.213952636542);
    eqk = eg[22] / eg[0] / eg[23] / pfac1;
    rb[190] = rf[190] / max(eqk,1e-200);
    /*     R191: O2 + C2H3 = CH2O + HCO */
    rf[191] = exp(67.3055959478895 - alogt * 5.312 - ti * 3272.468311329162);
    eqk = eg[10] * eg[11] / eg[3] / eg[22];
    rb[191] = rf[191] / max(eqk,1e-200);
    /*     R192: O2 + C2H3 = O + CH2CHO */
    rf[192] = exp(34.18210145097196 - alogt * .611 - ti * 2648.127391788313);
    eqk = eg[2] * eg[28] / eg[3] / eg[22];
    rb[192] = rf[192] / max(eqk,1e-200);
    /*     R193: O2 + C2H3 = H + CO + CH2O */
    rf[193] = exp(36.18551009208849 - alogt * 1.26 - ti * 1666.955533223998);
    rb[193] = 0.;
    /*     R194: CH3 + C2H3 = CH4 + C2H2 */
    rf[194] = 3.92e11;
    eqk = eg[15] * eg[23] / eg[16] / eg[22];
    rb[194] = rf[194] / max(eqk,1e-200);
    /*     R195: H + C2H3 = H2 + C2H2 */
    rf[195] = 9e13;
    eqk = eg[1] * eg[23] / eg[0] / eg[22];
    rb[195] = rf[195] / max(eqk,1e-200);
    /*     R196: H + C2H3 = H2 + H2CC */
    rf[196] = 6e13;
    eqk = eg[1] * eg[42] / eg[0] / eg[22];
    rb[196] = rf[196] / max(eqk,1e-200);
    /*     R197: OH + C2H3 = H2O + C2H2 */
    rf[197] = 3.011e13;
    eqk = eg[5] * eg[23] / eg[4] / eg[22];
    rb[197] = rf[197] / max(eqk,1e-200);
    /*     R198: 2C2H3 = C2H4 + C2H2 */
    rf[198] = 9.6e11;
    eqk = eg[21] * eg[23] / eg[22] / eg[22];
    rb[198] = rf[198] / max(eqk,1e-200);
    /*     R199: H + C2H = C2H2 */
    rf[199] = 1e17;
    eqk = eg[23] / eg[0] / eg[24] / pfac1;
    rb[199] = rf[199] / max(eqk,1e-200);
    /*     R200: OH + C2H = H + HCCO */
    rf[200] = 2e13;
    eqk = eg[0] * eg[30] / eg[4] / eg[24];
    rb[200] = rf[200] / max(eqk,1e-200);
    /*     R201: O2 + C2H = CO + HCO */
    rf[201] = exp(31.54304412135669 - ti * 754.8250014598796);
    eqk = eg[8] * eg[11] / eg[3] / eg[24];
    rb[201] = rf[201] / max(eqk,1e-200);
    /*     R202: H2 + C2H = H + C2H2 */
    rf[202] = exp(alogt * 2.5 + 13.10216067008681 - ti * 281.801333878355);
    eqk = eg[0] * eg[23] / eg[1] / eg[24];
    rb[202] = rf[202] / max(eqk,1e-200);
    /*     R203: C2H2 = H2CC */
    rf[203] = exp(34.31563284359648 - alogt * .52 - ti * 25538.24588272593);
    eqk = eg[42] / eg[23];
    rb[203] = rf[203] / max(eqk,1e-200);
    /*     R204: O + C2H2 = CO + CH2 */
    rf[204] = exp(alogt * 1.28 + 20.42148484011513 - ti * 1243.951602405881);
    eqk = eg[8] * eg[17] / eg[2] / eg[23];
    rb[204] = rf[204] / max(eqk,1e-200);
    /*     R205: O + C2H2 = H + HCCO */
    rf[205] = exp(alogt * 1.28 + 21.80777920123502 - ti * 1243.951602405881);
    eqk = eg[0] * eg[30] / eg[2] / eg[23];
    rb[205] = rf[205] / max(eqk,1e-200);
    /*     R206: OH + C2H2 = H2O + C2H */
    rf[206] = exp(alogt * 2.14 + 14.78325457142734 - ti * 8584.876349937029);
    eqk = eg[5] * eg[24] / eg[4] / eg[23];
    rb[206] = rf[206] / max(eqk,1e-200);
    /*     R207: OH + C2H2 = H + CH2CO */
    /*     Reaction of PLOG type */
    rf[207] = 1.;
    eqk = eg[0] * eg[29] / eg[4] / eg[23];
    rb[207] = rf[207] / max(eqk,1e-200);
    /*     R208: OH + C2H2 = CO + CH3 */
    /*     Reaction of PLOG type */
    rf[208] = 1.;
    eqk = eg[8] * eg[16] / eg[4] / eg[23];
    rb[208] = rf[208] / max(eqk,1e-200);
    /*     R209: HCO + C2H2 = CO + C2H3 */
    rf[209] = exp(alogt * 2. + 16.11809565095832 - ti * 3019.300005839518);
    eqk = eg[8] * eg[22] / eg[11] / eg[23];
    rb[209] = rf[209] / max(eqk,1e-200);
    /*     R210: CH2 + C2H2 = H + C3H3 */
    rf[210] = exp(30.11592776571655 - ti * 3331.294339776268);
    eqk = eg[0] * eg[41] / eg[17] / eg[23];
    rb[210] = rf[210] / max(eqk,1e-200);
    /*     R211: CH2(S) + C2H2 = H + C3H3 */
    rf[211] = 2e13;
    eqk = eg[0] * eg[41] / eg[18] / eg[23];
    rb[211] = rf[211] / max(eqk,1e-200);
    /*     R212: C2H2 + HCCO = CO + C3H3 */
    rf[212] = exp(25.3284360229345 - ti * 1509.650002919759);
    eqk = eg[8] * eg[41] / eg[23] / eg[30];
    rb[212] = rf[212] / max(eqk,1e-200);
    /*     R213: H2CC = C2H2 */
    rf[213] = 1e14;
    eqk = eg[23] / eg[42];
    rb[213] = rf[213] / max(eqk,1e-200);
    /*     R214: OH + H2CC = H + CH2CO */
    rf[214] = 2e13;
    eqk = eg[0] * eg[29] / eg[4] / eg[42];
    rb[214] = rf[214] / max(eqk,1e-200);
    /*     R215: O2 + H2CC = 2HCO */
    rf[215] = 1e13;
    eqk = eg[11] * eg[11] / eg[3] / eg[42];
    rb[215] = rf[215] / max(eqk,1e-200);
    /*     R216: HCO + C2H3 = C2H3CHO */
    rf[216] = 1.81e13;
    eqk = eg[33] / eg[11] / eg[22] / pfac1;
    rb[216] = rf[216] / max(eqk,1e-200);
    /*     R217: C3H8 = CH3 + C2H5 */
    rf[217] = exp(85.45029065915327 - alogt * 5.84 - ti * 49003.23909477538);
    eqk = eg[16] * eg[20] / eg[34] * pfac1;
    rb[217] = rf[217] / max(eqk,1e-200);
    /*     R218: H + NC3H7 = C3H8 */
    rf[218] = 1e14;
    eqk = eg[34] / eg[0] / eg[35] / pfac1;
    rb[218] = rf[218] / max(eqk,1e-200);
    /*     R219: O2 + C3H8 = HO2 + NC3H7 */
    rf[219] = exp(31.72536567815065 - ti * 26313.1995508914);
    eqk = eg[6] * eg[35] / eg[3] / eg[34];
    rb[219] = rf[219] / max(eqk,1e-200);
    /*     R220: H + C3H8 = H2 + NC3H7 */
    rf[220] = exp(alogt * 2.69 + 12.76282720118456 - ti * 3245.747506277482);
    eqk = eg[1] * eg[35] / eg[0] / eg[34];
    rb[220] = rf[220] / max(eqk,1e-200);
    /*     R221: O + C3H8 = OH + NC3H7 */
    rf[221] = exp(alogt * 2.4 + 15.12654243458362 - ti * 2770.207755357758);
    eqk = eg[4] * eg[35] / eg[2] / eg[34];
    rb[221] = rf[221] / max(eqk,1e-200);
    /*     R222: OH + C3H8 = H2O + NC3H7 */
    rf[222] = exp(alogt * .97 + 23.07844338005963 - ti * 798.1016348769126);
    eqk = eg[5] * eg[35] / eg[4] / eg[34];
    rb[222] = rf[222] / max(eqk,1e-200);
    /*     R223: HO2 + C3H8 = H2O2 + NC3H7 */
    rf[223] = exp(alogt * 3.59 + 3.708682081410116 - ti * 8635.198016701022);
    eqk = eg[7] * eg[35] / eg[6] / eg[34];
    rb[223] = rf[223] / max(eqk,1e-200);
    /*     R224: CH3 + C3H8 = CH4 + NC3H7 */
    rf[224] = exp(alogt * 3.65 - .1009259185899605 - ti * 3600.012040295986);
    eqk = eg[15] * eg[35] / eg[16] / eg[34];
    rb[224] = rf[224] / max(eqk,1e-200);
    /*     R225: C2H3 + C3H8 = C2H4 + NC3H7 */
    rf[225] = exp(25.3284360229345 - ti * 5233.453343455165);
    eqk = eg[21] * eg[35] / eg[22] / eg[34];
    rb[225] = rf[225] / max(eqk,1e-200);
    /*     R226: C2H5 + C3H8 = C2H6 + NC3H7 */
    rf[226] = exp(25.3284360229345 - ti * 5233.453343455165);
    eqk = eg[19] * eg[35] / eg[20] / eg[34];
    rb[226] = rf[226] / max(eqk,1e-200);
    /*     R227: C3H8 + C3H5-A = NC3H7 + C3H6 */
    rf[227] = exp(27.40034929819355 - ti * 10315.94168661835);
    eqk = eg[35] * eg[36] / eg[34] / eg[37];
    rb[227] = rf[227] / max(eqk,1e-200);
    /*     R228: CH3O + C3H8 = CH3OH + NC3H7 */
    rf[228] = exp(26.42704831160261 - ti * 3522.516673479438);
    eqk = eg[12] * eg[35] / eg[14] / eg[34];
    rb[228] = rf[228] / max(eqk,1e-200);
    /*     R229: CH3 + C2H4 = NC3H7 */
    rf[229] = exp(alogt * 2.48 + 9.775654181026242 - ti * 3084.718172632708);
    eqk = eg[35] / eg[16] / eg[21] / pfac1;
    rb[229] = rf[229] / max(eqk,1e-200);
    /*     R230: H + C3H6 = NC3H7 */
    rf[230] = exp(alogt * .51 + 26.24472675480866 - ti * 1318.42766921659);
    eqk = eg[35] / eg[0] / eg[36] / pfac1;
    rb[230] = rf[230] / max(eqk,1e-200);
    /*     R231: O2 + NC3H7 = HO2 + C3H6 */
    rf[231] = exp(-42.65050447821876 - ti * 1509.650002919759);
    eqk = eg[6] * eg[36] / eg[3] / eg[35];
    rb[231] = rf[231] / max(eqk,1e-200);
    /*     R232: CH3 + C2H3 = C3H6 */
    rf[232] = 2.5e13;
    eqk = eg[36] / eg[16] / eg[22] / pfac1;
    rb[232] = rf[232] / max(eqk,1e-200);
    /*     R233: H + C3H5-A = C3H6 */
    rf[233] = 2e14;
    eqk = eg[36] / eg[0] / eg[37] / pfac1;
    rb[233] = rf[233] / max(eqk,1e-200);
    /*     R234: C3H6 = H + C3H5-S */
    rf[234] = exp(160.9208896041644 - alogt * 16.09 - ti * 70450.33346958876);
    eqk = eg[0] * eg[38] / eg[36] * pfac1;
    rb[234] = rf[234] / max(eqk,1e-200);
    /*     R235: C3H6 = H + C3H5-T */
    rf[235] = exp(165.2098732664828 - alogt * 16.58 - ti * 70098.08180224081);
    eqk = eg[0] * eg[39] / eg[36] * pfac1;
    rb[235] = rf[235] / max(eqk,1e-200);
    /*     R236: O + C3H6 = HCO + C2H5 */
    rf[236] = exp(alogt * 1.76 + 16.5755204979972 + ti * 611.9114678501423);
    eqk = eg[11] * eg[20] / eg[2] / eg[36];
    rb[236] = rf[236] / max(eqk,1e-200);
    /*     R237: O + C3H6 = H + CH3 + CH2CO */
    rf[237] = exp(alogt * 1.76 + 17.03438638283248 - ti * 38.2444667406339);
    rb[237] = 0.;
    /*     R238: O + C3H6 = OH + C3H5-A */
    rf[238] = exp(alogt * .7 + 26.98475752126745 - ti * 2960.926872393287);
    eqk = eg[4] * eg[37] / eg[2] / eg[36];
    rb[238] = rf[238] / max(eqk,1e-200);
    /*     R239: O + C3H6 = OH + C3H5-S */
    rf[239] = exp(alogt * .7 + 25.51075757972846 - ti * 4508.31812538604);
    eqk = eg[4] * eg[38] / eg[2] / eg[36];
    rb[239] = rf[239] / max(eqk,1e-200);
    /*     R240: O + C3H6 = OH + C3H5-T */
    rf[240] = exp(alogt * .7 + 24.82259794067955 - ti * 3840.549607427867);
    eqk = eg[4] * eg[39] / eg[2] / eg[36];
    rb[240] = rf[240] / max(eqk,1e-200);
    /*     R241: OH + C3H6 = H2O + C3H5-A */
    rf[241] = exp(alogt * 2.2 + 14.49354410071417 - ti * 271.7370005255567);
    eqk = eg[5] * eg[37] / eg[4] / eg[36];
    rb[241] = rf[241] / max(eqk,1e-200);
    /*     R242: OH + C3H6 = H2O + C3H5-S */
    rf[242] = exp(alogt * 2. + 14.56219850545225 - ti * 1397.935902703697);
    eqk = eg[5] * eg[38] / eg[4] / eg[36];
    rb[242] = rf[242] / max(eqk,1e-200);
    /*     R243: OH + C3H6 = H2O + C3H5-T */
    rf[243] = exp(alogt * 2. + 13.91987057328852 - ti * 730.1673847455235);
    eqk = eg[5] * eg[39] / eg[4] / eg[36];
    rb[243] = rf[243] / max(eqk,1e-200);
    /*     R244: HO2 + C3H6 = H2O2 + C3H5-A */
    rf[244] = exp(alogt * 2.5 + 10.20359214498647 - ti * 6209.693678676609);
    eqk = eg[7] * eg[37] / eg[6] / eg[36];
    rb[244] = rf[244] / max(eqk,1e-200);
    /*     R245: HO2 + C3H6 = H2O2 + C3H5-S */
    rf[245] = exp(alogt * 2.5 + 9.798127036878302 - ti * 13898.84436021458);
    eqk = eg[7] * eg[38] / eg[6] / eg[36];
    rb[245] = rf[245] / max(eqk,1e-200);
    /*     R246: HO2 + C3H6 = H2O2 + C3H5-T */
    rf[246] = exp(alogt * 2.5 + 9.104979856318357 - ti * 11870.88118962571);
    eqk = eg[7] * eg[39] / eg[6] / eg[36];
    rb[246] = rf[246] / max(eqk,1e-200);
    /*     R247: H + C3H6 = H2 + C3H5-A */
    rf[247] = exp(alogt * 2.5 + 12.06104687347992 - ti * 1254.01593575868);
    eqk = eg[1] * eg[37] / eg[0] / eg[36];
    rb[247] = rf[247] / max(eqk,1e-200);
    /*     R248: H + C3H6 = H2 + C3H5-S */
    rf[248] = exp(alogt * 2.5 + 13.5973545481611 - ti * 6179.500678618214);
    eqk = eg[1] * eg[38] / eg[0] / eg[36];
    rb[248] = rf[248] / max(eqk,1e-200);
    /*     R249: H + C3H6 = H2 + C3H5-T */
    rf[249] = exp(alogt * 2.5 + 12.91164234608868 - ti * 4928.504042865374);
    eqk = eg[1] * eg[39] / eg[0] / eg[36];
    rb[249] = rf[249] / max(eqk,1e-200);
    /*     R250: H + C3H6 = CH3 + C2H4 */
    /*     Reaction of PLOG type */
    rf[250] = 1.;
    eqk = eg[16] * eg[21] / eg[0] / eg[36];
    rb[250] = rf[250] / max(eqk,1e-200);
    /*     R251: O2 + C3H6 = HO2 + C3H5-A */
    rf[251] = exp(29.01731547704844 - ti * 20078.3450388328);
    eqk = eg[6] * eg[37] / eg[3] / eg[36];
    rb[251] = rf[251] / max(eqk,1e-200);
    /*     R252: O2 + C3H6 = HO2 + C3H5-S */
    rf[252] = exp(28.32416829648849 - ti * 31652.32839455095);
    eqk = eg[6] * eg[38] / eg[3] / eg[36];
    rb[252] = rf[252] / max(eqk,1e-200);
    /*     R253: O2 + C3H6 = HO2 + C3H5-T */
    rf[253] = exp(27.96749335254976 - ti * 30545.25172574313);
    eqk = eg[6] * eg[39] / eg[3] / eg[36];
    rb[253] = rf[253] / max(eqk,1e-200);
    /*     R254: CH3 + C3H6 = CH4 + C3H5-A */
    rf[254] = exp(alogt * 3.5 + .7929925155296614 - ti * 2855.754588856544);
    eqk = eg[15] * eg[37] / eg[16] / eg[36];
    rb[254] = rf[254] / max(eqk,1e-200);
    /*     R255: CH3 + C3H6 = CH4 + C3H5-S */
    rf[255] = exp(alogt * 3.5 + .2986220124901153 - ti * 6466.334179172968);
    eqk = eg[15] * eg[38] / eg[16] / eg[36];
    rb[255] = rf[255] / max(eqk,1e-200);
    /*     R256: CH3 + C3H6 = CH4 + C3H5-T */
    rf[256] = exp(alogt * 3.5 - .1743533871447778 - ti * 5867.506344681464);
    eqk = eg[15] * eg[39] / eg[16] / eg[36];
    rb[256] = rf[256] / max(eqk,1e-200);
    /*     R257: C2H5 + C3H6 = C2H6 + C3H5-A */
    rf[257] = exp(25.3284360229345 - ti * 4931.523342871213);
    eqk = eg[19] * eg[37] / eg[20] / eg[36];
    rb[257] = rf[257] / max(eqk,1e-200);
    /*     R258: CH3 + C2H2 = C3H5-A */
    /*     Reaction of PLOG type */
    rf[258] = 1.;
    eqk = eg[37] / eg[16] / eg[23] / pfac1;
    rb[258] = rf[258] / max(eqk,1e-200);
    /*     R259: O + C3H5-A = H + C2H3CHO */
    rf[259] = 6e13;
    eqk = eg[0] * eg[33] / eg[2] / eg[37];
    rb[259] = rf[259] / max(eqk,1e-200);
    /*     R260: OH + C3H5-A = 2H + C2H3CHO */
    /*     Reaction of PLOG type */
    rf[260] = 1.;
    rb[260] = 0.;
    /*     R261: C2H5 + C3H5-A = C2H4 + C3H6 */
    rf[261] = 4e11;
    eqk = eg[21] * eg[36] / eg[20] / eg[37];
    rb[261] = rf[261] / max(eqk,1e-200);
    /*     R262: O2 + C3H5-A = CH2O + CH3CO */
    /*     Reaction of PLOG type */
    rf[262] = 1.;
    eqk = eg[10] * eg[27] / eg[3] / eg[37];
    rb[262] = rf[262] / max(eqk,1e-200);
    /*     R263: O2 + C3H5-A = OH + C2H3CHO */
    /*     Reaction of PLOG type */
    rf[263] = 1.;
    eqk = eg[4] * eg[33] / eg[3] / eg[37];
    rb[263] = rf[263] / max(eqk,1e-200);
    /*     R264: HCO + C3H5-A = CO + C3H6 */
    rf[264] = 6e13;
    eqk = eg[8] * eg[36] / eg[11] / eg[37];
    rb[264] = rf[264] / max(eqk,1e-200);
    /*     R265: CH3 + C2H3 = H + C3H5-A */
    rf[265] = exp(55.66750733996526 - alogt * 2.83 - ti * 9368.887918120025);
    eqk = eg[0] * eg[37] / eg[16] / eg[22];
    rb[265] = rf[265] / max(eqk,1e-200);
    /*     R266: C3H5-A = C3H5-T */
    /*     Reaction of PLOG type */
    rf[266] = 1.;
    eqk = eg[39] / eg[37];
    rb[266] = rf[266] / max(eqk,1e-200);
    /*     R267: C3H5-A = C3H5-S */
    /*     Reaction of PLOG type */
    rf[267] = 1.;
    eqk = eg[38] / eg[37];
    rb[267] = rf[267] / max(eqk,1e-200);
    /*     R268: CH3 + C2H2 = C3H5-S */
    /*     Reaction of PLOG type */
    rf[268] = 1.;
    eqk = eg[38] / eg[16] / eg[23] / pfac1;
    rb[268] = rf[268] / max(eqk,1e-200);
    /*     R269: O2 + C3H5-S = HCO + CH3CHO */
    rf[269] = 1e11;
    eqk = eg[11] * eg[25] / eg[3] / eg[38];
    rb[269] = rf[269] / max(eqk,1e-200);
    /*     R270: CH3 + C2H2 = C3H5-T */
    /*     Reaction of PLOG type */
    rf[270] = 1.;
    eqk = eg[39] / eg[16] / eg[23] / pfac1;
    rb[270] = rf[270] / max(eqk,1e-200);
    /*     R271: H + C3H4-P = C3H5-T */
    /*     Reaction of PLOG type */
    rf[271] = 1.;
    eqk = eg[39] / eg[0] / eg[40] / pfac1;
    rb[271] = rf[271] / max(eqk,1e-200);
    /*     R272: H + C3H5-S = H2 + C3H4-P */
    rf[272] = 3.34e12;
    eqk = eg[1] * eg[40] / eg[0] / eg[38];
    rb[272] = rf[272] / max(eqk,1e-200);
    /*     R273: O + C3H5-S = HCO + C2H4 */
    rf[273] = 6e13;
    eqk = eg[11] * eg[21] / eg[2] / eg[38];
    rb[273] = rf[273] / max(eqk,1e-200);
    /*     R274: OH + C3H5-S = H + HCO + C2H4 */
    rf[274] = 5e12;
    rb[274] = 0.;
    /*     R275: HO2 + C3H5-S = OH + HCO + C2H4 */
    rf[275] = 2e13;
    rb[275] = 0.;
    /*     R276: HCO + C3H5-S = CO + C3H6 */
    rf[276] = 9e13;
    eqk = eg[8] * eg[36] / eg[11] / eg[38];
    rb[276] = rf[276] / max(eqk,1e-200);
    /*     R277: CH3 + C3H5-S = CH4 + C3H4-P */
    rf[277] = 1e11;
    eqk = eg[15] * eg[40] / eg[16] / eg[38];
    rb[277] = rf[277] / max(eqk,1e-200);
    /*     R278: C3H5-T = C3H5-S */
    /*     Reaction of PLOG type */
    rf[278] = 1.;
    eqk = eg[38] / eg[39];
    rb[278] = rf[278] / max(eqk,1e-200);
    /*     R279: O2 + C3H5-T = CH2O + CH3CO */
    rf[279] = 1e11;
    eqk = eg[10] * eg[27] / eg[3] / eg[39];
    rb[279] = rf[279] / max(eqk,1e-200);
    /*     R280: H + C3H5-T = H2 + C3H4-P */
    rf[280] = 3.34e12;
    eqk = eg[1] * eg[40] / eg[0] / eg[39];
    rb[280] = rf[280] / max(eqk,1e-200);
    /*     R281: CH3 + C3H5-T = CH4 + C3H4-P */
    rf[281] = 1e11;
    eqk = eg[15] * eg[40] / eg[16] / eg[39];
    rb[281] = rf[281] / max(eqk,1e-200);
    /*     R282: O + C3H5-T = CH3 + CH2CO */
    rf[282] = 6e13;
    eqk = eg[16] * eg[29] / eg[2] / eg[39];
    rb[282] = rf[282] / max(eqk,1e-200);
    /*     R283: OH + C3H5-T = H + CH3 + CH2CO */
    rf[283] = 5e12;
    rb[283] = 0.;
    /*     R284: HO2 + C3H5-T = OH + CH3 + CH2CO */
    rf[284] = 2e13;
    rb[284] = 0.;
    /*     R285: HCO + C3H5-T = CO + C3H6 */
    rf[285] = 9e13;
    eqk = eg[8] * eg[36] / eg[11] / eg[39];
    rb[285] = rf[285] / max(eqk,1e-200);
    /*     R286: H + C3H4-P = C3H5-S */
    /*     Reaction of PLOG type */
    rf[286] = 1.;
    eqk = eg[38] / eg[0] / eg[40] / pfac1;
    rb[286] = rf[286] / max(eqk,1e-200);
    /*     R287: H + C3H4-P = C3H5-A */
    /*     Reaction of PLOG type */
    rf[287] = 1.;
    eqk = eg[37] / eg[0] / eg[40] / pfac1;
    rb[287] = rf[287] / max(eqk,1e-200);
    /*     R288: H + C3H4-P = H2 + C3H3 */
    rf[288] = exp(alogt * 2. + 14.07787482243177 - ti * 2767.691672019558);
    eqk = eg[1] * eg[41] / eg[0] / eg[40];
    rb[288] = rf[288] / max(eqk,1e-200);
    /*     R289: O + C3H4-P = CH3 + HCCO */
    rf[289] = exp(29.61889546408289 - ti * 1132.237502189819);
    eqk = eg[16] * eg[30] / eg[2] / eg[40];
    rb[289] = rf[289] / max(eqk,1e-200);
    /*     R290: O + C3H4-P = CO + C2H4 */
    rf[290] = exp(29.93360620892259 - ti * 1132.237502189819);
    eqk = eg[8] * eg[21] / eg[2] / eg[40];
    rb[290] = rf[290] / max(eqk,1e-200);
    /*     R291: OH + C3H4-P = H2O + C3H3 */
    rf[291] = exp(alogt * 2. + 13.81551055796427 - ti * 50.32166676399197);
    eqk = eg[5] * eg[41] / eg[4] / eg[40];
    rb[291] = rf[291] / max(eqk,1e-200);
    /*     R292: C2H + C3H4-P = C2H2 + C3H3 */
    rf[292] = 1e13;
    eqk = eg[23] * eg[41] / eg[24] / eg[40];
    rb[292] = rf[292] / max(eqk,1e-200);
    /*     R293: CH3 + C3H4-P = CH4 + C3H3 */
    rf[293] = exp(28.21880778083067 - ti * 3874.768340827382);
    eqk = eg[15] * eg[41] / eg[16] / eg[40];
    rb[293] = rf[293] / max(eqk,1e-200);
    /*     R294: CH3 + C2H2 = H + C3H4-P */
    /*     Reaction of PLOG type */
    rf[294] = 1.;
    eqk = eg[0] * eg[40] / eg[16] / eg[23];
    rb[294] = rf[294] / max(eqk,1e-200);
    /*     R295: H + C3H3 = C3H4-P */
    rf[295] = 1.5e13;
    eqk = eg[40] / eg[0] / eg[41] / pfac1;
    rb[295] = rf[295] / max(eqk,1e-200);
    /*     R296: CH3 + C2H = C3H4-P */
    rf[296] = 8e13;
    eqk = eg[40] / eg[16] / eg[24] / pfac1;
    rb[296] = rf[296] / max(eqk,1e-200);
    /*     R297: HO2 + C3H3 = O2 + C3H4-P */
    rf[297] = 2.5e12;
    eqk = eg[3] * eg[40] / eg[6] / eg[41];
    rb[297] = rf[297] / max(eqk,1e-200);
    /*     R298: HO2 + C3H4-P = OH + CO + C2H4 */
    rf[298] = exp(28.72963340459666 - ti * 9561.116685158473);
    rb[298] = 0.;
    /*     R299: OH + C3H4-P = CH3 + CH2CO */
    rf[299] = exp(alogt * 4.5 - 7.600902459542082 + ti * 503.2166676399197);
    eqk = eg[16] * eg[29] / eg[4] / eg[40];
    rb[299] = rf[299] / max(eqk,1e-200);
    /*     R300: O + C3H4-P = HCO + C2H3 */
    rf[300] = exp(28.79417192573423 - ti * 1011.465501956239);
    eqk = eg[11] * eg[22] / eg[2] / eg[40];
    rb[300] = rf[300] / max(eqk,1e-200);
    /*     R301: O + C3H4-P = OH + C3H3 */
    rf[301] = exp(alogt * 1.5 + 20.45538639179081 - ti * 4327.663341703309);
    eqk = eg[4] * eg[41] / eg[2] / eg[40];
    rb[301] = rf[301] / max(eqk,1e-200);
    /*     R302: C2H3 + C3H4-P = C2H4 + C3H3 */
    rf[302] = exp(27.63102111592855 - ti * 3874.768340827382);
    eqk = eg[21] * eg[41] / eg[22] / eg[40];
    rb[302] = rf[302] / max(eqk,1e-200);
    /*     R303: C3H5-A + C3H4-P = C3H6 + C3H3 */
    rf[303] = exp(27.63102111592855 - ti * 3874.768340827382);
    eqk = eg[36] * eg[41] / eg[37] / eg[40];
    rb[303] = rf[303] / max(eqk,1e-200);
    /*     R304: O + C3H3 = CH2O + C2H */
    rf[304] = 2e13;
    eqk = eg[10] * eg[24] / eg[2] / eg[41];
    rb[304] = rf[304] / max(eqk,1e-200);
    /*     R305: O2 + C3H3 = HCO + CH2CO */
    rf[305] = exp(24.12446321860857 - ti * 1443.22540279129);
    eqk = eg[11] * eg[29] / eg[3] / eg[41];
    rb[305] = rf[305] / max(eqk,1e-200);
    /*     R306: HO2 + C3H3 = OH + CO + C2H3 */
    rf[306] = 8e11;
    rb[306] = 0.;
    /*     R307: HCO + C3H3 = CO + C3H4-P */
    rf[307] = 2.5e13;
    eqk = eg[8] * eg[40] / eg[11] / eg[41];
    rb[307] = rf[307] / max(eqk,1e-200);
    /*     R308: C2H5 + C2H = CH3 + C3H3 */
    rf[308] = 1.81e13;
    eqk = eg[16] * eg[41] / eg[20] / eg[24];
    rb[308] = rf[308] / max(eqk,1e-200);
    /*     R309: C4H6 = H + C4H5-I */
    rf[309] = exp(84.63352952262615 - alogt * 6.27 - ti * 56537.9022593479);
    eqk = eg[0] * eg[48] / eg[43] * pfac1;
    rb[309] = rf[309] / max(eqk,1e-200);
    /*     R310: C4H6 = H2 + C4H4 */
    rf[310] = exp(35.45506712678484 - ti * 47654.61842550039);
    eqk = eg[1] * eg[44] / eg[43] * pfac1;
    rb[310] = rf[310] / max(eqk,1e-200);
    /*     R311: H + C4H6 = H2 + C4H5-I */
    rf[311] = exp(alogt * 2.53 + 13.40754231963799 - ti * 4649.722008992858);
    eqk = eg[1] * eg[48] / eg[0] / eg[43];
    rb[311] = rf[311] / max(eqk,1e-200);
    /*     R312: H + C4H6 = C2H4 + C2H3 */
    /*     Reaction of PLOG type */
    rf[312] = 1.;
    eqk = eg[21] * eg[22] / eg[0] / eg[43];
    rb[312] = rf[312] / max(eqk,1e-200);
    /*     R313: H + C4H6 = CH3 + C3H4-P */
    rf[313] = exp(28.32416829648849 - ti * 3522.516673479438);
    eqk = eg[16] * eg[40] / eg[0] / eg[43];
    rb[313] = rf[313] / max(eqk,1e-200);
    /*     R314: O + C4H6 = OH + C4H5-I */
    rf[314] = exp(alogt * 1.9 + 15.83041357850654 - ti * 1882.0303369733);
    eqk = eg[4] * eg[48] / eg[2] / eg[43];
    rb[314] = rf[314] / max(eqk,1e-200);
    /*     R315: O + C4H6 = C2H2 + C2H4O1-2 */
    rf[315] = exp(alogt * 1.45 + 18.42068074395237 + ti * 432.7663341703309);
    eqk = eg[23] * eg[31] / eg[2] / eg[43];
    rb[315] = rf[315] / max(eqk,1e-200);
    /*     R316: OH + C4H6 = H2O + C4H5-I */
    rf[316] = exp(alogt * 2. + 14.94691266945537 - ti * 216.3831670851655);
    eqk = eg[5] * eg[48] / eg[4] / eg[43];
    rb[316] = rf[316] / max(eqk,1e-200);
    /*     R317: CH3 + C4H6 = CH4 + C4H5-I */
    rf[317] = exp(32.23619130191664 - ti * 9963.690019270409);
    eqk = eg[15] * eg[48] / eg[16] / eg[43];
    rb[317] = rf[317] / max(eqk,1e-200);
    /*     R318: C2H3 + C4H6 = C2H4 + C4H5-I */
    rf[318] = exp(30.84989694079675 - ti * 9963.690019270409);
    eqk = eg[21] * eg[48] / eg[22] / eg[43];
    rb[318] = rf[318] / max(eqk,1e-200);
    /*     R319: C3H5-A + C4H6 = C3H6 + C4H5-I */
    rf[319] = exp(29.24045902836265 - ti * 9812.725018978434);
    eqk = eg[36] * eg[48] / eg[37] / eg[43];
    rb[319] = rf[319] / max(eqk,1e-200);
    /*     R320: C2H3 + C2H2 = H + C4H4 */
    /*     Reaction of PLOG type */
    rf[320] = 1.;
    eqk = eg[0] * eg[44] / eg[22] / eg[23];
    rb[320] = rf[320] / max(eqk,1e-200);
    /*     R321: C2H3 + C2H2 = C4H5-I */
    /*     Reaction of PLOG type */
    rf[321] = 1.;
    eqk = eg[48] / eg[22] / eg[23] / pfac1;
    rb[321] = rf[321] / max(eqk,1e-200);
    /*     R322: 2C2H3 = C4H6 */
    /*     Reaction of PLOG type */
    rf[322] = 1.;
    eqk = eg[43] / eg[22] / eg[22] / pfac1;
    rb[322] = rf[322] / max(eqk,1e-200);
    /*     R323: 2C2H3 = H + C4H5-I */
    /*     Reaction of PLOG type */
    rf[323] = 1.;
    eqk = eg[0] * eg[48] / eg[22] / eg[22];
    rb[323] = rf[323] / max(eqk,1e-200);
    /*     R324: H + C4H5-I = H2 + C4H4 */
    rf[324] = 3e13;
    eqk = eg[1] * eg[44] / eg[0] / eg[48];
    rb[324] = rf[324] / max(eqk,1e-200);
    /*     R325: H + C4H5-I = CH3 + C3H3 */
    rf[325] = exp(30.62675338948254 - ti * 1006.433335279839);
    eqk = eg[16] * eg[41] / eg[0] / eg[48];
    rb[325] = rf[325] / max(eqk,1e-200);
    /*     R326: OH + C4H5-I = H2O + C4H4 */
    rf[326] = 4e12;
    eqk = eg[5] * eg[44] / eg[4] / eg[48];
    rb[326] = rf[326] / max(eqk,1e-200);
    /*     R327: HCO + C4H5-I = CO + C4H6 */
    rf[327] = 5e12;
    eqk = eg[8] * eg[43] / eg[11] / eg[48];
    rb[327] = rf[327] / max(eqk,1e-200);
    /*     R328: HO2 + C4H5-I = O2 + C4H6 */
    rf[328] = 6e11;
    eqk = eg[3] * eg[43] / eg[6] / eg[48];
    rb[328] = rf[328] / max(eqk,1e-200);
    /*     R329: HO2 + C4H5-I = OH + C2H3 + CH2CO */
    rf[329] = 6.6e12;
    rb[329] = 0.;
    /*     R330: H2O2 + C4H5-I = HO2 + C4H6 */
    rf[330] = exp(ti * 299.9171339133922 + 23.21647128954911);
    eqk = eg[6] * eg[43] / eg[7] / eg[48];
    rb[330] = rf[330] / max(eqk,1e-200);
    /*     R331: O2 + C4H5-I = CH2CHO + CH2CO */
    rf[331] = exp(23.79595915163653 - ti * 1258.041669099799);
    eqk = eg[28] * eg[29] / eg[3] / eg[48];
    rb[331] = rf[331] / max(eqk,1e-200);
    /*     R332: C4H5-2 = C4H5-I */
    rf[332] = exp(154.6786663387092 - alogt * 16.89 - ti * 29740.10505751925);
    eqk = eg[48] / eg[49];
    rb[332] = rf[332] / max(eqk,1e-200);
    /*     R333: C4H5-2 = C4H5-I */
    rf[333] = exp(60.99861452933629 - alogt * 3.35 - ti * 8767.54400029032);
    eqk = eg[48] / eg[49];
    rb[333] = rf[333] / max(eqk,1e-200);
    /*     R334: HO2 + C4H5-2 = OH + C2H2 + CH3CO */
    rf[334] = 8e11;
    rb[334] = 0.;
    /*     R335: O2 + C4H5-2 = CH3CO + CH2CO */
    rf[335] = exp(23.79595915163653 - ti * 1258.041669099799);
    eqk = eg[27] * eg[29] / eg[3] / eg[49];
    rb[335] = rf[335] / max(eqk,1e-200);
    /*     R336: C4H6-2 = C4H6 */
    rf[336] = exp(31.03221849759071 - ti * 32709.08339659478);
    eqk = eg[43] / eg[50];
    rb[336] = rf[336] / max(eqk,1e-200);
    /*     R337: H + C4H6-2 = H2 + C4H5-2 */
    rf[337] = exp(alogt * 2.5 + 12.73670089659234 - ti * 1253.0095024234);
    eqk = eg[1] * eg[49] / eg[0] / eg[50];
    rb[337] = rf[337] / max(eqk,1e-200);
    /*     R338: H + C4H6-2 = CH3 + C3H4-P */
    rf[338] = exp(alogt * 2.5 + 12.46843690999767 - ti * 503.2166676399197);
    eqk = eg[16] * eg[40] / eg[0] / eg[50];
    rb[338] = rf[338] / max(eqk,1e-200);
    /*     R339: C4H6-2 = H + C4H5-2 */
    rf[339] = exp(36.14821430734479 - ti * 43930.81508496499);
    eqk = eg[0] * eg[49] / eg[50] * pfac1;
    rb[339] = rf[339] / max(eqk,1e-200);
    /*     R340: CH3 + C4H6-2 = CH4 + C4H5-2 */
    rf[340] = exp(32.57266353853785 - ti * 9309.508351338514);
    eqk = eg[15] * eg[49] / eg[16] / eg[50];
    rb[340] = rf[340] / max(eqk,1e-200);
    /*     R341: H + C4H4 = C4H5-I */
    /*     Reaction of PLOG type */
    rf[341] = 1.;
    eqk = eg[48] / eg[0] / eg[44] / pfac1;
    rb[341] = rf[341] / max(eqk,1e-200);
    /*     R342: H + C4H4 = H2 + C4H3-N */
    rf[342] = exp(alogt * 2.53 + 13.40754231963799 - ti * 6159.372011912617);
    eqk = eg[1] * eg[46] / eg[0] / eg[44];
    rb[342] = rf[342] / max(eqk,1e-200);
    /*     R343: H + C4H4 = H2 + C4H3-I */
    rf[343] = exp(alogt * 2.53 + 12.71589776896258 - ti * 4649.722008992858);
    eqk = eg[1] * eg[45] / eg[0] / eg[44];
    rb[343] = rf[343] / max(eqk,1e-200);
    /*     R344: OH + C4H4 = H2O + C4H3-N */
    rf[344] = exp(alogt * 2. + 17.24949776244942 - ti * 1726.033170004924);
    eqk = eg[5] * eg[46] / eg[4] / eg[44];
    rb[344] = rf[344] / max(eqk,1e-200);
    /*     R345: OH + C4H4 = H2O + C4H3-I */
    rf[345] = exp(alogt * 2. + 16.55635058188948 - ti * 216.3831670851655);
    eqk = eg[5] * eg[45] / eg[4] / eg[44];
    rb[345] = rf[345] / max(eqk,1e-200);
    /*     R346: O + C4H4 = HCO + C3H3 */
    rf[346] = exp(alogt * 1.45 + 20.21244021318042 + ti * 432.7663341703309);
    eqk = eg[11] * eg[41] / eg[2] / eg[44];
    rb[346] = rf[346] / max(eqk,1e-200);
    /*     R347: HCCO + C3H3 = CO + C4H4 */
    rf[347] = 2.5e13;
    eqk = eg[8] * eg[44] / eg[30] / eg[41];
    rb[347] = rf[347] / max(eqk,1e-200);
    /*     R348: CH2 + C3H3 = H + C4H4 */
    rf[348] = 5e13;
    eqk = eg[0] * eg[44] / eg[17] / eg[41];
    rb[348] = rf[348] / max(eqk,1e-200);
    /*     R349: C2H2 + C2H = C4H3-N */
    rf[349] = exp(alogt * .899 + 25.14210644474301 + ti * 182.6676503532908);
    eqk = eg[46] / eg[23] / eg[24] / pfac1;
    rb[349] = rf[349] / max(eqk,1e-200);
    /*     R350: C4H3-N = C4H3-I */
    rf[350] = exp(100.4221459724542 - alogt * 9.49 - ti * 26670.48338491574);
    eqk = eg[45] / eg[46];
    rb[350] = rf[350] / max(eqk,1e-200);
    /*     R351: C4H3-N = C4H3-I */
    rf[351] = exp(46.96799259175507 - alogt * 1.67 - ti * 5434.740010511133);
    eqk = eg[45] / eg[46];
    rb[351] = rf[351] / max(eqk,1e-200);
    /*     R352: H + C4H3-N = C2H2 + H2CC */
    rf[352] = exp(59.40517695824863 - alogt * 3.34 - ti * 5039.211709746156);
    eqk = eg[23] * eg[42] / eg[0] / eg[46];
    rb[352] = rf[352] / max(eqk,1e-200);
    /*     R353: H + C4H3-N = C4H4 */
    rf[353] = exp(108.9146465512801 - alogt * 10.26 - ti * 6577.04184605375);
    eqk = eg[44] / eg[0] / eg[46] / pfac1;
    rb[353] = rf[353] / max(eqk,1e-200);
    /*     R354: H + C4H3-N = H2 + C4H2 */
    rf[354] = 3e13;
    eqk = eg[1] * eg[47] / eg[0] / eg[46];
    rb[354] = rf[354] / max(eqk,1e-200);
    /*     R355: OH + C4H3-N = H2O + C4H2 */
    rf[355] = 2e12;
    eqk = eg[5] * eg[47] / eg[4] / eg[46];
    rb[355] = rf[355] / max(eqk,1e-200);
    /*     R356: C2H2 + C2H = C4H3-I */
    rf[356] = exp(alogt * .899 + 25.14210644474301 + ti * 182.6676503532908);
    eqk = eg[45] / eg[23] / eg[24] / pfac1;
    rb[356] = rf[356] / max(eqk,1e-200);
    /*     R357: H + C4H3-I = C2H2 + H2CC */
    rf[357] = exp(53.98907655604421 - alogt * 2.55 - ti * 5424.675677158334);
    eqk = eg[23] * eg[42] / eg[0] / eg[45];
    rb[357] = rf[357] / max(eqk,1e-200);
    /*     R358: H + C4H3-I = C4H4 */
    rf[358] = exp(100.2349344303661 - alogt * 9.01 - ti * 6098.986011795827);
    eqk = eg[44] / eg[0] / eg[45] / pfac1;
    rb[358] = rf[358] / max(eqk,1e-200);
    /*     R359: H + C4H3-I = H2 + C4H2 */
    rf[359] = 6e13;
    eqk = eg[1] * eg[47] / eg[0] / eg[45];
    rb[359] = rf[359] / max(eqk,1e-200);
    /*     R360: OH + C4H3-I = H2O + C4H2 */
    rf[360] = 4e12;
    eqk = eg[5] * eg[47] / eg[4] / eg[45];
    rb[360] = rf[360] / max(eqk,1e-200);
    /*     R361: O2 + C4H3-I = CH2CO + HCCO */
    rf[361] = exp(38.90314809434585 - alogt * 1.8);
    eqk = eg[29] * eg[30] / eg[3] / eg[45];
    rb[361] = rf[361] / max(eqk,1e-200);
    /*     R362: C2H2 + C2H = H + C4H2 */
    rf[362] = 9.6e13;
    eqk = eg[0] * eg[47] / eg[23] / eg[24];
    rb[362] = rf[362] / max(eqk,1e-200);
    /*     R363: H + C4H2 = C4H3-N */
    rf[363] = exp(96.80388408555425 - alogt * 8.720000000000001 - ti *
            7699.215014890771);
    eqk = eg[46] / eg[0] / eg[47] / pfac1;
    rb[363] = rf[363] / max(eqk,1e-200);
    /*     R364: H + C4H2 = C4H3-I */
    rf[364] = exp(69.1728629696257 - alogt * 4.92 - ti * 5434.740010511133);
    eqk = eg[45] / eg[0] / eg[47] / pfac1;
    rb[364] = rf[364] / max(eqk,1e-200);
    /*     R365: OH + C4H2 = H + H2C4O */
    rf[365] = exp(ti * 206.3188337323671 + 29.51809076496093);
    eqk = eg[0] * eg[51] / eg[4] / eg[47];
    rb[365] = rf[365] / max(eqk,1e-200);
    /*     R366: H + H2C4O = C2H2 + HCCO */
    rf[366] = exp(31.54304412135669 - ti * 1509.650002919759);
    eqk = eg[23] * eg[30] / eg[0] / eg[51];
    rb[366] = rf[366] / max(eqk,1e-200);
    /*     R367: OH + H2C4O = CH2CO + HCCO */
    rf[367] = exp(alogt * 2. + 16.11809565095832 - ti * 1006.433335279839);
    eqk = eg[29] * eg[30] / eg[4] / eg[51];
    rb[367] = rf[367] / max(eqk,1e-200);
    /*     R368: C2H2 + H2CC = C4H4 */
    rf[368] = exp(alogt * 2.055 + 12.7656884334656 + ti * 1207.720002335807);
    eqk = eg[44] / eg[23] / eg[42] / pfac1;
    rb[368] = rf[368] / max(eqk,1e-200);
    /*     R369: C2H4 + H2CC = C4H6 */
    rf[369] = 1e12;
    eqk = eg[43] / eg[21] / eg[42] / pfac1;
    rb[369] = rf[369] / max(eqk,1e-200);
    /*     R370: C2H3 + C3H5-A = 2H + c-C5H6 */
    rf[370] = exp(81.06048188403733 - alogt * 14. - ti * 30765.50966116912);
    rb[370] = 0.;
    /*     R371: C3H5-A + C3H3 = 2H + A1 */
    rf[371] = exp(47.77446845762201 - alogt * 2.54 - ti * 853.9083633181798);
    rb[371] = 0.;
    /*     R372: 2C3H3 = H + A1- */
    rf[372] = 5e12;
    eqk = eg[0] * eg[54] / eg[41] / eg[41];
    rb[372] = rf[372] / max(eqk,1e-200);
    /*     R373: 2C3H3 = A1 */
    rf[373] = 2e12;
    eqk = eg[52] / eg[41] / eg[41] / pfac1;
    rb[373] = rf[373] / max(eqk,1e-200);
    /*     R374: C2H3 + C4H6 = H + H2 + A1 */
    rf[374] = exp(27.0547676868401 - ti * 1630.42200315334);
    rb[374] = 0.;
    /*     R375: C2H2 + C4H5-2 = H + A1 */
    rf[375] = exp(33.84562921435074 - ti * 12580.41669099799);
    eqk = eg[0] * eg[52] / eg[23] / eg[49];
    rb[375] = rf[375] / max(eqk,1e-200);
    /*     R376: C2H4 + C4H5-2 = CH3 + c-C5H6 */
    rf[376] = exp(33.84562921435074 - ti * 12580.41669099799);
    eqk = eg[16] * eg[53] / eg[21] / eg[49];
    rb[376] = rf[376] / max(eqk,1e-200);
    /*     R377: C2H2 + C4H3-N = A1- */
    rf[377] = exp(163.442719608057 - alogt * 17.77 - ti * 15750.68169712949);
    eqk = eg[54] / eg[23] / eg[46] / pfac1;
    rb[377] = rf[377] / max(eqk,1e-200);
    /*     R378: CH3 + C4H3-I = c-C5H6 */
    rf[378] = 1e12;
    eqk = eg[53] / eg[16] / eg[45] / pfac1;
    rb[378] = rf[378] / max(eqk,1e-200);
    /*     R379: H + A1 = H2 + A1- */
    rf[379] = exp(33.1524820337908 - ti * 8051.466682238715);
    eqk = eg[1] * eg[54] / eg[0] / eg[52];
    rb[379] = rf[379] / max(eqk,1e-200);
    /*     R380: HO2 + A1 = H2O2 + A1- */
    rf[380] = exp(29.33576920816697 - ti * 14542.96169479368);
    eqk = eg[7] * eg[54] / eg[6] / eg[52];
    rb[380] = rf[380] / max(eqk,1e-200);
    /*     R381: OH + A1 = H2O + A1- */
    rf[381] = exp(alogt * 4.1 + .1823215567939546 + ti * 151.4682169596158);
    eqk = eg[5] * eg[54] / eg[4] / eg[52];
    rb[381] = rf[381] / max(eqk,1e-200);
    /*     R382: A1- = H + C2H2 + C4H2 */
    rf[382] = exp(alogt * .62 + 29.08963613862807 - ti * 38895.62910855995);
    rb[382] = 0.;
    /*     R383: CH2O + A1- = HCO + A1 */
    rf[383] = exp(alogt * 2.19 + 11.35627165492485 - ti * 19.12223337031695);
    eqk = eg[11] * eg[52] / eg[10] / eg[54];
    rb[383] = rf[383] / max(eqk,1e-200);
    /*     R384: HCO + A1- = CO + A1 */
    rf[384] = exp(alogt * 2.19 + 11.35627165492485 - ti * 19.12223337031695);
    eqk = eg[8] * eg[52] / eg[11] / eg[54];
    rb[384] = rf[384] / max(eqk,1e-200);
    /*     R385: HO2 + A1- = OH + C6H5O */
    rf[385] = 5e12;
    eqk = eg[4] * eg[55] / eg[6] / eg[54];
    rb[385] = rf[385] / max(eqk,1e-200);
    /*     R386: O2 + A1- = O + C6H5O */
    rf[386] = exp(30.88911765395003 - ti * 3079.686005956308);
    eqk = eg[2] * eg[55] / eg[3] / eg[54];
    rb[386] = rf[386] / max(eqk,1e-200);
    /*     R387: O + C6H5O = CO2 + c-C5H5 */
    rf[387] = 1e13;
    eqk = eg[9] * eg[56] / eg[2] / eg[55];
    rb[387] = rf[387] / max(eqk,1e-200);
    /*     R388: O2 + c-C5H6 = HO2 + c-C5H5 */
    rf[388] = exp(31.31990057004249 - ti * 18694.49920282302);
    eqk = eg[6] * eg[56] / eg[3] / eg[53];
    rb[388] = rf[388] / max(eqk,1e-200);
    /*     R389: HO2 + c-C5H6 = H2O2 + c-C5H5 */
    rf[389] = exp(alogt * 2.6 + 9.305650551780508 - ti * 6491.495012554964);
    eqk = eg[7] * eg[56] / eg[6] / eg[53];
    rb[389] = rf[389] / max(eqk,1e-200);
    /*     R390: H + c-C5H6 = C2H2 + C3H5-A */
    rf[390] = exp(84.93946503538729 - alogt * 6.18 - ti * 16550.79619867696);
    eqk = eg[23] * eg[37] / eg[0] / eg[53];
    rb[390] = rf[390] / max(eqk,1e-200);
    /*     R391: CH3 + c-C5H6 = CH4 + c-C5H5 */
    rf[391] = exp(alogt * 4. - 1.714798428091927);
    eqk = eg[15] * eg[56] / eg[16] / eg[53];
    rb[391] = rf[391] / max(eqk,1e-200);
    /*     R392: O + c-C5H5 = H + C5H4O */
    rf[392] = exp(31.69146412647497 - alogt * .02 - ti * 10.06433335279839);
    eqk = eg[0] * eg[57] / eg[2] / eg[56];
    rb[392] = rf[392] / max(eqk,1e-200);
    /*     R393: OH + c-C5H5 = CO + C4H6 */
    rf[393] = exp(33.62248566303653 - ti * 2264.475004379638);
    eqk = eg[8] * eg[43] / eg[4] / eg[56];
    rb[393] = rf[393] / max(eqk,1e-200);
    /*     R394: O + C5H4O = CO2 + C4H4 */
    rf[394] = exp(29.93360620892259 - ti * 1006.433335279839);
    eqk = eg[9] * eg[44] / eg[2] / eg[57];
    rb[394] = rf[394] / max(eqk,1e-200);
    /*     R395: HCO + c-C5H6 = CH2O + c-C5H5 */
    rf[395] = exp(alogt * 1.9 + 18.49764178508849 - ti * 8051.466682238715);
    eqk = eg[10] * eg[56] / eg[11] / eg[53];
    rb[395] = rf[395] / max(eqk,1e-200);
    /*     R396: C2H3 + c-C5H6 = CH3 + A1 */
    rf[396] = exp(155.0151385753304 - alogt * 16.08 - ti * 21366.57970799099);
    eqk = eg[16] * eg[52] / eg[22] / eg[53];
    rb[396] = rf[396] / max(eqk,1e-200);
    /*     R397: C4H5-I + c-C5H6 = C4H6 + c-C5H5 */
    rf[397] = 6e12;
    eqk = eg[43] * eg[56] / eg[48] / eg[53];
    rb[397] = rf[397] / max(eqk,1e-200);
    /*     R398: O2 + c-C5H5 = O + C5H5O(2,4) */
    rf[398] = exp(36.59033273310099 - alogt * .73 - ti * 24526.78038076968);
    eqk = eg[2] * eg[58] / eg[3] / eg[56];
    rb[398] = rf[398] / max(eqk,1e-200);
    /*     R399: O + c-C5H5 = C5H5O(2,4) */
    rf[399] = exp(alogt * 5.87 - 27.51769243062154 + ti * 8710.680516847009);
    eqk = eg[58] / eg[2] / eg[56] / pfac1;
    rb[399] = rf[399] / max(eqk,1e-200);
    /*     R400: OH + c-C5H5 = H + C5H5O(2,4) */
    rf[400] = exp(117.7393244424443 - alogt * 10.46 - ti * 28733.67172223942);
    eqk = eg[0] * eg[58] / eg[4] / eg[56];
    rb[400] = rf[400] / max(eqk,1e-200);
    /*     R401: C3H3 + C4H6 = H + C6H5CH3 */
    rf[401] = exp(alogt * 1.28 + 13.38933240825857 + ti * 2320.33205448767);
    eqk = eg[0] * eg[60] / eg[41] / eg[43];
    rb[401] = rf[401] / max(eqk,1e-200);
    /*     R402: H + C6H5CH2 = C6H5CH3 */
    rf[402] = 1e14;
    eqk = eg[60] / eg[0] / eg[59] / pfac1;
    rb[402] = rf[402] / max(eqk,1e-200);
    /*     R403: O + C6H5CH2 = H + A1CHO */
    rf[403] = 4e14;
    eqk = eg[0] * eg[62] / eg[2] / eg[59];
    rb[403] = rf[403] / max(eqk,1e-200);
    /*     R404: HO2 + C6H5CH2 = H + OH + A1CHO */
    rf[404] = 5e12;
    eqk = eg[0] * eg[4] * eg[62] / eg[6] / eg[59] * pfac1;
    rb[404] = rf[404] / max(eqk,1e-200);
    /*     R405: H + A1CHO = HCO + A1 */
    rf[405] = exp(alogt * 2.17 + 14.47303056088107 - ti * 2094.890987384986);
    eqk = eg[11] * eg[52] / eg[0] / eg[62];
    rb[405] = rf[405] / max(eqk,1e-200);
    /*     R406: O + A1 = H + C6H5O */
    rf[406] = exp(30.26290995606519 - ti * 2470.793838112006);
    eqk = eg[0] * eg[55] / eg[2] / eg[52];
    rb[406] = rf[406] / max(eqk,1e-200);
    /*     R407: H + C6H5O = CO + c-C5H6 */
    rf[407] = 5e13;
    eqk = eg[8] * eg[53] / eg[0] / eg[55];
    rb[407] = rf[407] / max(eqk,1e-200);
    /*     R408: C6H5O = CO + c-C5H5 */
    rf[408] = exp(125.6640139790803 - alogt * 12.06 - ti * 36634.17340418615);
    eqk = eg[8] * eg[56] / eg[55] * pfac1;
    rb[408] = rf[408] / max(eqk,1e-200);
    /*     R409: OH + c-C5H6 = H2O + c-C5H5 */
    rf[409] = exp(alogt * 2. + 14.94044015494976);
    eqk = eg[5] * eg[56] / eg[4] / eg[53];
    rb[409] = rf[409] / max(eqk,1e-200);
    /*     R410: O + c-C5H6 = OH + c-C5H5 */
    rf[410] = exp(alogt * 2.71 + 10.77268667687643 - ti * 556.5576344097511);
    eqk = eg[4] * eg[56] / eg[2] / eg[53];
    rb[410] = rf[410] / max(eqk,1e-200);
    /*     R411: C2H2 + C3H3 = c-C5H5 */
    rf[411] = exp(128.5693442209068 - alogt * 12.5 - ti * 21147.68045756763);
    eqk = eg[56] / eg[23] / eg[41] / pfac1;
    rb[411] = rf[411] / max(eqk,1e-200);
    /*     R412: O2 + c-C5H5 = OH + C5H4O */
    rf[412] = exp(alogt * .08 + 28.21880778083067 - ti * 9057.900017518554);
    eqk = eg[4] * eg[57] / eg[3] / eg[56];
    rb[412] = rf[412] / max(eqk,1e-200);
    /*     R413: HO2 + c-C5H5 = OH + C5H5O(2,4) */
    rf[413] = exp(68.61074405147215 - alogt * 4.69 - ti * 5862.474178005065);
    eqk = eg[4] * eg[58] / eg[6] / eg[56];
    rb[413] = rf[413] / max(eqk,1e-200);
    /*     R414: C5H5O(2,4) = H + C5H4O */
    rf[414] = exp(30.62675338948254 - ti * 15096.50002919759);
    eqk = eg[0] * eg[57] / eg[58] * pfac1;
    rb[414] = rf[414] / max(eqk,1e-200);
    /*     R415: O2 + C5H5O(2,4) = HO2 + C5H4O */
    rf[415] = 1e11;
    eqk = eg[6] * eg[57] / eg[3] / eg[58];
    rb[415] = rf[415] / max(eqk,1e-200);
    /*     R416: C5H4O = CO + 2C2H2 */
    rf[416] = exp(96.23053810480693 - alogt * 7.87 - ti * 49667.48509606007);
    eqk = eg[8] * eg[23] * eg[23] / eg[57] * pfac2;
    rb[416] = rf[416] / max(eqk,1e-200);
    /*     R417: C3H3 + C4H4 = C6H5CH2 */
    rf[417] = exp(alogt * 1.28 + 13.38933240825857 + ti * 2320.33205448767);
    eqk = eg[59] / eg[41] / eg[44] / pfac1;
    rb[417] = rf[417] / max(eqk,1e-200);
    /*     R418: H + C6H5CH3 = CH3 + A1 */
    rf[418] = exp(alogt * 2.17 + 14.47303056088107 - ti * 2094.890987384986);
    eqk = eg[16] * eg[52] / eg[0] / eg[60];
    rb[418] = rf[418] / max(eqk,1e-200);
    /*     R419: HO2 + A1- = O2 + A1 */
    rf[419] = 1e12;
    eqk = eg[3] * eg[52] / eg[6] / eg[54];
    rb[419] = rf[419] / max(eqk,1e-200);
    /*     R420: H + A1- = A1 */
    rf[420] = 1e14;
    eqk = eg[52] / eg[0] / eg[54] / pfac1;
    rb[420] = rf[420] / max(eqk,1e-200);
    /*     R421: HCO + A1- = A1CHO */
    rf[421] = 1e13;
    eqk = eg[62] / eg[11] / eg[54] / pfac1;
    rb[421] = rf[421] / max(eqk,1e-200);
    /*     R422: CH3 + A1- = C6H5CH3 */
    rf[422] = exp(30.25568970809171 - ti * 23.1479667114363);
    eqk = eg[60] / eg[16] / eg[54] / pfac1;
    rb[422] = rf[422] / max(eqk,1e-200);
    /*     R423: H + C6H5CH3 = H2 + C6H5CH2 */
    rf[423] = exp(32.46650905697885 - ti * 4206.388124802089);
    eqk = eg[1] * eg[59] / eg[0] / eg[60];
    rb[423] = rf[423] / max(eqk,1e-200);
    /*     R424: O2 + C6H5CH3 = HO2 + C6H5CH2 */
    rf[424] = exp(33.33480359058475 - ti * 21634.29097517543);
    eqk = eg[6] * eg[59] / eg[3] / eg[60];
    rb[424] = rf[424] / max(eqk,1e-200);
    /*     R425: OH + C6H5CH3 = H2O + C6H5CH2 */
    rf[425] = exp(30.41603235816689 - ti * 1393.910169362578);
    eqk = eg[5] * eg[59] / eg[4] / eg[60];
    rb[425] = rf[425] / max(eqk,1e-200);
    /*     R426: HO2 + C6H5CH3 = H2O2 + C6H5CH2 */
    rf[426] = exp(26.7084607710408 - ti * 7079.75529702603);
    eqk = eg[7] * eg[59] / eg[6] / eg[60];
    rb[426] = rf[426] / max(eqk,1e-200);
    /*     R427: CH4 + A1- = CH3 + A1 */
    rf[427] = exp(alogt * 4.57 - 5.549346121351782 - ti * 2644.906805115418);
    eqk = eg[16] * eg[52] / eg[15] / eg[54];
    rb[427] = rf[427] / max(eqk,1e-200);
    /*     R428: H + c-C5H6 = H2 + c-C5H5 */
    rf[428] = exp(alogt * 1.71 + 19.52924336347364 - ti * 2812.981172107151);
    eqk = eg[1] * eg[56] / eg[0] / eg[53];
    rb[428] = rf[428] / max(eqk,1e-200);
    /*     R429: H + C6H5CH2 = CH3 + A1- */
    rf[429] = exp(152.3760812457152 - alogt * 13.94 - ti * 32497.73239618601);
    eqk = eg[16] * eg[54] / eg[0] / eg[59];
    rb[429] = rf[429] / max(eqk,1e-200);
    /*     R430: A1- + C6H5CH3 = A1 + C6H5CH2 */
    rf[430] = exp(28.37438601264911 - ti * 2214.153337615647);
    eqk = eg[52] * eg[59] / eg[54] / eg[60];
    rb[430] = rf[430] / max(eqk,1e-200);
    /*     R431: CH3 + C6H5CH3 = CH4 + C6H5CH2 */
    rf[431] = exp(26.47900805053332 - ti * 4780.558342579237);
    eqk = eg[15] * eg[59] / eg[16] / eg[60];
    rb[431] = rf[431] / max(eqk,1e-200);
    /*     R432: H + c-C5H5 = c-C5H6 */
    rf[432] = 1e14;
    eqk = eg[53] / eg[0] / eg[56] / pfac1;
    rb[432] = rf[432] / max(eqk,1e-200);
    /*     R433: C2H + A1 = H + A1C2H */
    rf[433] = 5e13;
    eqk = eg[0] * eg[61] / eg[24] / eg[52];
    rb[433] = rf[433] / max(eqk,1e-200);
    /*     R434: C2H2 + A1- = H + A1C2H */
    /*     Reaction of PLOG type */
    rf[434] = 1.;
    eqk = eg[0] * eg[61] / eg[23] / eg[54];
    rb[434] = rf[434] / max(eqk,1e-200);
    /*     R435: OH + A1C2H = H2O + A1C2H* */
    rf[435] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[63] / eg[4] / eg[61];
    rb[435] = rf[435] / max(eqk,1e-200);
    /*     R436: H + A1C2H* = A1C2H */
    rf[436] = 1e14;
    eqk = eg[61] / eg[0] / eg[63] / pfac1;
    rb[436] = rf[436] / max(eqk,1e-200);
    /*     R437: C2H3 + A1 = H + A1C2H3 */
    rf[437] = exp(27.39529878240748 - ti * 3220.586672895486);
    eqk = eg[0] * eg[64] / eg[22] / eg[52];
    rb[437] = rf[437] / max(eqk,1e-200);
    /*     R438: C2H3 + A1- = A1C2H3 */
    /*     Reaction of PLOG type */
    rf[438] = 1.;
    eqk = eg[64] / eg[22] / eg[54] / pfac1;
    rb[438] = rf[438] / max(eqk,1e-200);
    /*     R439: C2H2 + A1C2H* = A1C2HC2H2 */
    rf[439] = exp(alogt * 2.53 + 10.04324949491129 - ti * 1439.19966945017);
    rb[439] = 0.;
    /*     R440: A1C2HC2H2 = C2H2 + A1C2H* */
    rf[440] = exp(alogt * .456 + 31.50222212683644 - ti * 21731.91500869757);
    rb[440] = 0.;
    /*     R441: A1C2HC2H2 = A1C2HC2H2u */
    rf[441] = exp(alogt * .04 + 28.62427288893883 - ti * 2562.379271622471);
    rb[441] = 0.;
    /*     R442: A1C2HC2H2u = A1C2HC2H2 */
    rf[442] = exp(alogt * .029 + 28.9894302735589 - ti * 2167.354187525134);
    rb[442] = 0.;
    /*     R443: A1C2HC2H2u = A2-1 */
    rf[443] = exp(27.8133426727225 - alogt * .003 - ti * 2256.4235376974);
    rb[443] = 0.;
    /*     R444: A2-1 = A1C2HC2H2u */
    rf[444] = exp(alogt * .6810000000000001 + 28.72294441644586 - ti *
            29988.69409133338);
    rb[444] = 0.;
    /*     R445: OH + A2 = H2O + A2-1 */
    rf[445] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[66] / eg[4] / eg[68];
    rb[445] = rf[445] / max(eqk,1e-200);
    /*     R446: H + A2-1 = A2 */
    rf[446] = 1e14;
    eqk = eg[68] / eg[0] / eg[66] / pfac1;
    rb[446] = rf[446] / max(eqk,1e-200);
    /*     R447: C2H2 + A2-1 = A2C2H2 */
    rf[447] = exp(alogt * 2.55 + 9.670293665368417 - ti * 1988.510983845907);
    eqk = eg[78] / eg[23] / eg[66] / pfac1;
    rb[447] = rf[447] / max(eqk,1e-200);
    /*     R448: C4H4 + A1- = H + A2 */
    rf[448] = exp(alogt * 2.61 + 9.441452092939569 - ti * 829.9049282717556);
    rb[448] = 0.;
    /*     R449: H + A2R5 = H2 + A2R5- */
    rf[449] = exp(33.1524820337908 - ti * 8051.466682238715);
    eqk = eg[1] * eg[71] / eg[0] / eg[69];
    rb[449] = rf[449] / max(eqk,1e-200);
    /*     R450: OH + A2R5 = H2O + A2R5- */
    rf[450] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[71] / eg[4] / eg[69];
    rb[450] = rf[450] / max(eqk,1e-200);
    /*     R451: H + A2R5- = A2R5 */
    rf[451] = 1e14;
    eqk = eg[69] / eg[0] / eg[71] / pfac1;
    rb[451] = rf[451] / max(eqk,1e-200);
    /*     R452: A2C2H2 = H + A2R5 */
    rf[452] = exp(alogt * .181 + 26.57546831672089 - ti * 7825.522398468391);
    eqk = eg[0] * eg[69] / eg[78] * pfac1;
    rb[452] = rf[452] / max(eqk,1e-200);
    /*     R453: C2H2 + A3-4 = H + A4 */
    rf[453] = exp(alogt * 2.05 + 14.0623706358958 - ti * 971.7969320461838);
    rb[453] = 0.;
    /*     R454: OH + A4 = H2O + A4-1 */
    rf[454] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[73] / eg[4] / eg[72];
    rb[454] = rf[454] / max(eqk,1e-200);
    /*     R455: H + A4-1 = A4 */
    rf[455] = 1e14;
    eqk = eg[72] / eg[0] / eg[73] / pfac1;
    rb[455] = rf[455] / max(eqk,1e-200);
    /*     R456: OH + A4 = H2O + A4-2 */
    rf[456] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[74] / eg[4] / eg[72];
    rb[456] = rf[456] / max(eqk,1e-200);
    /*     R457: H + A4-2 = A4 */
    rf[457] = 1e14;
    eqk = eg[72] / eg[0] / eg[74] / pfac1;
    rb[457] = rf[457] / max(eqk,1e-200);
    /*     R458: OH + A4 = H2O + A4-4 */
    rf[458] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[75] / eg[4] / eg[72];
    rb[458] = rf[458] / max(eqk,1e-200);
    /*     R459: H + A4-4 = A4 */
    rf[459] = 1e14;
    eqk = eg[72] / eg[0] / eg[75] / pfac1;
    rb[459] = rf[459] / max(eqk,1e-200);
    /*     R460: O + C6H5O = CO + HCO + 2C2H2 */
    rf[460] = 3e13;
    eqk = eg[8] * eg[11] * eg[23] * eg[23] / eg[2] / eg[55] * pfac2;
    rb[460] = rf[460] / max(eqk,1e-200);
    /*     R461: OH + A1C2H = CH2CO + A1- */
    rf[461] = exp(alogt * 4.5 - 8.431015495175185 + ti * 503.2166676399197);
    rb[461] = 0.;
    /*     R462: OH + A1C2H = C2H2 + C6H5O */
    rf[462] = exp(30.19597047339008 - ti * 5334.096676983148);
    rb[462] = 0.;
    /*     R463: OH + A1C2H3 = C2H4 + C6H5O */
    rf[463] = exp(30.19597047339008 - ti * 5334.096676983148);
    rb[463] = 0.;
    /*     R464: OH + A2 = H + CH2CO + A1C2H */
    rf[464] = exp(30.19597047339008 - ti * 5334.096676983148);
    rb[464] = 0.;
    /*     R465: OH + A4 = CH2CO + A3-4 */
    rf[465] = exp(30.19597047339008 - ti * 5334.096676983148);
    rb[465] = 0.;
    /*     R466: O + A1C2H = HCCO + A1- */
    rf[466] = exp(alogt * 2. + 16.83104545881444 - ti * 956.1116685158474);
    rb[466] = 0.;
    /*     R467: O + A1C2H3 = CO + CH3 + A1- */
    rf[467] = exp(alogt * 1.83 + 16.77042083699801 - ti * 110.7076668807823);
    rb[467] = 0.;
    /*     R468: O + A1C2H = C2H + C6H5O */
    rf[468] = exp(30.72206356928686 - ti * 2279.571504408836);
    rb[468] = 0.;
    /*     R469: O + A1C2H3 = C2H3 + C6H5O */
    rf[469] = exp(30.72206356928686 - ti * 2279.571504408836);
    rb[469] = 0.;
    /*     R470: O + A2 = CH2CO + A1C2H */
    rf[470] = exp(30.72206356928686 - ti * 2279.571504408836);
    rb[470] = 0.;
    /*     R471: O + A4 = HCCO + A3-4 */
    rf[471] = exp(30.72206356928686 - ti * 2279.571504408836);
    rb[471] = 0.;
    /*     R472: C2H3 + A1C2H = H + A2 */
    rf[472] = exp(40.42488042636084 - alogt * 1.4 - ti * 7929.486962002799);
    rb[472] = 0.;
    /*     R473: A1C2H = H + A1C2H* */
    rf[473] = exp(138.8970429243721 - alogt * 12.4 - ti * 74514.31127744875);
    eqk = eg[0] * eg[63] / eg[61] * pfac1;
    rb[473] = rf[473] / max(eqk,1e-200);
    /*     R474: CH3 + A1C2H4 = CH4 + A1C2H3 */
    rf[474] = exp(ti * 387.2755474156822 + 28.82796930531752);
    eqk = eg[15] * eg[64] / eg[16] / eg[76];
    rb[474] = rf[474] / max(eqk,1e-200);
    /*     R475: H + A1C2H4 = H2 + A1C2H3 */
    rf[475] = 1.67e12;
    eqk = eg[1] * eg[64] / eg[0] / eg[76];
    rb[475] = rf[475] / max(eqk,1e-200);
    /*     R476: H + A1C2H4 = A1C2H5 */
    rf[476] = 3.61e13;
    eqk = eg[77] / eg[0] / eg[76] / pfac1;
    rb[476] = rf[476] / max(eqk,1e-200);
    /*     R477: O2 + A1C2H4 = HO2 + A1C2H3 */
    rf[477] = exp(38.14969430755491 - alogt * 1.6 - ti * 1719.893926659718);
    eqk = eg[6] * eg[64] / eg[3] / eg[76];
    rb[477] = rf[477] / max(eqk,1e-200);
    /*     R478: OH + A1C2H4 = H2O + A1C2H3 */
    rf[478] = 2.41e13;
    eqk = eg[5] * eg[64] / eg[4] / eg[76];
    rb[478] = rf[478] / max(eqk,1e-200);
    /*     R479: A1C2H4 = H + A1C2H3 */
    rf[479] = exp(alogt * 2. + 15.14787657705861 - ti * 16156.0730445802);
    eqk = eg[0] * eg[64] / eg[76] * pfac1;
    rb[479] = rf[479] / max(eqk,1e-200);
    /*     R480: CH3 + A1C2H5 = CH4 + A1C2H4 */
    rf[480] = exp(alogt * 3.6 - .7940730991499059 - ti * 3599.710110295401);
    eqk = eg[15] * eg[76] / eg[16] / eg[77];
    rb[480] = rf[480] / max(eqk,1e-200);
    /*     R481: H + A1C2H5 = H2 + A1C2H4 */
    rf[481] = exp(alogt * 2.5 + 12.70684793344266 - ti * 3400.084058242645);
    eqk = eg[1] * eg[76] / eg[0] / eg[77];
    rb[481] = rf[481] / max(eqk,1e-200);
    /*     R482: H + A1C2H5 = C2H5 + A1 */
    rf[482] = exp(alogt * 2.2 + 14.65275808249798 - ti * 2095.142595718806);
    eqk = eg[20] * eg[52] / eg[0] / eg[77];
    rb[482] = rf[482] / max(eqk,1e-200);
    /*     R483: HO2 + A1C2H5 = H2O2 + A1C2H4 */
    rf[483] = exp(alogt * 2.6 + 8.480529207044645 - ti * 6999.794168538047);
    eqk = eg[7] * eg[76] / eg[6] / eg[77];
    rb[483] = rf[483] / max(eqk,1e-200);
    /*     R484: O + A1C2H5 = OH + A1C2H4 */
    rf[484] = exp(alogt * 2.7 + 11.47729828732708 - ti * 1870.204745283761);
    eqk = eg[4] * eg[76] / eg[2] / eg[77];
    rb[484] = rf[484] / max(eqk,1e-200);
    /*     R485: OH + A1C2H5 = H2O + A1C2H4 */
    rf[485] = exp(alogt * 1.5 + 19.11382792451231 - ti * 270.6299238567488);
    eqk = eg[5] * eg[76] / eg[4] / eg[77];
    rb[485] = rf[485] / max(eqk,1e-200);
    /*     R486: A1C2H5 = C2H5 + A1- */
    rf[486] = exp(67.5450759185234 - alogt * 3.6 - ti * 55436.81386888499);
    eqk = eg[20] * eg[54] / eg[77] * pfac1;
    rb[486] = rf[486] / max(eqk,1e-200);
    /*     R487: C2H4 + A1C2H* = H + A2 */
    rf[487] = exp(65.75885662967096 - alogt * 4.2 - ti * 12009.11480822639);
    rb[487] = 0.;
    /*     R488: O + A1 = OH + A1- */
    rf[488] = exp(31.31990057004249 - ti * 7396.731475972415);
    eqk = eg[4] * eg[54] / eg[2] / eg[52];
    rb[488] = rf[488] / max(eqk,1e-200);
    /*     R489: C2H4 + A1- = A1C2H4 */
    rf[489] = exp(alogt * 2. + 12.99680015442898 - ti * 970.6043085438771);
    eqk = eg[76] / eg[21] / eg[54] / pfac1;
    rb[489] = rf[489] / max(eqk,1e-200);
    /*     R490: C2H4 + A1- = C2H3 + A1 */
    rf[490] = exp(alogt * 4.5 - 4.662799298824727 - ti * 2249.076574349857);
    eqk = eg[22] * eg[52] / eg[21] / eg[54];
    rb[490] = rf[490] / max(eqk,1e-200);
    /*     R491: CH2O + A1- = H + A1CHO */
    rf[491] = exp(alogt * 2.1 + 10.27849345315958 + ti * 206.872372066771);
    eqk = eg[0] * eg[62] / eg[10] / eg[54];
    rb[491] = rf[491] / max(eqk,1e-200);
    /*     R492: CH3 + A1CHO = CO + CH4 + A1- */
    rf[492] = exp(alogt * 1.8 + 14.81614243827218 - ti * 2979.143315761853);
    rb[492] = 0.;
    /*     R493: H + A1CHO = H2 + CO + A1- */
    rf[493] = exp(alogt * 1.2 + 21.44110563009673 - ti * 1209.934155673423);
    rb[493] = 0.;
    /*     R494: HO2 + A1CHO = H2O2 + CO + A1- */
    rf[494] = exp(28.73296119468933 - ti * 6000.355544938402);
    rb[494] = 0.;
    /*     R495: O + A1CHO = OH + CO + A1- */
    rf[495] = exp(29.52111648587746 - ti * 910.4699167609067);
    rb[495] = 0.;
    /*     R496: O2 + A1CHO = HO2 + CO + A1- */
    rf[496] = exp(31.03554628768338 - ti * 19700.47964310198);
    rb[496] = 0.;
    /*     R497: OH + A1CHO = H2O + CO + A1- */
    rf[497] = exp(alogt * .7 + 23.87600185931007 + ti * 560.4827244173425);
    rb[497] = 0.;
    /*     R498: O + A1- = C6H5O */
    rf[498] = 1e14;
    eqk = eg[55] / eg[2] / eg[54] / pfac1;
    rb[498] = rf[498] / max(eqk,1e-200);
    /*     R499: C2H3 + A2 = H2 + A2C2H2 */
    rf[499] = exp(41.11802760692078 - alogt * 1.4 - ti * 7929.486962002799);
    eqk = eg[1] * eg[78] / eg[22] / eg[68];
    rb[499] = rf[499] / max(eqk,1e-200);
    /*     R500: H + A2CH3 = CH3 + A2 */
    rf[500] = exp(alogt * 2.2 + 15.16358370626397 - ti * 2095.142595718806);
    eqk = eg[16] * eg[68] / eg[0] / eg[80];
    rb[500] = rf[500] / max(eqk,1e-200);
    /*     R501: A2CH3 = CH3 + A2-1 */
    rf[501] = exp(79.45104397160324 - alogt * 5. - ti * 57493.46038952933);
    eqk = eg[16] * eg[66] / eg[80] * pfac1;
    rb[501] = rf[501] / max(eqk,1e-200);
    /*     R502: A2R5 = H + A2R5- */
    rf[502] = exp(140.3068677829022 - alogt * 12.5 - ti * 74514.31127744875);
    eqk = eg[0] * eg[71] / eg[69] * pfac1;
    rb[502] = rf[502] / max(eqk,1e-200);
    /*     R503: C2H3 + A2-1 = H + A2C2H2 */
    rf[503] = exp(alogt * 4.14 + 2.240709689275958 - ti * 11691.73605594589);
    eqk = eg[0] * eg[78] / eg[22] / eg[66];
    rb[503] = rf[503] / max(eqk,1e-200);
    /*     R504: C2H4 + A2-1 = H2 + A2C2H2 */
    rf[504] = exp(65.75885662967096 - alogt * 4.2 - ti * 12009.11480822639);
    eqk = eg[1] * eg[78] / eg[21] / eg[66];
    rb[504] = rf[504] / max(eqk,1e-200);
    /*     R505: A3-4 = C2H2 + A2R5- */
    rf[505] = exp(alogt * 1.6031 + 21.55617495988152 - ti * 31123.95089352903)
        ;
    rb[505] = 0.;
    /*     R506: C2H2 + A2R5- = A3-4 */
    rf[506] = exp(alogt * 2. + 15.00639812274156 - ti * 1591.22142474419);
    rb[506] = 0.;
    /*     R507: C2H2 + A4-2 = H + A4R5 */
    rf[507] = exp(alogt * 2.55 + 9.670293665368417 - ti * 1988.510983845907);
    eqk = eg[0] * eg[83] / eg[23] / eg[74];
    rb[507] = rf[507] / max(eqk,1e-200);
    /*     R508: C3H3 + C9H7 = H2 + A2R5 */
    rf[508] = exp(100.4221459724542 - alogt * 9.199999999999999 - ti *
            7625.242164747703);
    eqk = eg[1] * eg[69] / eg[41] / eg[79];
    rb[508] = rf[508] / max(eqk,1e-200);
    /*     R509: HO2 + C9H7 = H2O + C9H6O */
    rf[509] = exp(75.24297064405279 - alogt * 6.5 - ti * 6743.606563042564);
    eqk = eg[5] * eg[81] / eg[6] / eg[79];
    rb[509] = rf[509] / max(eqk,1e-200);
    /*     R510: HO2 + C9H7 = H + OH + C9H6O */
    rf[510] = exp(67.04499483404038 - alogt * 4.7 - ti * 5862.021283004188);
    rb[510] = 0.;
    /*     R511: O + C9H7 = H + C9H6O */
    rf[511] = 2.8e13;
    eqk = eg[0] * eg[81] / eg[2] / eg[79];
    rb[511] = rf[511] / max(eqk,1e-200);
    /*     R512: O2 + C9H7 = OH + C9H6O */
    rf[512] = exp(25.65051952210362 - ti * 12741.59698964306);
    eqk = eg[4] * eg[81] / eg[3] / eg[79];
    rb[512] = rf[512] / max(eqk,1e-200);
    /*     R513: O2 + C9H7 = H + O + C9H6O */
    rf[513] = exp(35.66694748582034 - alogt * .7 - ti * 24526.98166743674);
    rb[513] = 0.;
    /*     R514: HCO + C9H8 = CH2O + C9H7 */
    rf[514] = exp(alogt * 1.9 + 18.49764178508849 - ti * 8050.963465571075);
    eqk = eg[10] * eg[79] / eg[11] / eg[82];
    rb[514] = rf[514] / max(eqk,1e-200);
    /*     R515: HO2 + C9H8 = H2O2 + C9H7 */
    rf[515] = exp(alogt * 2.6 + 9.305650551780508 - ti * 6491.042117554088);
    eqk = eg[7] * eg[79] / eg[6] / eg[82];
    rb[515] = rf[515] / max(eqk,1e-200);
    /*     R516: O + C9H8 = OH + C9H7 */
    rf[516] = exp(alogt * 2.7 + 10.77268667687643 - ti * 1059.573015382615);
    eqk = eg[4] * eg[79] / eg[2] / eg[82];
    rb[516] = rf[516] / max(eqk,1e-200);
    /*     R517: O + C9H8 = 2H + C9H6O */
    rf[517] = exp(29.52563797059631 - alogt * .1 - ti * 181.610895351247);
    rb[517] = 0.;
    /*     R518: O2 + C9H8 = HO2 + C9H7 */
    rf[518] = exp(31.31990057004249 - ti * 24979.22248664474);
    eqk = eg[6] * eg[79] / eg[3] / eg[82];
    rb[518] = rf[518] / max(eqk,1e-200);
    /*     R519: OH + C9H8 = H2O + C9H7 */
    rf[519] = exp(alogt * 1.6 + 14.94044015494976 + ti * 18.06547836827312);
    eqk = eg[5] * eg[79] / eg[4] / eg[82];
    rb[519] = rf[519] / max(eqk,1e-200);
    /*     R520: C9H8 = H + C9H7 */
    rf[520] = exp(157.1239077321048 - alogt * 15.2 - ti * 58560.27972492597);
    eqk = eg[0] * eg[79] / eg[82] * pfac1;
    rb[520] = rf[520] / max(eqk,1e-200);
    /*     R521: C2H2 + C6H5CH2 = H + C9H8 */
    rf[521] = exp(alogt * 2.5 + 10.360912399575 - ti * 5566.18020409868);
    eqk = eg[0] * eg[82] / eg[23] / eg[59];
    rb[521] = rf[521] / max(eqk,1e-200);
    /*     R522: C6H5CH2 = C2H2 + c-C5H5 */
    rf[522] = exp(32.92933848247659 - ti * 35225.16673479438);
    eqk = eg[23] * eg[56] / eg[59] * pfac1;
    rb[522] = rf[522] / max(eqk,1e-200);
    /*     R523: O + A1C2H4 = CH2O + C6H5CH2 */
    rf[523] = 9.64e13;
    eqk = eg[10] * eg[59] / eg[2] / eg[76];
    rb[523] = rf[523] / max(eqk,1e-200);
    /*     R524: A1C2H5 = CH3 + C6H5CH2 */
    rf[524] = exp(84.75714347859332 - alogt * 5.8 - ti * 49007.1138631162);
    eqk = eg[16] * eg[59] / eg[77] * pfac1;
    rb[524] = rf[524] / max(eqk,1e-200);
    /*     R525: 2c-C5H5 = 2H + A2 */
    rf[525] = exp(68.62970196521677 - alogt * 4. - ti * 17715.99439259719);
    rb[525] = 0.;
    /*     R526: OH + c-C5H6 = HCO + C4H6 */
    rf[526] = exp(84.21481918776796 - alogt * 7.8 - ti * 3552.82541337139);
    eqk = eg[11] * eg[43] / eg[4] / eg[53];
    rb[526] = rf[526] / max(eqk,1e-200);
    /*     R527: HO2 + c-C5H5 = H2O + C5H4O */
    rf[527] = exp(76.15926137592695 - alogt * 6.52 - ti * 6743.631723875946);
    eqk = eg[5] * eg[57] / eg[6] / eg[56];
    rb[527] = rf[527] / max(eqk,1e-200);
    /*     R528: c-C5H5 + C6H5CH3 = c-C5H6 + C6H5CH2 */
    rf[528] = exp(28.37295846065793 - ti * 2214.153337615647);
    eqk = eg[53] * eg[59] / eg[56] / eg[60];
    rb[528] = rf[528] / max(eqk,1e-200);
    /*     R529: c-C5H6 + A1- = A1 + c-C5H5 */
    rf[529] = exp(26.4598381344256 - ti * 2767.691672019558);
    eqk = eg[52] * eg[56] / eg[53] / eg[54];
    rb[529] = rf[529] / max(eqk,1e-200);
    /*     R530: C2H3 + c-C5H6 = C2H4 + c-C5H5 */
    rf[530] = 6e12;
    eqk = eg[21] * eg[56] / eg[22] / eg[53];
    rb[530] = rf[530] / max(eqk,1e-200);
    /*     R531: c-C5H5 = C2H2 + C3H3 */
    rf[531] = exp(31.77574188547571 - alogt * .075 - ti * 31350.398393967);
    eqk = eg[23] * eg[41] / eg[56] * pfac1;
    rb[531] = rf[531] / max(eqk,1e-200);
    /*     R532: C2H3 + C4H6 = CH3 + c-C5H6 */
    rf[532] = exp(33.84562921435074 - ti * 12580.41669099799);
    eqk = eg[16] * eg[53] / eg[22] / eg[43];
    rb[532] = rf[532] / max(eqk,1e-200);
    /*     R533: CH2 + C6H5CH2 = H + A1C2H3 */
    rf[533] = 1e13;
    eqk = eg[0] * eg[64] / eg[17] / eg[59];
    rb[533] = rf[533] / max(eqk,1e-200);
    /*     R534: CH3 + C6H5CH2 = A1C2H5 */
    rf[534] = 2e13;
    eqk = eg[77] / eg[16] / eg[59] / pfac1;
    rb[534] = rf[534] / max(eqk,1e-200);
    /*     R535: C4H4 + A1- = C2H + A1C2H3 */
    rf[535] = exp(26.49158683274018 - ti * 956.1116685158474);
    eqk = eg[24] * eg[64] / eg[44] / eg[54];
    rb[535] = rf[535] / max(eqk,1e-200);
    /*     R536: C4H6 + A1- = C2H3 + A1C2H3 */
    rf[536] = exp(26.49158683274018 - ti * 956.1116685158474);
    eqk = eg[22] * eg[64] / eg[43] / eg[54];
    rb[536] = rf[536] / max(eqk,1e-200);
    /*     R537: 2C4H4 = A1C2H3 */
    rf[537] = exp(32.6416564100248 - ti * 19122.23337031695);
    eqk = eg[64] / eg[44] / eg[44] / pfac1;
    rb[537] = rf[537] / max(eqk,1e-200);
    /*     R538: A1C2H5 = H2 + A1C2H3 */
    rf[538] = exp(32.01304775060243 - ti * 4143.989258014739);
    eqk = eg[1] * eg[64] / eg[77] * pfac1;
    rb[538] = rf[538] / max(eqk,1e-200);
    /*     R539: A1C2H5 = H2 + A1C2H3 */
    rf[539] = exp(29.24245703102532 - ti * 32205.86672895486);
    eqk = eg[1] * eg[64] / eg[77] * pfac1;
    rb[539] = rf[539] / max(eqk,1e-200);
    /*     R540: C2H6 + A1- = H + A1C2H5 */
    rf[540] = exp(28.55130386907224 - ti * 3114.911172691103);
    eqk = eg[0] * eg[77] / eg[19] / eg[54];
    rb[540] = rf[540] / max(eqk,1e-200);
    /*     R541: C4H3-N + A1- = A2 */
    rf[541] = exp(173.1059916253803 - alogt * 17.9 - ti * 19927.38003854082);
    eqk = eg[68] / eg[46] / eg[54] / pfac1;
    rb[541] = rf[541] / max(eqk,1e-200);
    /*     R542: C4H3-N + A1- = H + A2-1 */
    rf[542] = exp(166.3958922671922 - alogt * 16.1 - ti * 29000.37655608857);
    eqk = eg[0] * eg[66] / eg[46] / eg[54];
    rb[542] = rf[542] / max(eqk,1e-200);
    /*     R543: C3H3 + A1- = C9H8 */
    rf[543] = exp(173.0993470826616 - alogt * 17.8 - ti * 19927.38003854082);
    eqk = eg[82] / eg[41] / eg[54] / pfac1;
    rb[543] = rf[543] / max(eqk,1e-200);
    /*     R544: C3H3 + A1 = H + C9H8 */
    rf[544] = exp(alogt * 2.61 + 22.55744602205842 - ti * 28431.74172165546);
    eqk = eg[0] * eg[82] / eg[41] / eg[52];
    rb[544] = rf[544] / max(eqk,1e-200);
    /*     R545: CH2 + C6H5CH3 = H2 + A1C2H3 */
    rf[545] = 3e13;
    eqk = eg[1] * eg[64] / eg[17] / eg[60];
    rb[545] = rf[545] / max(eqk,1e-200);
    /*     R546: H + C6H5CH3 = CH4 + A1- */
    rf[546] = exp(alogt * 5. - .5108256237659907 - ti * 15096.50002919759);
    eqk = eg[15] * eg[54] / eg[0] / eg[60];
    rb[546] = rf[546] / max(eqk,1e-200);
    /*     R547: O + C6H5CH3 = OH + C6H5CH2 */
    rf[547] = 6.3e11;
    eqk = eg[4] * eg[59] / eg[2] / eg[60];
    rb[547] = rf[547] / max(eqk,1e-200);
    /*     R548: C6H5CH2 + C9H7 = 2H2 + A4 */
    rf[548] = exp(28.32416829648849 - ti * 1006.433335279839);
    rb[548] = 0.;
    /*     R549: 2C9H7 = H2 + C2H2 + A4 */
    rf[549] = exp(68.62970196521677 - alogt * 4.03 - ti * 17715.99439259719);
    rb[549] = 0.;
    /*     R550: A1C2H + A1C2H* = H + A4 */
    rf[550] = exp(27.46967796551979 - ti * 2005.82163721272);
    rb[550] = 0.;
    /*     R551: CH3 + A2 = CH4 + A2-1 */
    rf[551] = exp(alogt * 3.933 - .9187938620922735 - ti * 5923.363394789494);
    rb[551] = 0.;
    /*     R552: CH4 + A2-1 = CH3 + A2 */
    rf[552] = exp(alogt * 4.248 - 3.105547139561198 - ti * 2152.257687495936);
    rb[552] = 0.;
    /*     R553: CH3 + A4 = CH4 + A4-1 */
    rf[553] = exp(alogt * 3.933 - .2256466815323282 - ti * 5923.363394789494);
    rb[553] = 0.;
    /*     R554: CH4 + A4-1 = CH3 + A4 */
    rf[554] = exp(alogt * 4.248 - 3.105547139561198 - ti * 2152.257687495936);
    rb[554] = 0.;
    /*     R555: CH3 + A4 = CH4 + A4-2 */
    rf[555] = exp(alogt * 3.933 - .2256466815323282 - ti * 5923.363394789494);
    rb[555] = 0.;
    /*     R556: CH4 + A4-2 = CH3 + A4 */
    rf[556] = exp(alogt * 4.248 - 3.105547139561198 - ti * 2152.257687495936);
    rb[556] = 0.;
    /*     R557: CH3 + A4 = CH4 + A4-4 */
    rf[557] = exp(alogt * 3.933 - .2256466815323282 - ti * 5923.363394789494);
    rb[557] = 0.;
    /*     R558: CH4 + A4-4 = CH3 + A4 */
    rf[558] = exp(alogt * 4.248 - 3.105547139561198 - ti * 2152.257687495936);
    rb[558] = 0.;
    /*     R559: H + A2 = H2 + A2-1 */
    rf[559] = exp(alogt * 1.884 + 19.316768768509 - ti * 4946.36823456659);
    rb[559] = 0.;
    /*     R560: H2 + A2-1 = H + A2 */
    rf[560] = exp(alogt * 2.467 + 10.79957557709276 - ti * 1472.613256181461);
    rb[560] = 0.;
    /*     R561: H + A4 = H2 + A4-1 */
    rf[561] = exp(alogt * 1.884 + 20.00991594906895 - ti * 4946.36823456659);
    rb[561] = 0.;
    /*     R562: H2 + A4-1 = H + A4 */
    rf[562] = exp(alogt * 2.467 + 10.79957557709276 - ti * 1472.613256181461);
    rb[562] = 0.;
    /*     R563: H + A4 = H2 + A4-2 */
    rf[563] = exp(alogt * 1.884 + 20.00991594906895 - ti * 4946.36823456659);
    rb[563] = 0.;
    /*     R564: H2 + A4-2 = H + A4 */
    rf[564] = exp(alogt * 2.467 + 10.79957557709276 - ti * 1472.613256181461);
    rb[564] = 0.;
    /*     R565: H + A4 = H2 + A4-4 */
    rf[565] = exp(alogt * 1.884 + 20.00991594906895 - ti * 4946.36823456659);
    rb[565] = 0.;
    /*     R566: H2 + A4-4 = H + A4 */
    rf[566] = exp(alogt * 2.467 + 10.79957557709276 - ti * 1472.613256181461);
    rb[566] = 0.;
    /*     R567: H + C9H8 = H2 + C9H7 */
    rf[567] = exp(alogt * 1.841 + 18.14755882283192 + ti * 63.65690845644984);
    rb[567] = 0.;
    /*     R568: H2 + C9H7 = H + C9H8 */
    rf[568] = exp(alogt * 2.378 + 10.907789161733 - ti * 13531.99940950508);
    rb[568] = 0.;
    /*     R569: CH3 + C9H8 = CH4 + C9H7 */
    rf[569] = exp(alogt * 3.614 + 1.010145307345779 - ti * 1700.922658289692);
    rb[569] = 0.;
    /*     R570: CH4 + C9H7 = CH3 + C9H8 */
    rf[570] = exp(alogt * 3.883 - 4.509860006183766 - ti * 14989.81809565793);
    rb[570] = 0.;
    /*     R571: H + A1C2H = H2 + A1C2H* */
    rf[571] = exp(alogt * 1.89 + 18.7207853364027 - ti * 5418.133860479015);
    rb[571] = 0.;
    /*     R572: H2 + A1C2H* = H + A1C2H */
    rf[572] = exp(alogt * 2.469 + 10.89302874614988 - ti * 1422.643841084817);
    rb[572] = 0.;
    /*     R573: CH3 + A1C2H = CH4 + A1C2H* */
    rf[573] = exp(alogt * 3.843 - .1266976530459575 - ti * 6166.316401926048);
    rb[573] = 0.;
    /*     R574: CH4 + A1C2H* = CH3 + A1C2H */
    rf[574] = exp(alogt * 4.154 - 1.63475572041839 - ti * 1862.354565268579);
    rb[574] = 0.;
    /*     R575: C4H4 + A4-1 = H + BAPYR */
    rf[575] = exp(alogt * 2.61 + 9.441452092939569 - ti * 829.9049282717556);
    rb[575] = 0.;
    /*     R576: C4H4 + A4-2 = H + BAPYR */
    rf[576] = exp(alogt * 2.61 + 9.441452092939569 - ti * 829.9049282717556);
    rb[576] = 0.;
    /*     R577: C4H4 + A4-4 = H + BEPYREN */
    rf[577] = exp(alogt * 2.61 + 9.441452092939569 - ti * 829.9049282717556);
    rb[577] = 0.;
    /*     R578: C2H2 + A4-1 = H + PYC2H-1 */
    rf[578] = exp(alogt * 5.71 - 11.74611935213794 - ti * 6707.878179640129);
    eqk = eg[0] * eg[86] / eg[23] / eg[73];
    rb[578] = rf[578] / max(eqk,1e-200);
    /*     R579: C2H2 + A4-2 = H + PYC2H-2 */
    rf[579] = exp(39.36709013221299 - alogt * .5600000000000001 - ti *
            11352.56802195659);
    eqk = eg[0] * eg[87] / eg[23] / eg[74];
    rb[579] = rf[579] / max(eqk,1e-200);
    /*     R580: C2H2 + A4-4 = H + PYC2H-4 */
    rf[580] = exp(alogt * 5.71 - 11.74611935213794 - ti * 6707.878179640129);
    eqk = eg[0] * eg[88] / eg[23] / eg[75];
    rb[580] = rf[580] / max(eqk,1e-200);
    /*     R581: H + PYC2H-1 = H2 + PYC2H-1JP */
    rf[581] = exp(alogt * 1.89 + 18.7207853364027 - ti * 5418.133860479015);
    rb[581] = 0.;
    /*     R582: H + PYC2H-2 = H2 + PYC2H-2JS */
    rf[582] = exp(alogt * 1.89 + 18.7207853364027 - ti * 5418.133860479015);
    rb[582] = 0.;
    /*     R583: H2 + PYC2H-1JP = H + PYC2H-1 */
    rf[583] = exp(alogt * 2.469 + 10.89302874614988 - ti * 1422.643841084817);
    rb[583] = 0.;
    /*     R584: H2 + PYC2H-2JS = H + PYC2H-2 */
    rf[584] = exp(alogt * 2.469 + 10.89302874614988 - ti * 1422.643841084817);
    rb[584] = 0.;
    /*     R585: CH3 + PYC2H-1 = CH4 + PYC2H-1JP */
    rf[585] = exp(alogt * 3.843 - .1266976530459575 - ti * 6166.316401926048);
    rb[585] = 0.;
    /*     R586: CH3 + PYC2H-2 = CH4 + PYC2H-2JS */
    rf[586] = exp(alogt * 3.843 - .1266976530459575 - ti * 6166.316401926048);
    rb[586] = 0.;
    /*     R587: CH4 + PYC2H-1JP = CH3 + PYC2H-1 */
    rf[587] = exp(alogt * 4.154 - 1.63475572041839 - ti * 1862.354565268579);
    rb[587] = 0.;
    /*     R588: CH4 + PYC2H-2JS = CH3 + PYC2H-2 */
    rf[588] = exp(alogt * 4.154 - 1.63475572041839 - ti * 1862.354565268579);
    rb[588] = 0.;
    /*     R589: OH + PYC2H-1 = H2O + PYC2H-1JP */
    rf[589] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[89] / eg[4] / eg[86];
    rb[589] = rf[589] / max(eqk,1e-200);
    /*     R590: OH + PYC2H-2 = H2O + PYC2H-2JS */
    rf[590] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[90] / eg[4] / eg[87];
    rb[590] = rf[590] / max(eqk,1e-200);
    /*     R591: C2H2 + PYC2H-1JP = BAPYRJS */
    rf[591] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[591] = 0.;
    /*     R592: C2H2 + PYC2H-2JS = BAPYRJS */
    rf[592] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[592] = 0.;
    /*     R593: H + BAPYRJS = BAPYR */
    rf[593] = 5e13;
    eqk = eg[84] / eg[0] / eg[91] / pfac1;
    rb[593] = rf[593] / max(eqk,1e-200);
    /*     R594: OH + BAPYR = H + CH2CO + PYC2H-2 */
    rf[594] = exp(29.50282329283014 - ti * 5334.096676983148);
    rb[594] = 0.;
    /*     R595: O + BAPYR = CH2CO + PYC2H-2 */
    rf[595] = exp(30.02891638872692 - ti * 2279.571504408836);
    eqk = eg[29] * eg[87] / eg[2] / eg[84];
    rb[595] = rf[595] / max(eqk,1e-200);
    /*     R596: O2 + BAPYRJS = CO + HCO + PYC2H-2 */
    rf[596] = exp(28.37295846065793 - ti * 3759.0285072702);
    rb[596] = 0.;
    /*     R597: OH + PYC2H-2 = CH2CO + A4-2 */
    rf[597] = exp(alogt * 4.5 - 8.431015495175185 + ti * 503.2166676399197);
    eqk = eg[29] * eg[74] / eg[4] / eg[87];
    rb[597] = rf[597] / max(eqk,1e-200);
    /*     R598: O + PYC2H-2 = HCCO + A4-2 */
    rf[598] = exp(alogt * 2. + 16.83104545881444 - ti * 956.1116685158474);
    eqk = eg[30] * eg[74] / eg[2] / eg[87];
    rb[598] = rf[598] / max(eqk,1e-200);
    /*     R599: H + PYC2H-4 = H2 + PYC2H-4JS */
    rf[599] = exp(alogt * 1.89 + 18.7207853364027 - ti * 5418.133860479015);
    rb[599] = 0.;
    /*     R600: H2 + PYC2H-4JS = H + PYC2H-4 */
    rf[600] = exp(alogt * 2.469 + 10.89302874614988 - ti * 1422.643841084817);
    rb[600] = 0.;
    /*     R601: CH3 + PYC2H-4 = CH4 + PYC2H-4JS */
    rf[601] = exp(alogt * 3.843 - .1266976530459575 - ti * 6166.316401926048);
    rb[601] = 0.;
    /*     R602: CH4 + PYC2H-4JS = CH3 + PYC2H-4 */
    rf[602] = exp(alogt * 4.154 - 1.63475572041839 - ti * 1862.354565268579);
    rb[602] = 0.;
    /*     R603: OH + PYC2H-4 = H2O + PYC2H-4JS */
    rf[603] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[92] / eg[4] / eg[88];
    rb[603] = rf[603] / max(eqk,1e-200);
    /*     R604: C2H2 + PYC2H-4JS = BEPYRENJS */
    rf[604] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[604] = 0.;
    /*     R605: H + BEPYRENJS = BEPYREN */
    rf[605] = 5e13;
    eqk = eg[85] / eg[0] / eg[93] / pfac1;
    rb[605] = rf[605] / max(eqk,1e-200);
    /*     R606: H + BEPYREN = H2 + BEPYRENJS */
    rf[606] = exp(alogt * 1.921 + 17.92965774748256 - ti * 4943.852151228391);
    rb[606] = 0.;
    /*     R607: H2 + BEPYRENJS = H + BEPYREN */
    rf[607] = exp(alogt * 2.498 + 10.29552964031215 - ti * 2391.084317957842);
    rb[607] = 0.;
    /*     R608: CH3 + BEPYREN = CH4 + BEPYRENJS */
    rf[608] = exp(alogt * 4.177 - 4.448166437178426 - ti * 6162.391311918456);
    rb[608] = 0.;
    /*     R609: CH4 + BEPYRENJS = CH3 + BEPYREN */
    rf[609] = exp(alogt * 4.212 - .9314043696842032 - ti * 2085.833087367467);
    rb[609] = 0.;
    /*     R610: OH + BEPYREN = H2O + BEPYRENJS */
    rf[610] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[93] / eg[4] / eg[85];
    rb[610] = rf[610] / max(eqk,1e-200);
    /*     R611: C2H2 + BEPYRENJS = H + BGHIPER */
    rf[611] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[611] = 0.;
    /*     R612: H + BAPYR = H2 + BAPYRJS */
    rf[612] = exp(alogt * 1.921 + 17.92965774748256 - ti * 4943.852151228391);
    rb[612] = 0.;
    /*     R613: H2 + BAPYRJS = H + BAPYR */
    rf[613] = exp(alogt * 2.498 + 10.29552964031215 - ti * 2391.084317957842);
    rb[613] = 0.;
    /*     R614: CH3 + BAPYR = CH4 + BAPYRJS */
    rf[614] = exp(alogt * 4.177 - 4.448166437178426 - ti * 6162.391311918456);
    rb[614] = 0.;
    /*     R615: CH4 + BAPYRJS = CH3 + BAPYR */
    rf[615] = exp(alogt * 4.212 - .9314043696842032 - ti * 2085.833087367467);
    rb[615] = 0.;
    /*     R616: OH + BAPYR = H2O + BAPYRJS */
    rf[616] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[91] / eg[4] / eg[84];
    rb[616] = rf[616] / max(eqk,1e-200);
    /*     R617: C2H2 + BAPYRJS = H + ANTHAN */
    rf[617] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[617] = 0.;
    /*     R618: OH + ANTHAN = CH2CO + BAPYRJS */
    rf[618] = exp(30.19597047339008 - ti * 5334.096676983148);
    rb[618] = 0.;
    /*     R619: O + ANTHAN = HCCO + BAPYRJS */
    rf[619] = exp(30.72206356928686 - ti * 2279.571504408836);
    rb[619] = 0.;
    /*     R620: H + BGHIPER = H2 + BGHIPEJS1 */
    rf[620] = exp(alogt * 1.921 + 17.92965774748256 - ti * 4943.852151228391);
    rb[620] = 0.;
    /*     R621: H2 + BGHIPEJS1 = H + BGHIPER */
    rf[621] = exp(alogt * 2.498 + 10.29552964031215 - ti * 2391.084317957842);
    rb[621] = 0.;
    /*     R622: CH3 + BGHIPER = CH4 + BGHIPEJS1 */
    rf[622] = exp(alogt * 4.177 - 4.448166437178426 - ti * 6162.391311918456);
    rb[622] = 0.;
    /*     R623: CH4 + BGHIPEJS1 = CH3 + BGHIPER */
    rf[623] = exp(alogt * 4.212 - .9314043696842032 - ti * 2085.833087367467);
    rb[623] = 0.;
    /*     R624: OH + BGHIPER = H2O + BGHIPEJS1 */
    rf[624] = exp(alogt * 1.42 + 18.8906843731981 - ti * 729.6641680778836);
    eqk = eg[5] * eg[96] / eg[4] / eg[94];
    rb[624] = rf[624] / max(eqk,1e-200);
    /*     R625: C2H2 + BGHIPEJS1 = H + CORONEN */
    rf[625] = exp(alogt * 1.787 + 16.74403408182482 - ti * 1641.492769841418);
    rb[625] = 0.;

    rklow[1] = exp(44.30127625414583 - alogt * 1.23);
    rklow[2] = exp(58.18788837794692 - alogt * 2.3 - ti * 24531.30933077844);
    rklow[3] = exp(56.17432494233371 - alogt * 2.3 - ti * 24531.30933077844);
    rklow[4] = exp(55.42160680152843 - alogt * 2.79 - ti * 2108.981054078904);
    rklow[5] = exp(55.56214682430743 - alogt * 2.57 - ti * 717.0837513868855);
    rklow[6] = exp(63.79313832844233 - alogt * 3.42 - ti * 42445.31948209195);
    rklow[7] = exp(58.18896018941073 - alogt * 3. - ti * 12231.68754032353);
    rklow[8] = exp(73.92173987627996 - alogt * 4.82 - ti * 3286.004839688676);
    rklow[9] = exp(99.41662410685213 - alogt * 6.995 - ti * 49311.30833870453)
        ;
    rklow[10] = exp(108.579173814992 - alogt * 8.227 - ti * 50028.34176842466)
        ;
    rklow[11] = exp(97.92940382714228 - alogt * 7.244 - ti *
            52953.64090074904);
    rklow[12] = exp(76.89235621931073 - alogt * 4.76 - ti * 1227.848669041404)
        ;
    rklow[13] = exp(63.33294832064492 - alogt * 3.14 - ti * 618.9565011971013)
        ;
    rklow[14] = exp(73.46630674524468 - alogt * 3.75 - ti * 493.9574809553452)
        ;
    rklow[15] = exp(95.09412345149228 - alogt * 7.08 - ti * 3364.003423172863)
        ;
    rklow[16] = exp(90.15077102494568 - alogt * 6.642 - ti *
            2903.056955614697);
    rklow[17] = exp(135.8820792888903 - alogt * 11.3 - ti * 48264.7686350138);
    rklow[18] = exp(133.6844662866123 - alogt * 11.3 - ti * 48264.7686350138);
    rklow[19] = exp(43.17818721905118 - alogt * .97 - ti * 7346.963347542827);
    rklow[20] = exp(68.56672716605539 - alogt * 3.8 - ti * 21851.63025392911);
    rklow[21] = exp(78.23870291760679 - alogt * 5.07 - ti * 20782.84837352868)
        ;
    rklow[22] = exp(76.97484926241725 - alogt * 5.11 - ti * 3570.32225690523);
    rklow[23] = exp(117.8479150299165 - alogt * 10.27 - ti *
            27873.17122057515);
    rklow[24] = exp(69.41402502644259 - alogt * 3.86 - ti * 1670.679336564533)
        ;
    rklow[25] = exp(117.0751647987576 - alogt * 9.310000000000001 - ti *
            50251.21643052238);
    rklow[26] = exp(73.22796257597642 - alogt * 4.664 - ti *
            1902.159003678896);
    rklow[27] = exp(77.30706390878582 - alogt * 4.8 - ti * 956.1116685158474);
    rklow[28] = exp(35.43486441946732 - alogt * .64 - ti * 25009.86838170401);
    rklow[29] = exp(172.1211809470694 - alogt * 15.74 - ti *
            49674.53012940703);
    rklow[30] = exp(135.0015492208952 - alogt * 11.94 - ti *
            4916.326199508487);
    rklow[31] = exp(138.4402845218764 - alogt * 12. - ti * 3003.096429141513);
    rklow[32] = exp(71.59524926243236 - alogt * 4.718 - ti *
            941.5183851542897);
    rklow[33] = exp(71.59524926243236 - alogt * 4.718 - ti *
            941.5183851542897);
    rklow[34] = exp(138.4915778162639 - alogt * 12.599 - ti *
            3732.358023885284);
    rklow[35] = exp(237.261574758191 - alogt * 24.63 - ti * 7341.931180866428)
        ;
    rklow[36] = exp(174.5809516235858 - alogt * 16.3 - ti * 3522.516673479438)
        ;
    rklow[37] = exp(185.6884119804479 - alogt * 18.28 - ti *
            6538.797379313116);
    rklow[38] = exp(174.5809516235858 - alogt * 16.3 - ti * 3522.516673479438)
        ;
    rklow[39] = exp(293.7633078769761 - alogt * 31.434 - ti *
            9398.074484843141);
    rklow[40] = exp(174.5809516235858 - alogt * 16.3 - ti * 3522.516673479438)
        ;
    rklow[41] = exp(141.9438303687264 - alogt * 13.55 - ti * 5715.30846355377)
        ;
    rklow[42] = exp(171.4280337665094 - alogt * 15.74 - ti *
            49677.01098757849);
    rklow[43] = exp(alogt * 8.4 + 103.6163291847321 - ti * 23902.79171289619);

    return 0;
} /* ratt_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int rdsmh_(double *t, double *smh)
{
    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double ti, tn[5], tlog;


    /* Parameter adjustments */
    --smh;

    /* Function Body */
    tlog = log(*t);
    ti = 1. / *t;

    tn[0] = tlog - 1.;
    tn[1] = *t;
    tn[2] = tn[1] * *t;
    tn[3] = tn[2] * *t;
    tn[4] = tn[3] * *t;
    /* H */
    if (*t > 1e3) {
        smh[1] = -.44668285 - ti * 25473.66 + tn[0] * 2.5;
    } else {
        smh[1] = -.44668285 - ti * 25473.66 + tn[0] * 2.5;
    }
    /* H2 */
    if (*t > 1e3) {
        smh[2] = ti * 813.065581 - 1.02432865 + tn[0] * 2.93286575 + tn[1] *
            4.13304013e-4 - tn[2] * 2.4400394e-8 + tn[3] *
            1.284170116666667e-12 - tn[4] * 3.444024e-17;
    } else {
        smh[2] = ti * 917.935173 + .6830102380000001 + tn[0] * 2.34433112 +
            tn[1] * .003990260375 - tn[2] * 3.2463585e-6 + tn[3] *
            1.67976745e-9 - tn[4] * 3.688058805e-13;
    }
    /* O */
    if (*t > 1e3) {
        smh[3] = 4.92229457 - ti * 29226.012 + tn[0] * 2.54363697 - tn[1] *
            1.36581243e-5 - tn[2] * 6.983825333333333e-10 + tn[3] *
            4.129015375e-13 - tn[4] * 2.39776847e-17;
    } else {
        smh[3] = 2.05193346 - ti * 29122.2592 + tn[0] * 3.1682671 - tn[1] *
            .00163965942 + tn[2] * 1.107177326666667e-6 - tn[3] *
            5.106721866666666e-10 + tn[4] * 1.056329855e-13;
    }
    /* O2 */
    if (*t > 1e3) {
        smh[4] = ti * 1215.97718 + 3.41536279 + tn[0] * 3.66096065 + tn[1] *
            3.281829055e-4 - tn[2] * 2.352493783333333e-8 + tn[3] *
            1.714982791666667e-12 - tn[4] * 6.495671800000001e-17;
    } else {
        smh[4] = ti * 1063.94356 + 3.65767573 + tn[0] * 3.78245636 - tn[1] *
            .00149836708 + tn[2] * 1.641217001666667e-6 - tn[3] *
            8.067745908333334e-10 + tn[4] * 1.621864185e-13;
    }
    /* OH */
    if (*t > 1e3) {
        smh[5] = 5.84494652 - ti * 3697.80808 + tn[0] * 2.83853033 + tn[1] *
            5.53706445e-4 - tn[2] * 4.900003483333333e-8 + tn[3] *
            3.505822741666666e-12 - tn[4] * 1.21144945e-16;
    } else {
        smh[5] = -.103998477 - ti * 3368.89836 + tn[0] * 3.99198424 - tn[1] *
            .001200533275 + tn[2] * 7.69440055e-7 - tn[3] *
            3.232635883333333e-10 + tn[4] * 6.8159751e-14;
    }
    /* H2O */
    if (*t > 1e3) {
        smh[6] = ti * 29885.894 + 6.88255 + tn[0] * 2.6770389 + tn[1] *
            .0014865908 - tn[2] * 1.289614816666667e-7 + tn[3] *
            7.869459500000001e-12 - tn[4] * 2.13449955e-16;
    } else {
        smh[6] = ti * 30293.726 - .84900901 + tn[0] * 4.1986352 - tn[1] *
            .00101820085 + tn[2] * 1.0867236e-6 - tn[3] *
            4.573272416666666e-10 + tn[4] * 8.859840000000001e-14;
    }
    /* HO2 */
    if (*t > 1e3) {
        smh[7] = 2.95767672 - ti * 31.0206839 + tn[0] * 4.17228741 + tn[1] *
            9.40588135e-4 - tn[2] * 5.7712881e-8 + tn[3] *
            1.622146241666667e-12 + tn[4] * 8.81284525e-18;
    } else {
        smh[7] = 3.7166622 - ti * 264.018485 + tn[0] * 4.30179807 - tn[1] *
            .002374560485 + tn[2] * 3.52638175e-6 - tn[3] *
            2.023032616666667e-9 + tn[4] * 4.646126125000001e-13;
    }
    /* H2O2 */
    if (*t > 1e3) {
        smh[8] = ti * 18007.1775 + .664970694 + tn[0] * 4.57977305 + tn[1] *
            .002026630015 - tn[2] * 2.164078833333333e-7 + tn[3] *
            1.651761666666667e-11 - tn[4] * 5.6984396e-16;
    } else {
        smh[8] = ti * 17706.7437 + 3.27373319 + tn[0] * 4.31515149 - tn[1] *
            4.23695311e-4 + tn[2] * 2.94007205e-6 - tn[3] * 1.8896912e-9
            + tn[4] * 4.54475079e-13;
    }
    /* CO */
    if (*t > 1e3) {
        smh[9] = ti * 14266.117 + 6.0170977 + tn[0] * 3.0484859 + tn[1] *
            6.7586405e-4 - tn[2] * 8.0965675e-8 + tn[3] * 6.571137e-12 -
            tn[4] * 2.3490373e-16;
    } else {
        smh[9] = ti * 14344.086 + 3.5084093 + tn[0] * 3.5795335 - tn[1] *
            3.05176845e-4 + tn[2] * 1.6946905e-7 + tn[3] *
            7.558382166666667e-11 - tn[4] * 4.52212245e-14;
    }
    /* CO2 */
    if (*t > 1e3) {
        smh[10] = ti * 49024.904 - 1.9348955 + tn[0] * 4.6365111 + tn[1] *
            .00137072845 - tn[2] * 1.659829316666667e-7 + tn[3] *
            1.3365555e-11 - tn[4] * 4.58099285e-16;
    } else {
        smh[10] = ti * 48371.971 + 9.9009035 + tn[0] * 2.356813 + tn[1] *
            .00449206495 - tn[2] * 1.187010533333333e-6 + tn[3] *
            2.047750666666667e-10 - tn[4] * 7.144273999999999e-15;
    }
    /* CH2O */
    if (*t > 1e3) {
        smh[11] = ti * 14548.6831 + 6.04207898 + tn[0] * 3.16952665 + tn[1] *
            .0030966028 - tn[2] * 3.750939433333333e-7 + tn[3] *
            3.049797166666667e-11 - tn[4] * 1.10074729e-15;
    } else {
        smh[11] = ti * 14379.1953 + .602798058 + tn[0] * 4.79372312 - tn[1] *
            .00495416661 + tn[2] * 6.220333166666666e-6 - tn[3] *
            3.160710308333333e-9 + tn[4] * 6.588632050000001e-13;
    }
    /* HCO */
    if (*t > 1e3) {
        smh[12] = 3.58077056 - ti * 3653.42928 + tn[0] * 3.92001542 + tn[1] *
            .00126139662 - tn[2] * 1.118340273333333e-7 + tn[3] *
            8.801329000000001e-12 - tn[4] * 3.718991305e-16;
    } else {
        smh[12] = 3.30834869 - ti * 3872.41185 + tn[0] * 4.2375461 - tn[1] *
            .001660376285 + tn[2] * 2.333837733333333e-6 - tn[3] *
            1.118666625e-9 + tn[4] * 2.18708104e-13;
    }
    /* CH3OH */
    if (*t > 1e3) {
        smh[13] = ti * 26002.8834 + 5.16758693 + tn[0] * 3.52726795 + tn[1] *
            .00515893915 - tn[2] * 6.048215733333333e-7 + tn[3] *
            4.8120668e-11 - tn[4] * 1.71091316e-15;
    } else {
        smh[13] = ti * 25611.9736 - .897330508 + tn[0] * 5.65851051 - tn[1] *
            .00814917095 + tn[2] * 1.15323026e-5 - tn[3] *
            6.319774383333334e-9 + tn[4] * 1.40213775e-12;
    }
    /* CH2OH */
    if (*t > 1e3) {
        smh[14] = ti * 4034.0964 - 1.84691493 + tn[0] * 5.0931437 + tn[1] *
            .0029738063 - tn[2] * 3.441624333333333e-7 + tn[3] *
            2.691734775e-11 - tn[4] * 9.406295100000002e-16;
    } else {
        smh[14] = ti * 3500.7289 + 3.309135 + tn[0] * 4.47834367 - tn[1] *
            6.7535155e-4 + tn[2] * 4.641416333333334e-6 - tn[3] *
            3.0405755e-9 + tn[4] * 7.395372499999999e-13;
    }
    /* CH3O */
    if (*t > 1e3) {
        smh[15] = -1.96680028 - ti * 378.11194 + tn[0] * 4.75779238 + tn[1] *
            .00372071237 - tn[2] * 4.495086266666667e-7 + tn[3] *
            3.6507542e-11 - tn[4] * 1.31768549e-15;
    } else {
        smh[15] = 6.57240864 - ti * 1295.6976 + tn[0] * 3.71180502 - tn[1] *
            .00140231653 + tn[2] * 6.275849516666667e-6 - tn[3] *
            3.942267408333333e-9 + tn[4] * 9.329421000000001e-13;
    }
    /* CH4 */
    if (*t > 1e3) {
        smh[16] = ti * 10009.5936 + 9.90506283 + tn[0] * 1.65326226 + tn[1] *
            .00501315495 - tn[2] * 5.5276873e-7 + tn[3] *
            4.470692816666667e-11 - tn[4] * 1.57348379e-15;
    } else {
        smh[16] = ti * 10246.5983 - 4.63848842 + tn[0] * 5.14911468 - tn[1] *
            .00683110045 + tn[2] * 8.190898683333333e-6 - tn[3] *
            4.035389725e-9 + tn[4] * 8.33017205e-13;
    }
    /* CH3 */
    if (*t > 1e3) {
        smh[17] = 4.7224799 - ti * 16509.513 + tn[0] * 2.9781206 + tn[1] *
            .002898926 - tn[2] * 3.292633333333333e-7 + tn[3] *
            2.560815833333334e-11 - tn[4] * 8.958708e-16;
    } else {
        smh[17] = 1.6735354 - ti * 16422.716 + tn[0] * 3.6571797 + tn[1] *
            .00106329895 + tn[2] * 9.097313833333333e-7 - tn[3] *
            5.515083583333334e-10 + tn[4] * 1.2328537e-13;
    }
    /* CH2 */
    if (*t > 1e3) {
        smh[18] = 4.72341711 - ti * 46041.2605 + tn[0] * 3.14631886 + tn[1] *
            .001518356295 - tn[2] * 1.660790731666667e-7 + tn[3] *
            1.254029833333333e-11 - tn[4] * 4.286677575e-16;
    } else {
        smh[18] = 1.75297945 - ti * 45872.3866 + tn[0] * 3.71757846 + tn[1] *
            6.369563e-4 + tn[2] * 3.622454183333334e-7 - tn[3] *
            2.907154166666667e-10 + tn[4] * 8.2604433e-14;
    }
    /* CH2(S) */
    if (*t > 1e3) {
        smh[19] = 4.06030621 - ti * 50504.0504 + tn[0] * 3.13501686 + tn[1] *
            .00144796963 - tn[2] * 1.361113483333333e-7 + tn[3] *
            9.464391416666667e-12 - tn[4] * 3.181314175e-16;
    } else {
        smh[19] = -.74673431 - ti * 50366.2246 + tn[0] * 4.19331325 - tn[1] *
            .00116552592 + tn[2] * 1.359460751666667e-6 - tn[3] *
            5.524883175e-10 + tn[4] * 9.661659949999999e-14;
    }
    /* C2H6 */
    if (*t > 1e3) {
        smh[20] = ti * 12447.3499 - .968698313 + tn[0] * 4.04666411 + tn[1] *
            .0076769401 - tn[2] * 9.11732475e-7 + tn[3] * 7.3152212e-11 -
            tn[4] * 2.615837655e-15;
    } else {
        smh[20] = ti * 11522.2056 + 2.66678994 + tn[0] * 4.29142572 - tn[1] *
            .002750774505 + tn[2] * 9.990640966666668e-6 - tn[3] *
            5.903887241666666e-9 + tn[4] * 1.34342918e-12;
    }
    /* C2H5 */
    if (*t > 1387.) {
        smh[21] = -8.496517709999999 - ti * 11506.5499 + tn[0] * 5.8878439 +
            tn[1] * .00515383965 - tn[2] * 5.780739933333333e-7 + tn[3] *
            4.437493808333333e-11 - tn[4] * 1.532563255e-15;
    } else {
        smh[21] = 17.1789216 - ti * 13428.4028 + tn[0] * 1.32730217 + tn[1] *
            .00883283765 - tn[2] * 1.024877596666667e-6 - tn[3] *
            2.509528883333333e-11 + tn[4] * 2.193088875e-14;
    }
    /* C2H4 */
    if (*t > 1e3) {
        smh[22] = -.269081762 - ti * 4268.65851 + tn[0] * 3.99182724 + tn[1] *
            .0052416954 - tn[2] * 6.1953557e-7 + tn[3] *
            4.955236383333334e-11 - tn[4] * 1.76815193e-15;
    } else {
        smh[22] = 4.09730213 - ti * 5089.77598 + tn[0] * 3.95920063 - tn[1] *
            .003785256865 + tn[2] * 9.516499883333334e-6 - tn[3] *
            5.763236266666667e-9 + tn[4] * 1.34942095e-12;
    }
    /* C2H3 */
    if (*t > 1e3) {
        smh[23] = 1.72812235 - ti * 33856.638 + tn[0] * 4.15026763 + tn[1] *
            .003770106705 - tn[2] * 4.38329745e-7 + tn[3] * 3.4664504e-11
            - tn[4] * 1.227037545e-15;
    } else {
        smh[23] = 7.91510092 - ti * 34474.9589 + tn[0] * 3.36377642 + tn[1] *
            1.32882861e-4 + tn[2] * 4.660345066666667e-6 - tn[3] *
            3.108224516666667e-9 + tn[4] * 7.579508800000001e-13;
    }
    /* C2H2 */
    if (*t > 1e3) {
        smh[24] = -3.99838194 - ti * 25759.4042 + tn[0] * 4.65878489 + tn[1] *
            .002441983335 - tn[2] * 2.680481466666667e-7 + tn[3] *
            2.0581212e-11 - tn[4] * 6.93029795e-16;
    } else {
        smh[24] = 13.9396761 - ti * 26428.9808 + tn[0] * .808679682 + tn[1] *
            .0116807881 - tn[2] * 5.919537233333333e-6 + tn[3] *
            2.334607983333333e-9 - tn[4] * 4.250375825000001e-13;
    }
    /* C2H */
    if (*t > 1e3) {
        smh[25] = 3.92205792 - ti * 67168.379 + tn[0] * 3.66270248 + tn[1] *
            .00191246126 - tn[2] * 2.277208333333333e-7 + tn[3] *
            1.778792e-11 - tn[4] * 6.1608424e-16;
    } else {
        smh[25] = 6.18547632 - ti * 67061.605 + tn[0] * 2.89867676 + tn[1] *
            .00664942445 - tn[2] * 4.678888783333334e-6 + tn[3] *
            2.412372958333333e-9 - tn[4] * 5.37511755e-13;
    }
    /* CH3CHO */
    if (*t > 1e3) {
        smh[26] = ti * 22593.122 - 3.4807917 + tn[0] * 5.4041108 + tn[1] *
            .0058615295 - tn[2] * 7.043856166666666e-7 + tn[3] *
            5.69770425e-11 - tn[4] * 2.04924315e-15;
    } else {
        smh[26] = ti * 21572.878 + 4.1030159 + tn[0] * 4.7294595 - tn[1] *
            .0015966429 + tn[2] * 7.922486833333334e-6 - tn[3] *
            4.788217583333333e-9 + tn[4] * 1.0965556e-12;
    }
    /* C2H3OH */
    if (*t > 1410.) {
        smh[27] = ti * 18322.1436 - 20.2080305 + tn[0] * 8.325981580000001 +
            tn[1] * .004016936405 - tn[2] * 4.39880675e-7 + tn[3] *
            3.320089383333334e-11 - tn[4] * 1.132755775e-15;
    } else {
        smh[27] = ti * 15991.4544 + 23.0438601 - tn[0] * .12797226 + tn[1] *
            .01692530365 - tn[2] * 5.510748916666667e-6 + tn[3] *
            1.373822825e-9 - tn[4] * 1.599677275e-13;
    }
    /* CH3CO */
    if (*t > 1e3) {
        smh[28] = ti * 3645.0414 - 1.6757558 + tn[0] * 5.3137165 + tn[1] *
            .00458688965 - tn[2] * 5.536731e-7 + tn[3] *
            4.495621333333334e-11 - tn[4] * 1.6226184e-15;
    } else {
        smh[28] = ti * 2682.0738 + 7.8617682 + tn[0] * 4.0358705 + tn[1] *
            4.38647435e-4 + tn[2] * 5.118335000000001e-6 - tn[3] *
            3.270630416666667e-9 + tn[4] * 7.6484345e-13;
    }
    /* CH2CHO */
    if (*t > 1e3) {
        smh[29] = ti * 1188.58659 - 8.72091393 + tn[0] * 6.53928338 + tn[1] *
            .003901193145 - tn[2] * 4.606893533333334e-7 + tn[3] *
            3.68415755e-11 - tn[4] * 1.31477145e-15;
    } else {
        smh[29] = 12.3646657 - ti * 162.944975 + tn[0] * 2.795026 + tn[1] *
            .0050549736 + tn[2] * 2.695844083333333e-6 - tn[3] *
            2.585859541666667e-9 + tn[4] * 6.97180695e-13;
    }
    /* CH2CO */
    if (*t > 1e3) {
        smh[30] = ti * 7902.94013 - 3.98525731 + tn[0] * 5.35869367 + tn[1] *
            .00347820793 - tn[2] * 4.413377283333333e-7 + tn[3] *
            3.875563266666666e-11 - tn[4] * 1.5432091e-15;
    } else {
        smh[30] = ti * 7053.94926 + 13.6079359 + tn[0] * 1.81422511 + tn[1] *
            .0099504295 - tn[2] * 3.6902668e-6 + tn[3] *
            1.208571008333333e-9 - tn[4] * 1.99438534e-13;
    }
    /* HCCO */
    if (*t > 1e3) {
        smh[31] = -5.50567269 - ti * 19359.6301 + tn[0] * 5.91479333 + tn[1] *
            .00185704365 - tn[2] * 2.168950166666666e-7 + tn[3] *
            1.720611208333333e-11 - tn[4] * 6.073837949999999e-16;
    } else {
        smh[31] = 13.696829 - ti * 20163.384 + tn[0] * 1.87607969 + tn[1] *
            .0110602709 - tn[2] * 5.981155416666667e-6 + tn[3] *
            2.545021175e-9 - tn[4] * 5.064053450000001e-13;
    }
    /* C2H4O1-2 */
    if (*t > 1e3) {
        smh[32] = ti * 9180.4251 - 7.0799605 + tn[0] * 5.4887641 + tn[1] *
            .006023095 - tn[2] * 7.222821833333333e-7 + tn[3] *
            5.835692583333333e-11 - tn[4] * 2.0974544e-15;
    } else {
        smh[32] = ti * 7560.8143 + 7.8497475 + tn[0] * 3.7590532 - tn[1] *
            .004720609 + tn[2] * 1.33849535e-5 - tn[3] *
            8.400656666666666e-9 + tn[4] * 2.00199605e-12;
    }
    /* C2H3O1-2 */
    if (*t > 1e3) {
        smh[33] = -5.47228512 - ti * 17144.6252 + tn[0] * 5.60158035 + tn[1] *
            .00458806981 - tn[2] * 5.467148366666667e-7 + tn[3] *
            4.399199066666667e-11 - tn[4] * 1.576811205e-15;
    } else {
        smh[33] = 9.59725926 - ti * 18568.1353 + tn[0] * 3.58349017 - tn[1] *
            .003011379025 + tn[2] * 1.054044778333333e-5 - tn[3] *
            6.821172558333333e-9 + tn[4] * 1.652222525e-12;
    }
    /* C2H3CHO */
    if (*t > 1393.) {
        smh[34] = ti * 14963.0281 - 30.7235061 + tn[0] * 10.4184959 + tn[1] *
            .004744816605 - tn[2] * 5.488508816666667e-7 + tn[3] *
            4.302326691666667e-11 - tn[4] * 1.507936455e-15;
    } else {
        smh[34] = ti * 11652.1584 + 22.887828 + tn[0] * .292355162 + tn[1] *
            .01771607085 - tn[2] * 4.9156054e-6 + tn[3] *
            1.067501033333333e-9 - tn[4] * 1.13072054e-13;
    }
    /* C3H8 */
    if (*t > 1e3) {
        smh[35] = ti * 16275.4066 - 13.1943379 + tn[0] * 6.6691976 + tn[1] *
            .01030543755 - tn[2] * 1.227520581666667e-6 + tn[3] *
            9.869521833333333e-11 - tn[4] * 3.53457315e-15;
    } else {
        smh[35] = ti * 14381.0883 + 5.61004451 + tn[0] * 4.21093013 + tn[1] *
            8.5443252e-4 + tn[2] * 1.177550273333333e-5 - tn[3] *
            7.667171375e-9 + tn[4] * 1.823092265e-12;
    }
    /* NC3H7 */
    if (*t > 1e3) {
        smh[36] = -8.5638971 - ti * 8859.73885 + tn[0] * 6.49636579 + tn[1] *
            .0088668996 - tn[2] * 1.041496743333333e-6 + tn[3] *
            8.294912458333333e-11 - tn[4] * 2.95099885e-15;
    } else {
        smh[36] = 8.395349189999999 - ti * 10407.4558 + tn[0] * 4.08211458 +
            tn[1] * .002616201705 + tn[2] * 8.5592411e-6 - tn[3] *
            5.827863316666667e-9 + tn[4] * 1.409097465e-12;
    }
    /* C3H6 */
    if (*t > 1e3) {
        smh[37] = ti * 741.715057 - 8.43825992 + tn[0] * 6.03870234 + tn[1] *
            .00814819655 - tn[2] * 9.702179999999999e-7 + tn[3] *
            7.799473575e-11 - tn[4] * 2.793015715e-15;
    } else {
        smh[37] = 7.53408013 - ti * 788.717123 + tn[0] * 3.83464468 + tn[1] *
            .00164539476 + tn[2] * 8.420466683333334e-6 - tn[3] *
            5.552093133333334e-9 + tn[4] * 1.318537365e-12;
    }
    /* C3H5-A */
    if (*t > 1e3) {
        smh[38] = -11.24305 - ti * 17482.449 + tn[0] * 6.5007877 + tn[1] *
            .0071623655 - tn[2] * 9.463605333333332e-7 + tn[3] *
            9.234000833333333e-11 - tn[4] * 4.518194349999999e-15;
    } else {
        smh[38] = 17.173214 - ti * 19245.629 + tn[0] * 1.3631835 + tn[1] *
            .009906910499999999 + tn[2] * 2.082843333333334e-6 - tn[3] *
            2.779629583333333e-9 + tn[4] * 7.9232855e-13;
    }
    /* C3H5-S */
    if (*t > 1e3) {
        smh[39] = -3.4186478 - ti * 29614.76 + tn[0] * 5.3725281 + tn[1] *
            .007890254500000001 - tn[2] * 9.987141666666666e-7 + tn[3] *
            7.757471999999999e-11 - tn[4] * 1.8275483e-15;
    } else {
        smh[39] = 19.989269 - ti * 30916.867 + tn[0] * .91372931 + tn[1] *
            .0132161715 - tn[2] * 1.959825e-6 - tn[3] *
            1.919639833333333e-10 + tn[4] * 1.3857744e-13;
    }
    /* C3H5-T */
    if (*t > 1e3) {
        smh[40] = -3.3527184 - ti * 27843.027 + tn[0] * 5.4255528 + tn[1] *
            .007755536 - tn[2] * 9.446391666666667e-7 + tn[3] *
            6.602032333333333e-11 - tn[4] * 8.439017000000001e-16;
    } else {
        smh[40] = 16.568878 - ti * 29040.498 + tn[0] * 1.7329209 + tn[1] *
            .01119731 - tn[2] * 8.581768500000001e-7 - tn[3] *
            5.633038833333333e-10 + tn[4] * 1.91266055e-13;
    }
    /* C3H4-P */
    if (*t > 1e3) {
        smh[41] = -8.604378499999999 - ti * 19620.942 + tn[0] * 6.02524 + tn[
            1] * .005668271 - tn[2] * 6.703898500000001e-7 + tn[3] *
                5.364671916666667e-11 - tn[4] * 1.91498175e-15;
    } else {
        smh[41] = 9.876935100000001 - ti * 20802.374 + tn[0] * 2.6803869 + tn[
            1] * .007899825500000001 + tn[2] * 4.178432666666667e-7 - tn[
                3] * 1.13813525e-9 + tn[4] * 3.30771425e-13;
    }
    /* C3H3 */
    if (*t > 1e3) {
        smh[42] = -12.584869 - ti * 39570.9594 + tn[0] * 7.14221719 + tn[1] *
            .003809511055 - tn[2] * 4.457667166666667e-7 + tn[3] *
            3.540957533333334e-11 - tn[4] * 1.257377215e-15;
    } else {
        smh[42] = 15.2058598 - ti * 40767.9941 + tn[0] * 1.35110873 + tn[1] *
            .01637056455 - tn[2] * 7.89712345e-6 + tn[3] * 3.1359185e-9 -
            tn[4] * 5.9270564e-13;
    }
    /* H2CC */
    if (*t > 1e3) {
        smh[43] = .64023701 - ti * 48316.688 + tn[0] * 4.278034 + tn[1] *
            .0023781402 - tn[2] * 2.716834833333333e-7 + tn[3] *
            2.1219005e-11 - tn[4] * 7.443189499999999e-16;
    } else {
        smh[43] = 5.920391 - ti * 48621.794 + tn[0] * 3.2815483 + tn[1] *
            .00348823955 - tn[2] * 3.975874e-7 - tn[3] *
            1.008702666666667e-10 + tn[4] * 4.90947725e-14;
    }
    /* C4H6 */
    if (*t > 1e3) {
        smh[44] = -23.328171 - ti * 9133.8516 + tn[0] * 8.867313400000001 +
            tn[1] * .007459335 - tn[2] * 5.258119333333333e-7 - tn[3] *
            3.4867775e-11 + tn[4] * 7.880629000000001e-15;
    } else {
        smh[44] = 23.089996 - ti * 11802.27 + tn[0] * .11284465 + tn[1] *
            .017184511 - tn[2] * 1.851232e-6 - tn[3] * 7.675555e-10 + tn[
            4] * 3.10325895e-13;
    }
    /* C4H4 */
    if (*t > 1e3) {
        smh[45] = -9.795211800000001 - ti * 31195.992 + tn[0] * 6.6507092 +
            tn[1] * .008064717000000001 - tn[2] * 1.19898125e-6 + tn[3] *
            1.24848225e-10 - tn[4] * 5.932055e-15;
    } else {
        smh[45] = 31.419983 - ti * 32978.504 - tn[0] * 1.9152479 + tn[1] *
            .026375439 - tn[2] * 1.194265733333333e-5 + tn[3] *
            4.589368583333333e-9 - tn[4] * 8.643114e-13;
    }
    /* C4H3-I */
    if (*t > 1e3) {
        smh[46] = -19.802597 - ti * 56600.574 + tn[0] * 9.0978165 + tn[1] *
            .00461035595 - tn[2] * 5.646406833333334e-7 + tn[3] *
            4.096708166666666e-11 - tn[4] * 7.264890000000001e-16;
    } else {
        smh[46] = 13.617462 - ti * 58005.129 + tn[0] * 2.0830412 + tn[1] *
            .020417137 - tn[2] * 1.03599475e-5 + tn[3] *
            4.306613166666666e-9 - tn[4] * 8.514591999999999e-13;
    }
    /* C4H3-N */
    if (*t > 1e3) {
        smh[47] = -1.5673981 - ti * 61600.68 + tn[0] * 5.4328279 + tn[1] *
            .0084304905 - tn[2] * 1.57188515e-6 + tn[3] * 2.14199125e-10
            - tn[4] * 1.37281545e-14;
    } else {
        smh[47] = 24.622559 - ti * 62476.199 - tn[0] * .31684113 + tn[1] *
            .02345605 - tn[2] * 1.134896833333333e-5 + tn[3] *
            4.431660083333334e-9 - tn[4] * 8.2615025e-13;
    }
    /* C4H2 */
    if (*t > 1e3) {
        smh[48] = -23.71146 - ti * 52588.039 + tn[0] * 9.1576328 + tn[1] *
            .0027715259 - tn[2] * 2.265267333333333e-7 + tn[3] *
            1.56500625e-12 + tn[4] * 1.1594768e-15;
    } else {
        smh[48] = 14.866591 - ti * 54185.211 + tn[0] * 1.0543978 + tn[1] *
            .02081348 - tn[2] * 1.097863066666667e-5 + tn[3] *
            4.438089583333333e-9 - tn[4] * 8.341581000000001e-13;
    }
    /* C4H5-I */
    if (*t > 1e3) {
        smh[49] = -28.564529 - ti * 34642.812 + tn[0] * 10.229092 + tn[1] *
            .0047425069 - tn[2] * 1.506774083333333e-8 - tn[3] *
            1.049675e-10 + tn[4] * 1.2390734e-14;
    } else {
        smh[49] = 24.394241 - ti * 37496.223 - tn[0] * .0199329 + tn[1] *
            .019002836 - tn[2] * 4.593241666666666e-6 + tn[3] *
            6.486295916666666e-10 + tn[4] * 2.01046915e-14;
    }
    /* C4H5-2 */
    if (*t > 1e3) {
        smh[50] = -45.369497 - ti * 33259.095 + tn[0] * 14.538171 - tn[1] *
            .0042838528 + tn[2] * 3.926587333333333e-6 - tn[3] *
            1.13969825e-9 + tn[4] * 1.22184635e-13;
    } else {
        smh[50] = 12.036051 - ti * 35503.316 + tn[0] * 2.969628 + tn[1] *
            .0122211225 - tn[2] * 1.520857066666667e-6 - tn[3] *
            3.538905916666667e-19 + tn[4] * 8.152364e-23;
    }
    /* C4H6-2 */
    if (*t > 1e3) {
        smh[51] = -20.985762 - ti * 14335.068 + tn[0] * 9.0338133 + tn[1] *
            .0041062255 + tn[2] * 1.1958992e-6 - tn[3] *
            4.902861166666666e-10 + tn[4] * 5.1719575e-14;
    } else {
        smh[51] = 13.529426 - ti * 15710.902 + tn[0] * 2.1373338 + tn[1] *
            .0132431145 - tn[2] * 1.509478516666667e-6 - tn[3] *
            4.615533083333334e-20 + tn[4] * 1.0640942e-23;
    }
    /* H2C4O */
    if (*t > 1e3) {
        smh[52] = -28.15985 - ti * 23469.03 + tn[0] * 10.26888 + tn[1] *
            .002448082 - tn[2] * 8.141801666666667e-8 - tn[3] *
            2.257138333333333e-11 + tn[4] * 2.5535065e-15;
    } else {
        smh[52] = 2.113424 - ti * 25458.03 + tn[0] * 4.810971 + tn[1] *
            .006569995 + tn[2] * 1.644178833333333e-7 - tn[3] *
            5.1006e-10 + tn[4] * 8.200015e-14;
    }
    /* A1 */
    if (*t > 1e3) {
        smh[53] = -40.041331 - ti * 4306.41035 + tn[0] * 11.0809576 + tn[1] *
            .0103588373 - tn[2] * 1.253576651666667e-6 + tn[3] *
            1.019341533333333e-10 - tn[4] * 3.680456395e-15;
    } else {
        smh[53] = 21.6412893 - ti * 8552.47913 + tn[0] * .504818632 + tn[1] *
            .009251032100000001 + tn[2] * 1.230576468333333e-5 - tn[3] *
            9.844645083333332e-9 + tn[4] * 2.536052145e-12;
    }
    /* c-C5H6 */
    if (*t > 1e3) {
        smh[54] = -32.209454 - ti * 11081.693 + tn[0] * 9.9757848 + tn[1] *
            .0094527715 - tn[2] * 1.140191016666667e-6 + tn[3] *
            9.24945e-11 - tn[4] * 3.3340118e-15;
    } else {
        smh[54] = 21.353453 - ti * 14801.755 + tn[0] * .86108957 + tn[1] *
            .0074020155 + tn[2] * 1.201814916666667e-5 - tn[3] *
            9.448379166666666e-9 + tn[4] * 2.4344986e-12;
    }
    /* A1- */
    if (*t > 1e3) {
        smh[55] = -35.3735134 - ti * 35559.8475 + tn[0] * 10.8444762 + tn[1] *
            .008660623649999999 - tn[2] * 1.048722081666667e-6 + tn[3] *
            8.530830083333334e-11 - tn[4] * 3.08108414e-15;
    } else {
        smh[55] = 25.2910455 - ti * 39546.8722 + tn[0] * .210306633 + tn[1] *
            .01023727535 + tn[2] * 9.8290501e-6 - tn[3] *
            8.461187916666666e-9 + tn[4] * 2.2355283e-12;
    }
    /* C6H5O */
    if (*t > 1e3) {
        smh[56] = -48.818168 - ti * 287.274751 + tn[0] * 13.722172 + tn[1] *
            .008734438549999999 - tn[2] * 1.0591742e-6 + tn[3] *
            8.624359000000001e-11 - tn[4] * 3.11705252e-15;
    } else {
        smh[56] = 27.6990274 - ti * 4778.58391 - tn[0] * .466204455 + tn[1] *
            .02067219875 + tn[2] * 2.206883183333333e-6 - tn[3] *
            4.773939741666667e-9 + tn[4] * 1.448818535e-12;
    }
    /* c-C5H5 */
    if (*t > 969.35f) {
        smh[57] = 16.0315806 - ti * 30073.0524 + tn[0] * 1.33675715 + tn[1] *
            .0162396956 - tn[2] * 2.793129566666667e-6 + tn[3] *
            3.362617808333333e-10 - tn[4] * 1.85369518e-14;
    } else {
        smh[57] = 36.7153636 - ti * 30176.9405 - tn[0] * 3.97555452 + tn[1] *
            .03706854955 - tn[2] * 1.863389083333334e-5 + tn[3] *
            7.538573133333332e-9 - tn[4] * 1.404998735e-12;
    }
    /* C5H4O */
    if (*t > 1e3) {
        smh[58] = -29.4521623 - ti * 1943.64771 + tn[0] * 10.0806824 + tn[1] *
            .008057173250000001 - tn[2] * 9.721908483333333e-7 + tn[3] *
            7.889661e-11 - tn[4] * 2.84486103e-15;
    } else {
        smh[58] = 23.5409513 - ti * 5111.59287 + tn[0] * .264576497 + tn[1] *
            .01674369135 + tn[2] * 2.795641166666667e-7 - tn[3] *
            2.468395458333334e-9 + tn[4] * 7.721573800000001e-13;
    }
    /* C5H5O(2,4) */
    if (*t > 1e3) {
        smh[59] = -20.818825 - ti * 22263.699 + tn[0] * 8.5405312 + tn[1] *
            .011494755 - tn[2] * 1.59062605e-6 + tn[3] * 1.421801e-10 -
            tn[4] * 4.872968e-15;
    } else {
        smh[59] = 39.591522 - ti * 25510.455 - tn[0] * 3.07776 + tn[1] *
            .0262908395 - tn[2] * 4.809418833333333e-6 - tn[3] *
            2.823789916666667e-10 + tn[4] * 3.16806995e-13;
    }
    /* C6H5CH2 */
    if (*t > 1e3) {
        smh[60] = -51.665589 - ti * 18564.203 + tn[0] * 14.04398 + tn[1] *
            .0117469365 - tn[2] * 1.422922783333333e-6 + tn[3] *
            1.157570083333333e-10 - tn[4] * 4.180721e-15;
    } else {
        smh[60] = 23.54882 - ti * 23307.027 + tn[0] * .4811154 + tn[1] *
            .019256416 + tn[2] * 5.476915333333333e-6 - tn[3] *
            6.414393416666667e-9 + tn[4] * 1.7711534e-12;
    }
    /* C6H5CH3 */
    if (*t > 1e3) {
        smh[61] = ti * 697.64908 - 46.728785 + tn[0] * 12.940034 + tn[1] *
            .0133456435 - tn[2] * 1.613975083333333e-6 + tn[3] *
            1.311552416666667e-10 - tn[4] * 4.73318005e-15;
    } else {
        smh[61] = 20.28221 - ti * 4075.63 + tn[0] * 1.6152663 + tn[1] *
            .010549719 + tn[2] * 1.422766966666667e-5 - tn[3] *
            1.105088833333333e-8 + tn[4] * 2.7978302e-12;
    }
    /* A1C2H */
    if (*t > 1e3) {
        smh[62] = -104.99631 - ti * 27429.445 + tn[0] * 24.090759 + tn[1] *
            3.91162e-4 + tn[2] * 1.908994e-6 - tn[3] * 5.135042e-10 + tn[
            4] * 4.667334250000001e-14;
    } else {
        smh[62] = 46.378815 - ti * 35566.242 - tn[0] * 5.2645016 + tn[1] *
            .042255521 - tn[2] * 1.2766308e-5 + tn[3] * 2.7680815e-9 - tn[
            4] * 2.38365315e-13;
    }
    /* A1CHO */
    if (*t > 1e3) {
        smh[63] = ti * 11019.744 - 47.965796 + tn[0] * 13.650737 + tn[1] *
            .0128402095 - tn[2] * 1.744454833333333e-6 + tn[3] *
            1.617785833333333e-10 - tn[4] * 6.741896e-15;
    } else {
        smh[63] = ti * 6116.9349 + 40.231735 - tn[0] * 3.1627334 + tn[1] *
            .0331846225 - tn[2] * 5.8027255e-6 - tn[3] *
            5.249948083333333e-10 + tn[4] * 4.29035505e-13;
    }
    /* A1C2H* */
    if (*t > 1e3) {
        smh[64] = -11.7512654 - ti * 64952.8135 + tn[0] * 7.23812069 + tn[1] *
            .01919060545 - tn[2] * 3.647512183333333e-6 + tn[3] *
            4.976343725e-10 - tn[4] * 3.151757335e-14;
    } else {
        smh[64] = 44.8118287 - ti * 67330.2359 - tn[0] * 4.42757639 + tn[1] *
            .04183343225 - tn[2] * 1.45017727e-5 + tn[3] * 3.919047175e-9
            - tn[4] * 5.09084925e-13;
    }
    /* A1C2H3 */
    if (*t > 1e3) {
        smh[65] = 21.4502681 - ti * 15041.317 + tn[0] * .540554217 + tn[1] *
            .0308651181 - tn[2] * 6.232455083333334e-6 + tn[3] *
            8.920548916666666e-10 - tn[4] * 5.8652492e-14;
    } else {
        smh[65] = 50.1104513 - ti * 16085.7559 - tn[0] * 5.38499941 + tn[1] *
            .04101825775 - tn[2] * 8.907697966666667e-6 + tn[3] *
            4.659125058333333e-10 + tn[4] * 2.80569525e-13;
    }
    /* A2-1 */
    if (*t > 1e3) {
        smh[67] = 5.82016697 - ti * 47840.084 + tn[0] * 3.22892303 + tn[1] *
            .0315632243 - tn[2] * 6.343039683333334e-6 + tn[3] *
            9.037839083333333e-10 - tn[4] * 5.917125599999999e-14;
    } else {
        smh[67] = 60.8902264 - ti * 50136.3344 - tn[0] * 8.027180339999999 +
            tn[1] * .051462259 - tn[2] * 1.39045335e-5 + tn[3] *
            2.267794858333334e-9 - tn[4] * 3.62279777e-14;
    }
    /* A2 */
    if (*t > 1e3) {
        smh[69] = 10.6256608 - ti * 12688.3657 + tn[0] * 1.76826275 + tn[1] *
            .0344571753 - tn[2] * 6.905369600000001e-6 + tn[3] *
            9.826192416666666e-10 - tn[4] * 6.42985305e-14;
    } else {
        smh[69] = 61.982754 - ti * 14805.9774 - tn[0] * 8.724345850000001 +
            tn[1] * .052688004 - tn[2] * 1.336184483333333e-5 + tn[3] *
            1.82121645e-9 + tn[4] * 7.103330300000001e-14;
    }
    /* A2R5 */
    if (*t > 1e3) {
        smh[70] = .723303392 - ti * 26522.3472 + tn[0] * 3.65432884 + tn[1] *
            .0376323618 - tn[2] * 7.581082516666666e-6 + tn[3] *
            1.081627841666667e-9 - tn[4] * 7.08654135e-14;
    } else {
        smh[70] = 70.2667419 - ti * 29442.6605 - tn[0] * 10.5497902 + tn[1] *
            .06276839500000001 - tn[2] * 1.727434083333333e-5 + tn[3] *
            2.941576083333333e-9 - tn[4] * 8.225419200000001e-14;
    }
    /* A2R5- */
    if (*t > 1e3) {
        smh[72] = -3.69961369 - ti * 59439.114 + tn[0] * 4.90108932 + tn[1] *
            .0349465809 - tn[2] * 7.070431e-6 + tn[3] *
            1.011219133333333e-9 - tn[4] * 6.6341855e-14;
    } else {
        smh[72] = 68.2743775 - ti * 62484.0181 - tn[0] * 9.796990170000001 +
            tn[1] * .0611386065 - tn[2] * 1.74887515e-5 + tn[3] *
            3.232886741666666e-9 - tn[4] * 1.581448305e-13;
    }
    /* A4 */
    if (*t > 1e3) {
        smh[73] = -6.19295231 - ti * 21275.589 + tn[0] * 4.54060055 + tn[1] *
            .04990576035 - tn[2] * 1.010627168333333e-5 + tn[3] *
            1.446458491666667e-9 - tn[4] * 9.495115900000001e-14;
    } else {
        smh[73] = 80.7618418 - ti * 24967.3872 - tn[0] * 13.1524443 + tn[1] *
            .0804394215 - tn[2] * 2.12866195e-5 + tn[3] *
            3.257657483333333e-9 + tn[4] * 3.719955625e-15;
    }
    /* A4-1 */
    if (*t > 1401.) {
        smh[74] = -177.494305 - ti * 40514.3093 + tn[0] * 36.3345177 + tn[1] *
            .015698401 - tn[2] * 1.817411e-6 + tn[3] *
            1.426046641666667e-10 - tn[4] * 5.00281775e-15;
    } else {
        smh[74] = 76.8451211 - ti * 55647.3533 - tn[0] * 12.0603441 + tn[1] *
            .07962377700000001 - tn[2] * 2.3593767e-5 + tn[3] *
            5.21726375e-9 - tn[4] * 5.46525805e-13;
    }
    /* A4-2 */
    if (*t > 1401.) {
        smh[75] = -178.212405 - ti * 40497.8572 + tn[0] * 36.3627635 + tn[1] *
            .0157074288 - tn[2] * 1.819840833333333e-6 + tn[3] *
            1.4286054e-10 - tn[4] * 5.013248e-15;
    } else {
        smh[75] = 76.0615466 - ti * 55577.7705 - tn[0] * 12.0714675 + tn[1] *
            .080119254 - tn[2] * 2.394946516666667e-5 + tn[3] *
            5.347921583333334e-9 - tn[4] * 5.65562605e-13;
    }
    /* A4-4 */
    if (*t > 1402.) {
        smh[76] = -177.288275 - ti * 40451.7725 + tn[0] * 36.3232136 + tn[1] *
            .01569703225 - tn[2] * 1.81673605e-6 + tn[3] *
            1.425225558333333e-10 - tn[4] * 4.999187955e-15;
    } else {
        smh[76] = 76.41929930000001 - ti * 55523.9715 - tn[0] * 11.9738229 +
            tn[1] * .0796719315 - tn[2] * 2.368199466666666e-5 + tn[3] *
            5.255467266666667e-9 - tn[4] * 5.52402935e-13;
    }
    /* A1C2H4 */
    if (*t > 1e3) {
        smh[77] = -60.0115413 - ti * 20879.1061 + tn[0] * 16.1326962 + tn[1] *
            .01414521365 - tn[2] * 1.696697933333333e-6 + tn[3] *
            1.368138641666667e-10 - tn[4] * 4.906876645e-15;
    } else {
        smh[77] = 25.0411074 - ti * 26157.2945 + tn[0] * .733299107 + tn[1] *
            .0229526579 + tn[2] * 6.304287183333333e-6 - tn[3] *
            7.603061758333333e-9 + tn[4] * 2.12794839e-12;
    }
    /* A1C2H5 */
    if (*t > 1e3) {
        smh[78] = ti * 502.4922 + 3.837099 + tn[0] * 3.878978 + tn[1] *
            .029050295 - tn[2] * 5.3273e-6 + tn[3] * 7.0408275e-10 - tn[4]
            * 4.3474125e-14;
    } else {
        smh[78] = 58.57746 - ti * 1987.29 - tn[0] * 7.266845 + tn[1] *
            .05015445 - tn[2] * 1.608619166666667e-5 + tn[3] *
            4.638256666666667e-9 - tn[4] * 7.266850000000001e-13;
    }
    /* A2C2H2 */
    if (*t > 1e3) {
        smh[79] = -19.768417 - ti * 52758.3345 + tn[0] * 8.53385239 + tn[1] *
            .03437713985 - tn[2] * 6.695835016666667e-6 + tn[3] *
            9.2900955e-10 - tn[4] * 5.95311215e-14;
    } else {
        smh[79] = 67.1406964 - ti * 56483.2554 - tn[0] * 9.26784872 + tn[1] *
            .06752156500000001 - tn[2] * 2.163197366666667e-5 + tn[3] *
            5.226836091666668e-9 - tn[4] * 5.81742885e-13;
    }
    /* C9H7 */
    if (*t > 1e3) {
        smh[80] = 2.57680216 - ti * 30684.3457 + tn[0] * 3.65597547 + tn[1] *
            .02874042315 - tn[2] * 5.714510000000001e-6 + tn[3] *
            8.085656608333334e-10 - tn[4] * 5.2693206e-14;
    } else {
        smh[80] = 62.8218291 - ti * 33164.1009 - tn[0] * 8.73685384 + tn[1] *
            .051710818 - tn[2] * 1.539038988333333e-5 + tn[3] *
            3.130191316666667e-9 - tn[4] * 2.20302635e-13;
    }
    /* A2CH3 */
    if (*t > 1e3) {
        smh[81] = 6.01699287 - ti * 9791.041499999999 + tn[0] * 3.05206674 +
            tn[1] * .0384085034 - tn[2] * 7.538687066666666e-6 + tn[3] *
            1.056115308333333e-9 - tn[4] * 6.832961549999999e-14;
    } else {
        smh[81] = 61.6822357 - ti * 12078.5125 - tn[0] * 8.360184009999999 +
            tn[1] * .0588522945 - tn[2] * 1.540497628333333e-5 + tn[3] *
            2.504592666666667e-9 - tn[4] * 5.1295817e-14;
    }
    /* C9H6O */
    if (*t > 1e3) {
        smh[82] = -.703868477 - ti * 4578.5714 + tn[0] * 4.65659248 + tn[1] *
            .0285027911 - tn[2] * 5.719569983333333e-6 + tn[3] *
            8.134812016666667e-10 - tn[4] * 5.31670185e-14;
    } else {
        smh[82] = 54.0945996 - ti * 6888.83578 - tn[0] * 6.53928778 + tn[1] *
            .0484661643 - tn[2] * 1.362831093333333e-5 + tn[3] *
            2.472495616666667e-9 - tn[4] * 1.12496696e-13;
    }
    /* C9H8 */
    if (*t > 1e3) {
        smh[83] = 15.3482362 - ti * 16816.6108 + tn[0] * 1.15459802 + tn[1] *
            .0327112098 - tn[2] * 6.541751783333334e-6 + tn[3] *
            9.297403416666667e-10 - tn[4] * 6.07963375e-14;
    } else {
        smh[83] = 60.6774763 - ti * 18658.9996 - tn[0] * 8.12447817 + tn[1] *
            .04888285335 - tn[2] * 1.21739329e-5 + tn[3] *
            1.569125083333333e-9 + tn[4] * 9.20166065e-14;
    }
    /* A4R5 */
    if (*t > 1e3) {
        smh[84] = -14.2387586 - ti * 33443.9422 + tn[0] * 6.20190827 + tn[1] *
            .0533281375 - tn[2] * 1.084715501666667e-5 + tn[3] *
            1.55646775e-9 - tn[4] * 1.02323675e-13;
    } else {
        smh[84] = 88.818403 - ti * 37846.7972 - tn[0] * 14.7695663 + tn[1] *
            .089829011 - tn[2] * 2.469831666666667e-5 + tn[3] *
            4.147684183333334e-9 - tn[4] * 1.032186675e-13;
    }
    /* BAPYR */
    if (*t > 1396.) {
        smh[85] = -219.04074 - ti * 16316.1947 + tn[0] * 43.9390181 + tn[1] *
            .02146427815 - tn[2] * 2.486653983333333e-6 + tn[3] *
            1.9518595e-10 - tn[4] * 6.8487222e-15;
    } else {
        smh[85] = 77.3099666 - ti * 34529.8004 - tn[0] * 11.9841377 + tn[1] *
            .09205345550000001 - tn[2] * 2.5460174e-5 + tn[3] *
            5.330713008333333e-9 - tn[4] * 5.35981215e-13;
    }
    /* BEPYREN */
    if (*t > 1396.) {
        smh[86] = -218.993441 - ti * 15218.8378 + tn[0] * 43.9081229 + tn[1] *
            .0214622456 - tn[2] * 2.4853468e-6 + tn[3] *
            1.950297416666667e-10 - tn[4] * 6.841981250000001e-15;
    } else {
        smh[86] = 76.9031762 - ti * 33411.1148 - tn[0] * 11.9222172 + tn[1] *
            .091873361 - tn[2] * 2.537338966666667e-5 + tn[3] *
            5.304148716666666e-9 - tn[4] * 5.32458555e-13;
    }
    /* PYC2H-1 */
    if (*t > 1380.) {
        smh[87] = -195.526559 - ti * 35449.181 + tn[0] * 40.0549856 + tn[1] *
            .0185390285 - tn[2] * 2.18869815e-6 + tn[3] *
            1.739284108333334e-10 - tn[4] * 6.15448435e-15;
    } else {
        smh[87] = 30.3979898 - ti * 50437.0421 - tn[0] * 1.66684897 + tn[1] *
            .0645734905 - tn[2] * 1.480535228333333e-5 + tn[3] *
            2.470921766666667e-9 - tn[4] * 1.938004235e-13;
    }
    /* PYC2H-2 */
    if (*t > 1380.) {
        smh[88] = -196.221074 - ti * 35449.181 + tn[0] * 40.0549856 + tn[1] *
            .0185390285 - tn[2] * 2.18869815e-6 + tn[3] *
            1.739284108333334e-10 - tn[4] * 6.15448435e-15;
    } else {
        smh[88] = 29.7034755 - ti * 50437.0421 - tn[0] * 1.66684897 + tn[1] *
            .0645734905 - tn[2] * 1.480535228333333e-5 + tn[3] *
            2.470921766666667e-9 - tn[4] * 1.938004235e-13;
    }
    /* PYC2H-4 */
    if (*t > 1380.) {
        smh[89] = -195.526559 - ti * 35449.181 + tn[0] * 40.0549856 + tn[1] *
            .0185390285 - tn[2] * 2.18869815e-6 + tn[3] *
            1.739284108333334e-10 - tn[4] * 6.15448435e-15;
    } else {
        smh[89] = 30.3979898 - ti * 50437.0421 - tn[0] * 1.66684897 + tn[1] *
            .0645734905 - tn[2] * 1.480535228333333e-5 + tn[3] *
            2.470921766666667e-9 - tn[4] * 1.938004235e-13;
    }
    /* PYC2H-1JP */
    if (*t > 1376.) {
        smh[90] = -191.256136 - ti * 66776.4045 + tn[0] * 39.5503456 + tn[1] *
            .01754904365 - tn[2] * 2.080250466666667e-6 + tn[3] *
            1.657542983333333e-10 - tn[4] * 5.876126600000001e-15;
    } else {
        smh[90] = 25.9457844 - ti * 81324.6792 - tn[0] * .446464976 + tn[1] *
            .0608389695 - tn[2] * 1.360838016666667e-5 + tn[3] *
            2.178414525e-9 - tn[4] * 1.606582045e-13;
    }
    /* PYC2H-2JS */
    if (*t > 1376.) {
        smh[91] = -191.256136 - ti * 66021.4976 + tn[0] * 39.5503456 + tn[1] *
            .01754904365 - tn[2] * 2.080250466666667e-6 + tn[3] *
            1.657542983333333e-10 - tn[4] * 5.876126600000001e-15;
    } else {
        smh[91] = 25.9457844 - ti * 80569.7723 - tn[0] * .446464976 + tn[1] *
            .0608389695 - tn[2] * 1.360838016666667e-5 + tn[3] *
            2.178414525e-9 - tn[4] * 1.606582045e-13;
    }
    /* BAPYRJS */
    if (*t > 1395.) {
        smh[92] = -212.799296 - ti * 47023.4412 + tn[0] * 43.0858151 + tn[1] *
            .02072365285 - tn[2] * 2.413771e-6 + tn[3] *
            1.901235516666667e-10 - tn[4] * 6.686803049999999e-15;
    } else {
        smh[92] = 71.1484565 - ti * 64606.5317 - tn[0] * 10.3965058 + tn[1] *
            .0875479515 - tn[2] * 2.39467615e-5 + tn[3] *
            4.961913991666666e-9 - tn[4] * 4.94965696e-13;
    }
    /* PYC2H-4JS */
    if (*t > 1376.) {
        smh[93] = -191.256136 - ti * 66021.4976 + tn[0] * 39.5503456 + tn[1] *
            .01754904365 - tn[2] * 2.080250466666667e-6 + tn[3] *
            1.657542983333333e-10 - tn[4] * 5.876126600000001e-15;
    } else {
        smh[93] = 25.9457844 - ti * 80569.7723 - tn[0] * .446464976 + tn[1] *
            .0608389695 - tn[2] * 1.360838016666667e-5 + tn[3] *
            2.178414525e-9 - tn[4] * 1.606582045e-13;
    }
    /* BEPYRENJS */
    if (*t > 1395.) {
        smh[94] = -212.695279 - ti * 45929.3949 + tn[0] * 43.0444483 + tn[1] *
            .02072879015 - tn[2] * 2.41345485e-6 + tn[3] *
            1.900520475e-10 - tn[4] * 6.683175600000001e-15;
    } else {
        smh[94] = 70.8784003 - ti * 63492.907 - tn[0] * 10.362817 + tn[1] *
            .087416878 - tn[2] * 2.388148383333334e-5 + tn[3] *
            4.941508708333334e-9 - tn[4] * 4.92205212e-13;
    }
    /* BGHIPER */
    if (*t > 1369.) {
        smh[95] = -245.189878 - ti * 12295.1825 + tn[0] * 48.125628 + tn[1] *
            .02318756935 - tn[2] * 2.756869133333334e-6 + tn[3] *
            2.200976508333333e-10 - tn[4] * 7.813143049999999e-15;
    } else {
        smh[95] = 25.6313149 - ti * 30919.5278 - tn[0] * 1.3673045 + tn[1] *
            .074071659 - tn[2] * 1.527899893333333e-5 + tn[3] *
            2.148452666666667e-9 - tn[4] * 1.27219431e-13;
    }
    /* BGHIPEJS1 */
    if (*t > 1365.) {
        smh[97] = -239.621324 - ti * 43100.1224 + tn[0] * 47.4614538 + tn[1] *
            .02217007545 - tn[2] * 2.637896616666667e-6 + tn[3] *
            2.1070548e-10 - tn[4] * 7.482352500000001e-15;
    } else {
        smh[97] = 13.2343772 - ti * 60758.9096 + tn[0] * 1.49345169 + tn[1] *
            .067513695 - tn[2] * 1.289483925e-5 + tn[3] *
            1.514787783333333e-9 - tn[4] * 5.0036207e-14;
    }

    return 0;
} /* rdsmh_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int ratx_(double *t, double *c__, double *rf,
        double *rb, double *rklow)
{
    /* System generated locals */
    double d__1;

    /* Builtin functions */
    //    double log(double), log10(double), exp(double), pow(
    //	    double *, double *);

    /* Local variables */
    static long int k;
    static double p, fc, pl, rl, pr, xn, p1l, r1l, r2l, ctb, dpl, flog,
                  pcor, ctot, fclog, fcent, alogt, prlog, cprlog;


    /* Parameter adjustments */
    --rklow;
    --rb;
    --rf;
    --c__;

    /* Function Body */
    alogt = log(*t);
    ctot = 0.f;
    for (k = 1; k <= 86; ++k) {
        ctot += c__[k];
    }

    p = ctot * 83145100. * *t;
    pl = log(p);
    /*     R1: H + O2 = O + OH */
    rf[1] = rf[1] * c__[1] * c__[4];
    rb[1] = rb[1] * c__[3] * c__[5];
    /*     R2: H2 + O = H + OH */
    rf[2] = rf[2] * c__[2] * c__[3];
    rb[2] = rb[2] * c__[1] * c__[5];
    /*     R3: H2 + OH = H + H2O */
    rf[3] = rf[3] * c__[2] * c__[5];
    rb[3] = rb[3] * c__[1] * c__[6];
    /*     R4: O + H2O = 2OH */
    rf[4] = rf[4] * c__[3] * c__[6];
    rb[4] = rb[4] * c__[5] * c__[5];
    /*     R5: H2 = 2H */
    ctb = ctot + c__[2] * 1.5 + c__[6] * 11. + c__[9] * .8999999999999999 +
        c__[10] * 2.8 + c__[14] + c__[17] * 2.;
    rf[5] = rf[5] * ctb * c__[2];
    rb[5] = rb[5] * ctb * c__[1] * c__[1];
    /*     R6: 2O = O2 */
    ctb = ctot + c__[2] * 1.5 + c__[6] * 11. + c__[9] * .8999999999999999 +
        c__[10] * 2.8 + c__[14] + c__[17] * 2.;
    rf[6] = rf[6] * ctb * c__[3] * c__[3];
    rb[6] = rb[6] * ctb * c__[4];
    /*     R7: H + O = OH */
    ctb = ctot + c__[2] * 1.5 + c__[6] * 11. + c__[9] * .5 + c__[10] + c__[14]
        + c__[17] * 2.;
    rf[7] = rf[7] * ctb * c__[1] * c__[3];
    rb[7] = rb[7] * ctb * c__[5];
    /*     R8: H + OH = H2O */
    ctb = ctot - c__[2] * .27 + c__[6] * 2.65 + c__[14] + c__[17] * 2.;
    rf[8] = rf[8] * ctb * c__[1] * c__[5];
    rb[8] = rb[8] * ctb * c__[6];
    /*     R9: H + O2 = HO2 */
    ctb = ctot + c__[2] * .3 + c__[9] * .8999999999999999 + c__[10] * 2.8 +
        c__[6] * 9. + c__[14] + c__[17] * 2.;
    pr = rklow[1] * ctb / rf[9];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 1e-30) * .33 + exp(-(*t) / 1e30) * .67 + exp(-1e30 / *
            t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[9] *= pcor;
    rb[9] *= pcor;
    rf[9] = rf[9] * c__[1] * c__[4];
    rb[9] *= c__[7];
    /*     R10: H + HO2 = 2OH */
    rf[10] = rf[10] * c__[1] * c__[7];
    rb[10] = rb[10] * c__[5] * c__[5];
    /*     R11: H2 + O2 = H + HO2 */
    rf[11] = rf[11] * c__[2] * c__[4];
    rb[11] = rb[11] * c__[1] * c__[7];
    /*     R12: O + HO2 = O2 + OH */
    rf[12] = rf[12] * c__[3] * c__[7];
    rb[12] = rb[12] * c__[4] * c__[5];
    /*     R13: OH + HO2 = O2 + H2O */
    rf[13] = rf[13] * c__[5] * c__[7];
    rb[13] = rb[13] * c__[4] * c__[6];
    /*     R14: 2HO2 = O2 + H2O2 */
    rf[14] = rf[14] * c__[7] * c__[7];
    rb[14] = rb[14] * c__[4] * c__[8];
    /*     R15: 2HO2 = O2 + H2O2 */
    rf[15] = rf[15] * c__[7] * c__[7];
    rb[15] = rb[15] * c__[4] * c__[8];
    /*     R16: H2O2 = 2OH */
    pr = rklow[2] * c__[6] / rf[16];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fc = (prlog - .2040720179856174) / (1.121385876355621 - (prlog -
                .2040720179856174) * .14);
    fc = exp(-.6733445532637656 / (fc * fc + 1.));
    pcor = fc * pcor;
    rf[16] *= pcor;
    rb[16] *= pcor;
    rf[16] *= c__[8];
    rb[16] = rb[16] * c__[5] * c__[5];
    /*     R17: H2O2 = 2OH */
    ctb = ctot - c__[6] + c__[10] * .6000000000000001 + c__[86] * .5 + c__[4]
        * .2 + c__[8] * 6.7 + c__[2] * 2.7 + c__[9] * 1.8;
    pr = rklow[3] * ctb / rf[17];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fc = (prlog - .154423865238323) / (1.215495061413925 - (prlog -
                .154423865238323) * .14);
    fc = exp(-.843970070294529 / (fc * fc + 1.));
    pcor = fc * pcor;
    rf[17] *= pcor;
    rb[17] *= pcor;
    rf[17] *= c__[8];
    rb[17] = rb[17] * c__[5] * c__[5];
    /*     R18: H + H2O2 = OH + H2O */
    rf[18] = rf[18] * c__[1] * c__[8];
    rb[18] = rb[18] * c__[5] * c__[6];
    /*     R19: H + H2O2 = H2 + HO2 */
    rf[19] = rf[19] * c__[1] * c__[8];
    rb[19] = rb[19] * c__[2] * c__[7];
    /*     R20: O + H2O2 = OH + HO2 */
    rf[20] = rf[20] * c__[3] * c__[8];
    rb[20] = rb[20] * c__[5] * c__[7];
    /*     R21: OH + H2O2 = H2O + HO2 */
    rf[21] = rf[21] * c__[5] * c__[8];
    rb[21] = rb[21] * c__[6] * c__[7];
    /*     R22: OH + H2O2 = H2O + HO2 */
    rf[22] = rf[22] * c__[5] * c__[8];
    rb[22] = rb[22] * c__[6] * c__[7];
    /*     R23: O + CO = CO2 */
    ctb = ctot + c__[2] + c__[6] * 11. + c__[9] * .75 + c__[10] * 2.6;
    pr = rklow[4] * ctb / rf[23];
    pcor = pr / (pr + 1.f);
    rf[23] *= pcor;
    rb[23] *= pcor;
    rf[23] = rf[23] * c__[3] * c__[9];
    rb[23] *= c__[10];
    /*     R24: O2 + CO = O + CO2 */
    rf[24] = rf[24] * c__[4] * c__[9];
    rb[24] = rb[24] * c__[3] * c__[10];
    /*     R25: OH + CO = H + CO2 */
    rf[25] = rf[25] * c__[5] * c__[9];
    rb[25] = rb[25] * c__[1] * c__[10];
    /*     R26: OH + CO = H + CO2 */
    rf[26] = rf[26] * c__[5] * c__[9];
    rb[26] = rb[26] * c__[1] * c__[10];
    /*     R27: HO2 + CO = OH + CO2 */
    rf[27] = rf[27] * c__[7] * c__[9];
    rb[27] = rb[27] * c__[5] * c__[10];
    /*     R28: HCO = H + CO */
    ctb = ctot + c__[2] + c__[6] * 11. + c__[9] * .5 + c__[10] + c__[14] +
        c__[17] * 2.;
    rf[28] = rf[28] * ctb * c__[12];
    rb[28] = rb[28] * ctb * c__[1] * c__[9];
    /*     R29: O2 + HCO = HO2 + CO */
    rf[29] = rf[29] * c__[4] * c__[12];
    rb[29] = rb[29] * c__[7] * c__[9];
    /*     R30: H + HCO = H2 + CO */
    rf[30] = rf[30] * c__[1] * c__[12];
    rb[30] = rb[30] * c__[2] * c__[9];
    /*     R31: O + HCO = OH + CO */
    rf[31] = rf[31] * c__[3] * c__[12];
    rb[31] = rb[31] * c__[5] * c__[9];
    /*     R32: O + HCO = H + CO2 */
    rf[32] = rf[32] * c__[3] * c__[12];
    rb[32] = rb[32] * c__[1] * c__[10];
    /*     R33: OH + HCO = H2O + CO */
    rf[33] = rf[33] * c__[5] * c__[12];
    rb[33] = rb[33] * c__[6] * c__[9];
    /*     R34: HO2 + HCO = H + OH + CO2 */
    rf[34] = rf[34] * c__[7] * c__[12];
    /*     R35: 2HCO = H2 + 2CO */
    rf[35] = rf[35] * c__[12] * c__[12];
    /*     R36: HCO + CH3 = CO + CH4 */
    rf[36] = rf[36] * c__[12] * c__[15];
    rb[36] = rb[36] * c__[9] * c__[14];
    /*     R37: O2 + CH2O = HO2 + HCO */
    rf[37] = rf[37] * c__[4] * c__[11];
    rb[37] = rb[37] * c__[7] * c__[12];
    /*     R38: 2HCO = CO + CH2O */
    rf[38] = rf[38] * c__[12] * c__[12];
    rb[38] = rb[38] * c__[9] * c__[11];
    /*     R39: H + HCO = CH2O */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[5] * ctb / rf[39];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 271.) * .2176 + exp(-(*t) / 2755.) * .7824 + exp(
            -6570. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[39] *= pcor;
    rb[39] *= pcor;
    rf[39] = rf[39] * c__[1] * c__[12];
    rb[39] *= c__[11];
    /*     R40: H2 + CO = CH2O */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[6] * ctb / rf[40];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 197.) * .06799999999999995 + exp(-(*t) / 1540.) *
        .9320000000000001 + exp(-10300. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[40] *= pcor;
    rb[40] *= pcor;
    rf[40] = rf[40] * c__[2] * c__[9];
    rb[40] *= c__[11];
    /*     R41: OH + CH2O = H2O + HCO */
    rf[41] = rf[41] * c__[5] * c__[11];
    rb[41] = rb[41] * c__[6] * c__[12];
    /*     R42: H + CH2O = H2 + HCO */
    rf[42] = rf[42] * c__[1] * c__[11];
    rb[42] = rb[42] * c__[2] * c__[12];
    /*     R43: O + CH2O = OH + HCO */
    rf[43] = rf[43] * c__[3] * c__[11];
    rb[43] = rb[43] * c__[5] * c__[12];
    /*     R44: CH2O + CH3 = HCO + CH4 */
    rf[44] = rf[44] * c__[11] * c__[15];
    rb[44] = rb[44] * c__[12] * c__[14];
    /*     R45: HO2 + CH2O = H2O2 + HCO */
    rf[45] = rf[45] * c__[7] * c__[11];
    rb[45] = rb[45] * c__[8] * c__[12];
    /*     R46: CH3O = H + CH2O */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[7] * ctb / rf[46];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 2500.) * .09999999999999998 + exp(-(*t) / 1300.) * .9
        + exp(-1e99 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[46] *= pcor;
    rb[46] *= pcor;
    rb[46] = rb[46] * c__[1] * c__[11];
    /*     R47: O2 + CH3O = HO2 + CH2O */
    rf[47] *= c__[4];
    rb[47] = rb[47] * c__[7] * c__[11];
    /*     R48: CH2O + CH3O = HCO + CH3OH */
    rf[48] *= c__[11];
    rb[48] = rb[48] * c__[12] * c__[13];
    /*     R49: CH3OH + CH3 = CH3O + CH4 */
    rf[49] = rf[49] * c__[13] * c__[15];
    rb[49] *= c__[14];
    /*     R50: CH3O + CH3 = CH2O + CH4 */
    rf[50] *= c__[15];
    rb[50] = rb[50] * c__[11] * c__[14];
    /*     R51: H + CH3O = H2 + CH2O */
    rf[51] *= c__[1];
    rb[51] = rb[51] * c__[2] * c__[11];
    /*     R52: HO2 + CH3O = H2O2 + CH2O */
    rf[52] *= c__[7];
    rb[52] = rb[52] * c__[8] * c__[11];
    /*     R53: H + CH2O = CH2OH */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[8] * ctb / rf[53];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 103.) * .2813 + exp(-(*t) / 1291.) * .7187 + exp(
            -4160. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[53] *= pcor;
    rb[53] *= pcor;
    rf[53] = rf[53] * c__[1] * c__[11];
    /*     R54: O2 + CH2OH = HO2 + CH2O */
    rf[54] *= c__[4];
    rb[54] = rb[54] * c__[7] * c__[11];
    /*     R55: O2 + CH2OH = HO2 + CH2O */
    rf[55] *= c__[4];
    rb[55] = rb[55] * c__[7] * c__[11];
    /*     R56: H + CH2OH = H2 + CH2O */
    rf[56] *= c__[1];
    rb[56] = rb[56] * c__[2] * c__[11];
    /*     R57: HO2 + CH2OH = H2O2 + CH2O */
    rf[57] *= c__[7];
    rb[57] = rb[57] * c__[8] * c__[11];
    /*     R58: HCO + CH2OH = 2CH2O */
    rf[58] *= c__[12];
    rb[58] = rb[58] * c__[11] * c__[11];
    /*     R59: CH2OH + CH3O = CH2O + CH3OH */
    rb[59] = rb[59] * c__[11] * c__[13];
    /*     R60: HCO + CH3OH = CH2O + CH2OH */
    rf[60] = rf[60] * c__[12] * c__[13];
    rb[60] *= c__[11];
    /*     R61: OH + CH2OH = H2O + CH2O */
    rf[61] *= c__[5];
    rb[61] = rb[61] * c__[6] * c__[11];
    /*     R62: O + CH2OH = OH + CH2O */
    rf[62] *= c__[3];
    rb[62] = rb[62] * c__[5] * c__[11];
    /*     R63: 2CH2OH = CH2O + CH3OH */
    rb[63] = rb[63] * c__[11] * c__[13];
    /*     R64: CH3OH = OH + CH3 */
    ctb = ctot;
    pr = rklow[9] * ctb / rf[64];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 35580.) * 1.4748 - exp(-(*t) / 1116.) * .4748 + exp(
            -9023. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[64] *= pcor;
    rb[64] *= pcor;
    rf[64] *= c__[13];
    rb[64] = rb[64] * c__[5] * c__[15];
    /*     R65: CH3OH = H2O + CH2(S) */
    ctb = ctot;
    pr = rklow[10] * ctb / rf[65];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 3290.) * -1.545 + exp(-(*t) / 47320.) * 2.545 + exp(
            -47110. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[65] *= pcor;
    rb[65] *= pcor;
    rf[65] *= c__[13];
    rb[65] *= c__[6];
    /*     R66: CH3OH = H + CH2OH */
    ctb = ctot;
    pr = rklow[11] * ctb / rf[66];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 37050.) * 74.91 - exp(-(*t) / 41500.) * 73.91 + exp(
            -5220. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[66] *= pcor;
    rb[66] *= pcor;
    rf[66] *= c__[13];
    rb[66] *= c__[1];
    /*     R67: H + CH3OH = H2 + CH2OH */
    rf[67] = rf[67] * c__[1] * c__[13];
    rb[67] *= c__[2];
    /*     R68: H + CH3OH = H2 + CH3O */
    rf[68] = rf[68] * c__[1] * c__[13];
    rb[68] *= c__[2];
    /*     R69: O + CH3OH = OH + CH2OH */
    rf[69] = rf[69] * c__[3] * c__[13];
    rb[69] *= c__[5];
    /*     R70: OH + CH3OH = H2O + CH2OH */
    rf[70] = rf[70] * c__[5] * c__[13];
    rb[70] *= c__[6];
    /*     R71: OH + CH3OH = H2O + CH3O */
    rf[71] = rf[71] * c__[5] * c__[13];
    rb[71] *= c__[6];
    /*     R72: O2 + CH3OH = HO2 + CH2OH */
    rf[72] = rf[72] * c__[4] * c__[13];
    rb[72] *= c__[7];
    /*     R73: HO2 + CH3OH = H2O2 + CH2OH */
    rf[73] = rf[73] * c__[7] * c__[13];
    rb[73] *= c__[8];
    /*     R74: CH3OH + CH3 = CH2OH + CH4 */
    rf[74] = rf[74] * c__[13] * c__[15];
    rb[74] *= c__[14];
    /*     R75: CH3O = CH2OH */
    rf[75] *= c__[13];
    rb[75] *= c__[13];
    /*     R76: 2CH3O = CH2O + CH3OH */
    rb[76] = rb[76] * c__[11] * c__[13];
    /*     R77: H + CH3 = CH4 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[12] * ctb / rf[77];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 74.) * .217 + exp(-(*t) / 2941.) * .783 + exp(-6964. /
            *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[77] *= pcor;
    rb[77] *= pcor;
    rf[77] = rf[77] * c__[1] * c__[15];
    rb[77] *= c__[14];
    /*     R78: H + CH4 = H2 + CH3 */
    rf[78] = rf[78] * c__[1] * c__[14];
    rb[78] = rb[78] * c__[2] * c__[15];
    /*     R79: OH + CH4 = H2O + CH3 */
    rf[79] = rf[79] * c__[5] * c__[14];
    rb[79] = rb[79] * c__[6] * c__[15];
    /*     R80: O + CH4 = OH + CH3 */
    rf[80] = rf[80] * c__[3] * c__[14];
    rb[80] = rb[80] * c__[5] * c__[15];
    /*     R81: HO2 + CH4 = H2O2 + CH3 */
    rf[81] = rf[81] * c__[7] * c__[14];
    rb[81] = rb[81] * c__[8] * c__[15];
    /*     R82: CH4 + CH2 = 2CH3 */
    rf[82] = rf[82] * c__[14] * c__[16];
    rb[82] = rb[82] * c__[15] * c__[15];
    /*     R83: OH + CH3 = H2O + CH2(S) */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = 33.83274658851973 - alogt * .669 + 224.3339904338762 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = 33.83274658851973 - alogt * .669 + 224.3339904338762 / *t;
        r2l = 34.72691433702608 - alogt * .778 + 88.3648468375699 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = 34.72691433702608 - alogt * .778 + 88.3648468375699 / *t;
        r2l = 40.80825139477372 - alogt * 1.518 - 891.6999350579376 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 40.80825139477372 - alogt * 1.518 - 891.6999350579376 / *t;
        r2l = 54.52556992655878 - alogt * 3.155 - 3524.026323482357 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = 54.52556992655878 - alogt * 3.155 - 3524.026323482357 / *t;
        r2l = 45.88126934747938 - alogt * 1.962 - 4148.518208023498 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 45.88126934747938 - alogt * 1.962 - 4148.518208023498 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[83] *= pcor;
    rb[83] *= pcor;
    rf[83] = rf[83] * c__[5] * c__[15];
    rb[83] *= c__[6];
    /*     R84: OH + CH3 = H2 + CH2O */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 1.441 + 12.76625969883389 + 1632.434869823899 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = alogt * 1.441 + 12.76625969883389 + 1632.434869823899 / *t;
        r2l = alogt * 1.327 + 13.69379479928018 + 1497.069586228761 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = alogt * 1.327 + 13.69379479928018 + 1497.069586228761 / *t;
        r2l = alogt * .973 + 16.61887093887081 + 1011.465501956239 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * .973 + 16.61887093887081 + 1011.465501956239 / *t;
        r2l = alogt * .287 + 22.40483834713971 - 140.9006669391775 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .287 + 22.40483834713971 - 140.9006669391775 / *t;
        r2l = 43.69719169402195 - alogt * 2.199 - 4915.923626174375 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 43.69719169402195 - alogt * 2.199 - 4915.923626174375 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[84] *= pcor;
    rb[84] *= pcor;
    rf[84] = rf[84] * c__[5] * c__[15];
    rb[84] = rb[84] * c__[2] * c__[11];
    /*     R85: OH + CH3 = H + CH2OH */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * .965 + 23.508894172694 - 1617.338369794702 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = alogt * .965 + 23.508894172694 - 1617.338369794702 / *t;
        r2l = alogt * .95 + 23.61751894155055 - 1633.944519826819 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = alogt * .95 + 23.61751894155055 - 1633.944519826819 / *t;
        r2l = alogt * .833 + 24.57043027002606 - 1794.470636803954 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * .833 + 24.57043027002606 - 1794.470636803954 / *t;
        r2l = alogt * .134 + 30.35560061898197 - 2838.645222156787 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .134 + 30.35560061898197 - 2838.645222156787 / *t;
        r2l = 33.51434350441683 - alogt * .186 - 4328.166558370949 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 33.51434350441683 - alogt * .186 - 4328.166558370949 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[85] *= pcor;
    rb[85] *= pcor;
    rf[85] = rf[85] * c__[5] * c__[15];
    rb[85] *= c__[1];
    /*     R86: OH + CH3 = H + CH3O */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 1.016 + 20.89385213752195 - 6008.407011620641 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = alogt * 1.016 + 20.89385213752195 - 6008.407011620641 / *t;
        r2l = alogt * 1.016 + 20.89553705788687 - 6008.407011620641 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = alogt * 1.016 + 20.89553705788687 - 6008.407011620641 / *t;
        r2l = alogt * 1.011 + 20.93028000633074 - 6013.439178297041 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 1.011 + 20.93028000633074 - 6013.439178297041 / *t;
        r2l = alogt * .965 + 21.30994077299584 - 6068.793011737432 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .965 + 21.30994077299584 - 6068.793011737432 / *t;
        r2l = alogt * .5510000000000001 + 24.68255403484167 -
            6577.04184605375 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * .5510000000000001 + 24.68255403484167 -
            6577.04184605375 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[86] *= pcor;
    rb[86] *= pcor;
    rf[86] = rf[86] * c__[5] * c__[15];
    rb[86] *= c__[1];
    /*     R87: HO2 + CH3 = OH + CH3O */
    rf[87] = rf[87] * c__[7] * c__[15];
    rb[87] *= c__[5];
    /*     R88: HO2 + CH3 = O2 + CH4 */
    rf[88] = rf[88] * c__[7] * c__[15];
    rb[88] = rb[88] * c__[4] * c__[14];
    /*     R89: O + CH3 = H + CH2O */
    rf[89] = rf[89] * c__[3] * c__[15];
    rb[89] = rb[89] * c__[1] * c__[11];
    /*     R90: O2 + CH3 = O + CH3O */
    rf[90] = rf[90] * c__[4] * c__[15];
    rb[90] *= c__[3];
    /*     R91: O2 + CH3 = OH + CH2O */
    rf[91] = rf[91] * c__[4] * c__[15];
    rb[91] = rb[91] * c__[5] * c__[11];
    /*     R92: CH2(S) = CH2 */
    rf[92] *= c__[86];
    rb[92] = rb[92] * c__[16] * c__[86];
    /*     R93: O + CH2(S) = H2 + CO */
    rf[93] *= c__[3];
    rb[93] = rb[93] * c__[2] * c__[9];
    /*     R94: O + CH2(S) = H + HCO */
    rf[94] *= c__[3];
    rb[94] = rb[94] * c__[1] * c__[12];
    /*     R95: OH + CH2(S) = H + CH2O */
    rf[95] *= c__[5];
    rb[95] = rb[95] * c__[1] * c__[11];
    /*     R96: H2 + CH2(S) = H + CH3 */
    rf[96] *= c__[2];
    rb[96] = rb[96] * c__[1] * c__[15];
    /*     R97: O2 + CH2(S) = H + OH + CO */
    rf[97] *= c__[4];
    /*     R98: O2 + CH2(S) = H2O + CO */
    rf[98] *= c__[4];
    rb[98] = rb[98] * c__[6] * c__[9];
    /*     R99: CH2(S) = CH2 */
    rf[99] *= c__[6];
    rb[99] = rb[99] * c__[6] * c__[16];
    /*     R100: CH2(S) = CH2 */
    rf[100] *= c__[9];
    rb[100] = rb[100] * c__[9] * c__[16];
    /*     R101: CH2(S) = CH2 */
    rf[101] *= c__[10];
    rb[101] = rb[101] * c__[10] * c__[16];
    /*     R102: CO2 + CH2(S) = CO + CH2O */
    rf[102] *= c__[10];
    rb[102] = rb[102] * c__[9] * c__[11];
    /*     R103: H + CH2 = CH3 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[13] * ctb / rf[103];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 78.) * .32 + exp(-(*t) / 1995.) * .6800000000000001 +
        exp(-5590. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[103] *= pcor;
    rb[103] *= pcor;
    rf[103] = rf[103] * c__[1] * c__[16];
    rb[103] *= c__[15];
    /*     R104: O2 + CH2 = OH + HCO */
    rf[104] = rf[104] * c__[4] * c__[16];
    rb[104] = rb[104] * c__[5] * c__[12];
    /*     R105: O2 + CH2 = 2H + CO2 */
    rf[105] = rf[105] * c__[4] * c__[16];
    /*     R106: O + CH2 = 2H + CO */
    rf[106] = rf[106] * c__[3] * c__[16];
    /*     R107: 2CH3 = C2H6 */
    ctb = ctot + c__[6] * 4. + c__[9] + c__[10] * 2.;
    pr = rklow[14] * ctb / rf[107];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 570.) * 1. + exp(-1e30 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[107] *= pcor;
    rb[107] *= pcor;
    rf[107] = rf[107] * c__[15] * c__[15];
    rb[107] *= c__[17];
    /*     R108: H + C2H5 = C2H6 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[15] * ctb / rf[108];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 125.) * .158 + exp(-(*t) / 2219.) * .842 + exp(-6882.
            / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[108] *= pcor;
    rb[108] *= pcor;
    rf[108] = rf[108] * c__[1] * c__[18];
    rb[108] *= c__[17];
    /*     R109: H + C2H6 = H2 + C2H5 */
    rf[109] = rf[109] * c__[1] * c__[17];
    rb[109] = rb[109] * c__[2] * c__[18];
    /*     R110: O + C2H6 = OH + C2H5 */
    rf[110] = rf[110] * c__[3] * c__[17];
    rb[110] = rb[110] * c__[5] * c__[18];
    /*     R111: OH + C2H6 = H2O + C2H5 */
    rf[111] = rf[111] * c__[5] * c__[17];
    rb[111] = rb[111] * c__[6] * c__[18];
    /*     R112: O2 + C2H6 = HO2 + C2H5 */
    rf[112] = rf[112] * c__[4] * c__[17];
    rb[112] = rb[112] * c__[7] * c__[18];
    /*     R113: CH3 + C2H6 = CH4 + C2H5 */
    rf[113] = rf[113] * c__[15] * c__[17];
    rb[113] = rb[113] * c__[14] * c__[18];
    /*     R114: HO2 + C2H6 = H2O2 + C2H5 */
    rf[114] = rf[114] * c__[7] * c__[17];
    rb[114] = rb[114] * c__[8] * c__[18];
    /*     R115: CH3O + C2H6 = CH3OH + C2H5 */
    rf[115] *= c__[17];
    rb[115] = rb[115] * c__[13] * c__[18];
    /*     R116: CH2(S) + C2H6 = CH3 + C2H5 */
    rf[116] *= c__[17];
    rb[116] = rb[116] * c__[15] * c__[18];
    /*     R117: H + C2H4 = C2H5 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[16] * ctb / rf[117];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 299.) * 1.569 - exp(*t / 9147.) * .569 + exp(-152.4 /
            *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[117] *= pcor;
    rb[117] *= pcor;
    rf[117] = rf[117] * c__[1] * c__[19];
    rb[117] *= c__[18];
    /*     R118: 2C2H4 = C2H5 + C2H3 */
    rf[118] = rf[118] * c__[19] * c__[19];
    rb[118] = rb[118] * c__[18] * c__[20];
    /*     R119: CH3 + C2H5 = CH4 + C2H4 */
    rf[119] = rf[119] * c__[15] * c__[18];
    rb[119] = rb[119] * c__[14] * c__[19];
    /*     R120: 2CH3 = H + C2H5 */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * .105 + 29.18705825163553 - 5366.453508712395 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = alogt * .105 + 29.18705825163553 - 5366.453508712395 / *t;
        r2l = 30.87751210782972 - alogt * .096 - 5739.739632767688 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = 30.87751210782972 - alogt * .096 - 5739.739632767688 / *t;
        r2l = 33.36759341340774 - alogt * .362 - 6729.264888014826 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 33.36759341340774 - alogt * .362 - 6729.264888014826 / *t;
        r2l = alogt * .885 + 23.79131877208003 - 6809.779554837213 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .885 + 23.79131877208003 - 6809.779554837213 / *t;
        r2l = alogt * 3.23 + 4.636668853047462 - 5654.192799268902 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * 3.23 + 4.636668853047462 - 5654.192799268902 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[120] *= pcor;
    rb[120] *= pcor;
    rf[120] = rf[120] * c__[15] * c__[15];
    rb[120] = rb[120] * c__[1] * c__[18];
    /*     R121: H + C2H5 = H2 + C2H4 */
    rf[121] = rf[121] * c__[1] * c__[18];
    rb[121] = rb[121] * c__[2] * c__[19];
    /*     R122: O + C2H5 = H + CH3CHO */
    rf[122] = rf[122] * c__[3] * c__[18];
    rb[122] = rb[122] * c__[1] * c__[22];
    /*     R123: O2 + C2H5 = HO2 + C2H4 */
    /*     Reaction of PLOG type */
    if (p <= 41300.) {
        r1l = alogt * .49 + 21.46234194939476 + 196.9590037142646 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = alogt * .49 + 21.46234194939476 + 196.9590037142646 / *t;
        r2l = alogt * 1.13 + 16.72949032964601 + 362.6179307013262 / *t;
        dpl = 3.218875824868201;
        p1l = 10.62861777894912;
    } else if (p < 1.0325e7) {
        r1l = alogt * 1.13 + 16.72949032964601 + 362.6179307013262 / *t;
        r2l = 34.25919475849278 - alogt * 1.01 - 2389.775954621979 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 34.25919475849278 - alogt * 1.01 - 2389.775954621979 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[123] *= pcor;
    rb[123] *= pcor;
    rf[123] = rf[123] * c__[4] * c__[18];
    rb[123] = rb[123] * c__[7] * c__[19];
    /*     R124: O2 + C2H5 = HO2 + C2H4 */
    rf[124] = rf[124] * c__[4] * c__[18];
    rb[124] = rb[124] * c__[7] * c__[19];
    /*     R125: O2 + C2H5 = OH + C2H4O1-2 */
    /*     Reaction of PLOG type */
    if (p <= 41300.) {
        r1l = alogt * 1.93 + 7.172424577124845 + 252.9670188225876 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = alogt * 1.93 + 7.172424577124845 + 252.9670188225876 / *t;
        r2l = alogt * 2.18 + 5.496348217047172 + 31.45104172749498 / *t;
        dpl = 3.218875824868201;
        p1l = 10.62861777894912;
    } else if (p < 1.0325e7) {
        r1l = alogt * 2.18 + 5.496348217047172 + 31.45104172749498 / *t;
        r2l = alogt * .15 + 22.25387696883454 - 2721.898955264326 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = alogt * .15 + 22.25387696883454 - 2721.898955264326 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[125] *= pcor;
    rb[125] *= pcor;
    rf[125] = rf[125] * c__[4] * c__[18];
    rb[125] = rb[125] * c__[5] * c__[27];
    /*     R126: O2 + C2H5 = OH + CH3CHO */
    /*     Reaction of PLOG type */
    if (p <= 41300.) {
        r1l = alogt * 4.76 - 12.22464403111561 - 127.9679985808316 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = alogt * 4.76 - 12.22464403111561 - 127.9679985808316 / *t;
        r2l = alogt * 3.57 - 2.687806494625168 - 1330.001652572308 / *t;
        dpl = 3.218875824868201;
        p1l = 10.62861777894912;
    } else if (p < 1.0325e7) {
        r1l = alogt * 3.57 - 2.687806494625168 - 1330.001652572308 / *t;
        r2l = alogt * 2.41 + 6.717199917261079 - 2659.500088476976 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = alogt * 2.41 + 6.717199917261079 - 2659.500088476976 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[126] *= pcor;
    rb[126] *= pcor;
    rf[126] = rf[126] * c__[4] * c__[18];
    rb[126] = rb[126] * c__[5] * c__[22];
    /*     R127: C2H4O1-2 = HCO + CH3 */
    rf[127] *= c__[27];
    rb[127] = rb[127] * c__[12] * c__[15];
    /*     R128: C2H4O1-2 = CH3CHO */
    rf[128] *= c__[27];
    rb[128] *= c__[22];
    /*     R129: OH + C2H4O1-2 = H2O + C2H3O1-2 */
    rf[129] = rf[129] * c__[5] * c__[27];
    rb[129] *= c__[6];
    /*     R130: H + C2H4O1-2 = H2 + C2H3O1-2 */
    rf[130] = rf[130] * c__[1] * c__[27];
    rb[130] *= c__[2];
    /*     R131: HO2 + C2H4O1-2 = H2O2 + C2H3O1-2 */
    rf[131] = rf[131] * c__[7] * c__[27];
    rb[131] *= c__[8];
    /*     R132: CH3 + C2H4O1-2 = CH4 + C2H3O1-2 */
    rf[132] = rf[132] * c__[15] * c__[27];
    rb[132] *= c__[14];
    /*     R133: CH3O + C2H4O1-2 = CH3OH + C2H3O1-2 */
    rf[133] *= c__[27];
    rb[133] *= c__[13];
    /*     R134: C2H3O1-2 = CH3CO */
    /*     R135: C2H3O1-2 = CH2CHO */
    rb[135] *= c__[24];
    /*     R136: CH3CHO = HCO + CH3 */
    ctb = ctot;
    pr = rklow[17] * ctb / rf[136];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 718.1) * .99751 + exp(-(*t) / 6.089) * .00249 + exp(
            -3780. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[136] *= pcor;
    rb[136] *= pcor;
    rf[136] *= c__[22];
    rb[136] = rb[136] * c__[12] * c__[15];
    /*     R137: CH3CHO = CO + CH4 */
    ctb = ctot;
    pr = rklow[18] * ctb / rf[137];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 718.1) * .99751 + exp(-(*t) / 6.089) * .00249 + exp(
            -3780. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[137] *= pcor;
    rb[137] *= pcor;
    rf[137] *= c__[22];
    rb[137] = rb[137] * c__[9] * c__[14];
    /*     R138: H + CH3CHO = H2 + CH3CO */
    rf[138] = rf[138] * c__[1] * c__[22];
    rb[138] *= c__[2];
    /*     R139: H + CH3CHO = H2 + CH2CHO */
    rf[139] = rf[139] * c__[1] * c__[22];
    rb[139] = rb[139] * c__[2] * c__[24];
    /*     R140: O + CH3CHO = OH + CH3CO */
    rf[140] = rf[140] * c__[3] * c__[22];
    rb[140] *= c__[5];
    /*     R141: OH + CH3CHO = H2O + CH3CO */
    rf[141] = rf[141] * c__[5] * c__[22];
    rb[141] *= c__[6];
    /*     R142: O2 + CH3CHO = HO2 + CH3CO */
    rf[142] = rf[142] * c__[4] * c__[22];
    rb[142] *= c__[7];
    /*     R143: CH3 + CH3CHO = CH4 + CH3CO */
    rf[143] = rf[143] * c__[15] * c__[22];
    rb[143] *= c__[14];
    /*     R144: HO2 + CH3CHO = H2O2 + CH3CO */
    rf[144] = rf[144] * c__[7] * c__[22];
    rb[144] *= c__[8];
    /*     R145: OH + CH3CHO = H2O + CH2CHO */
    rf[145] = rf[145] * c__[5] * c__[22];
    rb[145] = rb[145] * c__[6] * c__[24];
    /*     R146: CH3CO = CO + CH3 */
    ctb = ctot;
    pr = rklow[19] * ctb / rf[146];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 8.73e9) * .371 + exp(-(*t) / 5.52) * .629 + exp(
            -7.6e7 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[146] *= pcor;
    rb[146] *= pcor;
    rb[146] = rb[146] * c__[9] * c__[15];
    /*     R147: H + CH3CO = H2 + CH2CO */
    rf[147] *= c__[1];
    rb[147] = rb[147] * c__[2] * c__[25];
    /*     R148: O + CH3CO = OH + CH2CO */
    rf[148] *= c__[3];
    rb[148] = rb[148] * c__[5] * c__[25];
    /*     R149: CH3 + CH3CO = CH4 + CH2CO */
    rf[149] *= c__[15];
    rb[149] = rb[149] * c__[14] * c__[25];
    /*     R150: CH2CHO = H + CH2CO */
    ctb = ctot;
    pr = rklow[20] * ctb / rf[150];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 393.) * .01500000000000001 + exp(-(*t) / 9.8e9) *
        .985 + exp(-5e9 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[150] *= pcor;
    rb[150] *= pcor;
    rf[150] *= c__[24];
    rb[150] = rb[150] * c__[1] * c__[25];
    /*     R151: CH2CHO = CO + CH3 */
    ctb = ctot;
    pr = rklow[21] * ctb / rf[151];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 1150.) * .9999999999999999 + exp(-(*t) / 4.99e9) *
        7.13e-17 + exp(-1.79e9 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[151] *= pcor;
    rb[151] *= pcor;
    rf[151] *= c__[24];
    rb[151] = rb[151] * c__[9] * c__[15];
    /*     R152: O2 + CH2CHO = HO2 + CH2CO */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 2.37 + 12.14419724181209 - 11941.33152309529 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = alogt * 2.37 + 12.14419724181209 - 11941.33152309529 / *t;
        r2l = alogt * 2.37 + 12.14419724181209 - 13773.0401933046 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = alogt * 2.37 + 12.14419724181209 - 13773.0401933046 / *t;
        r2l = alogt * 2.33 + 12.43320821811392 - 11976.55668983009 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 2.33 + 12.43320821811392 - 11976.55668983009 / *t;
        r2l = alogt * 1.63 + 18.0711232677825 - 12726.34952461357 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = alogt * 1.63 + 18.0711232677825 - 12726.34952461357 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[152] *= pcor;
    rb[152] *= pcor;
    rf[152] = rf[152] * c__[4] * c__[24];
    rb[152] = rb[152] * c__[7] * c__[25];
    /*     R153: O2 + CH2CHO = OH + CO + CH2O */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = 40.12976337542154 - alogt * 1.84 - 3286.004839688676 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 103250.) {
        r1l = 40.12976337542154 - alogt * 1.84 - 3286.004839688676 / *t;
        r2l = 46.4704121947391 - alogt * 2.58 - 4518.885675406479 / *t;
        dpl = 2.302585092994046;
        p1l = 9.242323417829233;
    } else if (p < 1032500.) {
        r1l = 46.4704121947391 - alogt * 2.58 - 4518.885675406479 / *t;
        r2l = 44.24989205479936 - alogt * 2.22 - 5203.26034339677 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 44.24989205479936 - alogt * 2.22 - 5203.26034339677 / *t;
        r2l = 32.12559488057461 - alogt * .6 - 5092.552676515988 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 32.12559488057461 - alogt * .6 - 5092.552676515988 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[153] *= pcor;
    rb[153] *= pcor;
    rf[153] = rf[153] * c__[4] * c__[24];
    /*     R154: CO + CH2 = CH2CO */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[22] * ctb / rf[154];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 275.) * .4093 + exp(-(*t) / 1226.) * .5907 + exp(
            -5185. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[154] *= pcor;
    rb[154] *= pcor;
    rf[154] = rf[154] * c__[9] * c__[16];
    rb[154] *= c__[25];
    /*     R155: CH3CO = H + CH2CO */
    ctb = ctot;
    pr = rklow[23] * ctb / rf[155];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 8.103e9) * .3991 + exp(-(*t) / 667.7000000000001) *
        .6009 + exp(-5e9 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[155] *= pcor;
    rb[155] *= pcor;
    rb[155] = rb[155] * c__[1] * c__[25];
    /*     R156: H + CH2CO = H2 + HCCO */
    rf[156] = rf[156] * c__[1] * c__[25];
    rb[156] = rb[156] * c__[2] * c__[26];
    /*     R157: H + CH2CO = CO + CH3 */
    rf[157] = rf[157] * c__[1] * c__[25];
    rb[157] = rb[157] * c__[9] * c__[15];
    /*     R158: O + CH2CO = CO2 + CH2 */
    rf[158] = rf[158] * c__[3] * c__[25];
    rb[158] = rb[158] * c__[10] * c__[16];
    /*     R159: O + CH2CO = OH + HCCO */
    rf[159] = rf[159] * c__[3] * c__[25];
    rb[159] = rb[159] * c__[5] * c__[26];
    /*     R160: OH + CH2CO = H2O + HCCO */
    rf[160] = rf[160] * c__[5] * c__[25];
    rb[160] = rb[160] * c__[6] * c__[26];
    /*     R161: OH + CH2CO = CO + CH2OH */
    rf[161] = rf[161] * c__[5] * c__[25];
    rb[161] *= c__[9];
    /*     R162: CH3 + CH2CO = CO + C2H5 */
    rf[162] = rf[162] * c__[15] * c__[25];
    rb[162] = rb[162] * c__[9] * c__[18];
    /*     R163: CH2(S) + CH2CO = CO + C2H4 */
    rf[163] *= c__[25];
    rb[163] = rb[163] * c__[9] * c__[19];
    /*     R164: OH + HCCO = H2 + 2CO */
    rf[164] = rf[164] * c__[5] * c__[26];
    /*     R165: O + HCCO = H + 2CO */
    rf[165] = rf[165] * c__[3] * c__[26];
    /*     R166: H + HCCO = CO + CH2(S) */
    rf[166] = rf[166] * c__[1] * c__[26];
    rb[166] *= c__[9];
    /*     R167: O2 + HCCO = OH + 2CO */
    rf[167] = rf[167] * c__[4] * c__[26];
    /*     R168: O2 + HCCO = H + CO + CO2 */
    rf[168] = rf[168] * c__[4] * c__[26];
    /*     R169: H + C2H3 = C2H4 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[24] * ctb / rf[169];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 207.5) * .218 + exp(-(*t) / 2663.) * .782 + exp(
            -6095. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[169] *= pcor;
    rb[169] *= pcor;
    rf[169] = rf[169] * c__[1] * c__[20];
    rb[169] *= c__[19];
    /*     R170: C2H4 = H2 + H2CC */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[25] * ctb / rf[170];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 180.) * .2655 + exp(-(*t) / 1035.) * .7345 + exp(
            -5417. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[170] *= pcor;
    rb[170] *= pcor;
    rf[170] *= c__[19];
    rb[170] *= c__[2];
    /*     R171: H + C2H4 = H2 + C2H3 */
    rf[171] = rf[171] * c__[1] * c__[19];
    rb[171] = rb[171] * c__[2] * c__[20];
    /*     R172: O + C2H4 = HCO + CH3 */
    rf[172] = rf[172] * c__[3] * c__[19];
    rb[172] = rb[172] * c__[12] * c__[15];
    /*     R173: O + C2H4 = H + CH2CHO */
    rf[173] = rf[173] * c__[3] * c__[19];
    rb[173] = rb[173] * c__[1] * c__[24];
    /*     R174: OH + C2H4 = H2O + C2H3 */
    rf[174] = rf[174] * c__[5] * c__[19];
    rb[174] = rb[174] * c__[6] * c__[20];
    /*     R175: OH + C2H4 = CH2O + CH3 */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 2.92 + 1.677096560907915 + 871.9235200196888 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 25812.5) {
        r1l = alogt * 2.92 + 1.677096560907915 + 871.9235200196888 / *t;
        r2l = alogt * 2.71 + 3.462606009790799 + 589.9208994742778 / *t;
        dpl = .9162907318741551;
        p1l = 9.242323417829233;
    } else if (p < 103250.) {
        r1l = alogt * 2.71 + 3.462606009790799 + 589.9208994742778 / *t;
        r2l = alogt * 2.36 + 6.318968113746434 + 90.98157350929749 / *t;
        dpl = 1.386294361119891;
        p1l = 10.15861414970339;
    } else if (p < 1032500.) {
        r1l = alogt * 2.36 + 6.318968113746434 + 90.98157350929749 / *t;
        r2l = alogt * 1.68 + 12.08953882927422 - 1036.877943672055 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 1.68 + 12.08953882927422 - 1036.877943672055 / *t;
        r2l = alogt * .5600000000000001 + 21.58615579209345 -
            3022.671557512706 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .5600000000000001 + 21.58615579209345 -
            3022.671557512706 / *t;
        r2l = 30.94883688865165 - alogt * .5 - 5764.397249482044 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 30.94883688865165 - alogt * .5 - 5764.397249482044 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[175] *= pcor;
    rb[175] *= pcor;
    rf[175] = rf[175] * c__[5] * c__[19];
    rb[175] = rb[175] * c__[11] * c__[15];
    /*     R176: OH + C2H4 = H + CH3CHO */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 5.3 - 15.25520569581128 + 1031.896098662419 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 25812.5) {
        r1l = alogt * 5.3 - 15.25520569581128 + 1031.896098662419 / *t;
        r2l = alogt * 4.57 - 9.346160095118718 + 310.9879006014704 / *t;
        dpl = .9162907318741551;
        p1l = 9.242323417829233;
    } else if (p < 103250.) {
        r1l = alogt * 4.57 - 9.346160095118718 + 310.9879006014704 / *t;
        r2l = alogt * 3.54 - .908818717035454 - 946.9028034980369 / *t;
        dpl = 1.386294361119891;
        p1l = 10.15861414970339;
    } else if (p < 1032500.) {
        r1l = alogt * 3.54 - .908818717035454 - 946.9028034980369 / *t;
        r2l = alogt * 3.91 - 3.738069698304708 - 866.8913533432897 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 3.91 - 3.738069698304708 - 866.8913533432897 / *t;
        r2l = alogt * 1.01 + 20.53089394429896 - 5287.448491892927 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * 1.01 + 20.53089394429896 - 5287.448491892927 / *t;
        r2l = alogt * .8100000000000001 + 22.64018844912847 -
            6978.256495163058 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * .8100000000000001 + 22.64018844912847 -
            6978.256495163058 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[176] *= pcor;
    rb[176] *= pcor;
    rf[176] = rf[176] * c__[5] * c__[19];
    rb[176] = rb[176] * c__[1] * c__[22];
    /*     R177: OH + C2H4 = H + C2H3OH */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 2.6 + 9.249561085129464 - 2073.755887344109 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 25812.5) {
        r1l = alogt * 2.6 + 9.249561085129464 - 2073.755887344109 / *t;
        r2l = alogt * 2.6 + 9.277999020449997 - 2077.781620685229 / *t;
        dpl = .9162907318741551;
        p1l = 9.242323417829233;
    } else if (p < 103250.) {
        r1l = alogt * 2.6 + 9.277999020449997 - 2077.781620685229 / *t;
        r2l = alogt * 2.56 + 9.629050706834368 - 2132.783202458272 / *t;
        dpl = 1.386294361119891;
        p1l = 10.15861414970339;
    } else if (p < 1032500.) {
        r1l = alogt * 2.56 + 9.629050706834368 - 2132.783202458272 / *t;
        r2l = alogt * 2.19 + 12.67294638176698 - 2644.705518448362 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 2.19 + 12.67294638176698 - 2644.705518448362 / *t;
        r2l = alogt * 1.43 + 19.0833687170276 - 3939.582647619403 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * 1.43 + 19.0833687170276 - 3939.582647619403 / *t;
        r2l = alogt * .75 + 25.17178221288913 - 5782.362084516789 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * .75 + 25.17178221288913 - 5782.362084516789 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[177] *= pcor;
    rb[177] *= pcor;
    rf[177] = rf[177] * c__[5] * c__[19];
    rb[177] = rb[177] * c__[1] * c__[23];
    /*     R178: O2 + C2H3OH = HO2 + CH2CHO */
    rf[178] = rf[178] * c__[4] * c__[23];
    rb[178] = rb[178] * c__[7] * c__[24];
    /*     R179: O + C2H3OH = OH + CH2CHO */
    rf[179] = rf[179] * c__[3] * c__[23];
    rb[179] = rb[179] * c__[5] * c__[24];
    /*     R180: OH + C2H3OH = H2O + CH2CHO */
    rf[180] = rf[180] * c__[5] * c__[23];
    rb[180] = rb[180] * c__[6] * c__[24];
    /*     R181: CH3 + C2H3OH = CH4 + CH2CHO */
    rf[181] = rf[181] * c__[15] * c__[23];
    rb[181] = rb[181] * c__[14] * c__[24];
    /*     R182: H + C2H3OH = H2 + CH2CHO */
    rf[182] = rf[182] * c__[1] * c__[23];
    rb[182] = rb[182] * c__[2] * c__[24];
    /*     R183: C2H3OH = CH3CHO */
    rf[183] = rf[183] * c__[7] * c__[23];
    rb[183] = rb[183] * c__[7] * c__[22];
    /*     R184: C2H3OH = CH3CHO */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 107.9230933349054 - alogt * 10.56 - 33926.86773228338 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 107.9230933349054 - alogt * 10.56 - 33926.86773228338 / *t;
        r2l = 98.19471360183952 - alogt * 9.09 - 33750.3393252753 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e8) {
        r1l = 98.19471360183952 - alogt * 9.09 - 33750.3393252753 / *t;
        r2l = 63.23450824783166 - alogt * 4.35 - 31004.63822163161 / *t;
        dpl = 4.605170185988092;
        p1l = 13.84749360381733;
    } else {
        r1l = 63.23450824783166 - alogt * 4.35 - 31004.63822163161 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[184] *= pcor;
    rb[184] *= pcor;
    rf[184] *= c__[23];
    rb[184] *= c__[22];
    /*     R185: CH3 + C2H4 = CH4 + C2H3 */
    rf[185] = rf[185] * c__[15] * c__[19];
    rb[185] = rb[185] * c__[14] * c__[20];
    /*     R186: O2 + C2H4 = HO2 + C2H3 */
    rf[186] = rf[186] * c__[4] * c__[19];
    rb[186] = rb[186] * c__[7] * c__[20];
    /*     R187: CH3O + C2H4 = CH3OH + C2H3 */
    rf[187] *= c__[19];
    rb[187] = rb[187] * c__[13] * c__[20];
    /*     R188: HO2 + C2H4 = OH + C2H4O1-2 */
    rf[188] = rf[188] * c__[7] * c__[19];
    rb[188] = rb[188] * c__[5] * c__[27];
    /*     R189: CH3 + CH2(S) = H + C2H4 */
    rf[189] *= c__[15];
    rb[189] = rb[189] * c__[1] * c__[19];
    /*     R190: H + C2H2 = C2H3 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[26] * ctb / rf[190];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(*t / 10200.) * .212 + exp(-(*t) / 1e-30) * .788;
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[190] *= pcor;
    rb[190] *= pcor;
    rf[190] = rf[190] * c__[1] * c__[21];
    rb[190] *= c__[20];
    /*     R191: O2 + C2H3 = CH2O + HCO */
    rf[191] = rf[191] * c__[4] * c__[20];
    rb[191] = rb[191] * c__[11] * c__[12];
    /*     R192: O2 + C2H3 = O + CH2CHO */
    rf[192] = rf[192] * c__[4] * c__[20];
    rb[192] = rb[192] * c__[3] * c__[24];
    /*     R193: O2 + C2H3 = H + CO + CH2O */
    rf[193] = rf[193] * c__[4] * c__[20];
    /*     R194: CH3 + C2H3 = CH4 + C2H2 */
    rf[194] = rf[194] * c__[15] * c__[20];
    rb[194] = rb[194] * c__[14] * c__[21];
    /*     R195: H + C2H3 = H2 + C2H2 */
    rf[195] = rf[195] * c__[1] * c__[20];
    rb[195] = rb[195] * c__[2] * c__[21];
    /*     R196: H + C2H3 = H2 + H2CC */
    rf[196] = rf[196] * c__[1] * c__[20];
    rb[196] *= c__[2];
    /*     R197: OH + C2H3 = H2O + C2H2 */
    rf[197] = rf[197] * c__[5] * c__[20];
    rb[197] = rb[197] * c__[6] * c__[21];
    /*     R198: 2C2H3 = C2H4 + C2H2 */
    rf[198] = rf[198] * c__[20] * c__[20];
    rb[198] = rb[198] * c__[19] * c__[21];
    /*     R199: H + C2H = C2H2 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[27] * ctb / rf[199];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 132.) * .354 + exp(-(*t) / 1315.) * .646 + exp(-5566.
            / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[199] *= pcor;
    rb[199] *= pcor;
    rf[199] *= c__[1];
    rb[199] *= c__[21];
    /*     R200: OH + C2H = H + HCCO */
    rf[200] *= c__[5];
    rb[200] = rb[200] * c__[1] * c__[26];
    /*     R201: O2 + C2H = CO + HCO */
    rf[201] *= c__[4];
    rb[201] = rb[201] * c__[9] * c__[12];
    /*     R202: H2 + C2H = H + C2H2 */
    rf[202] *= c__[2];
    rb[202] = rb[202] * c__[1] * c__[21];
    /*     R203: C2H2 = H2CC */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2. + c__[21] * 1.5 + c__[19] * 1.5;
    pr = rklow[28] * ctb / rf[203];
    pcor = pr / (pr + 1.f);
    rf[203] *= pcor;
    rb[203] *= pcor;
    rf[203] *= c__[21];
    /*     R204: O + C2H2 = CO + CH2 */
    rf[204] = rf[204] * c__[3] * c__[21];
    rb[204] = rb[204] * c__[9] * c__[16];
    /*     R205: O + C2H2 = H + HCCO */
    rf[205] = rf[205] * c__[3] * c__[21];
    rb[205] = rb[205] * c__[1] * c__[26];
    /*     R206: OH + C2H2 = H2O + C2H */
    rf[206] = rf[206] * c__[5] * c__[21];
    rb[206] *= c__[6];
    /*     R207: OH + C2H2 = H + CH2CO */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 2.56 + 7.363913501405819 + 424.9664758219122 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 25812.5) {
        r1l = alogt * 2.56 + 7.363913501405819 + 424.9664758219122 / *t;
        r2l = alogt * 2.28 + 9.627734050949622 + 146.9895886176205 / *t;
        dpl = .9162907318741551;
        p1l = 9.242323417829233;
    } else if (p < 103250.) {
        r1l = alogt * 2.28 + 9.627734050949622 + 146.9895886176205 / *t;
        r2l = alogt * 1.92 + 12.61718842514715 - 300.973888915436 / *t;
        dpl = 1.386294361119891;
        p1l = 10.15861414970339;
    } else if (p < 1032500.) {
        r1l = alogt * 1.92 + 12.61718842514715 - 300.973888915436 / *t;
        r2l = alogt * 1.55 + 15.83413996024735 - 1059.774302049671 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * 1.55 + 15.83413996024735 - 1059.774302049671 / *t;
        r2l = alogt * 1.65 + 15.44494715690506 - 1710.936669975727 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * 1.65 + 15.44494715690506 - 1710.936669975727 / *t;
        r2l = alogt * 2.45 + 9.58671989918925 - 2252.90102102392 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * 2.45 + 9.58671989918925 - 2252.90102102392 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[207] *= pcor;
    rb[207] *= pcor;
    rf[207] = rf[207] * c__[5] * c__[21];
    rb[207] = rb[207] * c__[1] * c__[25];
    /*     R208: OH + C2H2 = CO + CH3 */
    /*     Reaction of PLOG type */
    if (p <= 10325.) {
        r1l = alogt * 1.68 + 13.07254268242037 + 165.9608569876455 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 25812.5) {
        r1l = alogt * 1.68 + 13.07254268242037 + 165.9608569876455 / *t;
        r2l = alogt * 1.4 + 15.29073112827857 - 113.9785752204418 / *t;
        dpl = .9162907318741551;
        p1l = 9.242323417829233;
    } else if (p < 103250.) {
        r1l = alogt * 1.4 + 15.29073112827857 - 113.9785752204418 / *t;
        r2l = alogt * 1.05 + 18.15253982670742 - 561.0865844185105 / *t;
        dpl = 1.386294361119891;
        p1l = 10.15861414970339;
    } else if (p < 1032500.) {
        r1l = alogt * 1.05 + 18.15253982670742 - 561.0865844185105 / *t;
        r2l = alogt * .73 + 20.96777941399681 - 1297.795785843353 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = alogt * .73 + 20.96777941399681 - 1297.795785843353 / *t;
        r2l = alogt * .92 + 19.88208257755906 - 1880.01747030274 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = alogt * .92 + 19.88208257755906 - 1880.01747030274 / *t;
        r2l = alogt * 1.77 + 13.62313866531682 - 2363.608687904703 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * 1.77 + 13.62313866531682 - 2363.608687904703 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[208] *= pcor;
    rb[208] *= pcor;
    rf[208] = rf[208] * c__[5] * c__[21];
    rb[208] = rb[208] * c__[9] * c__[15];
    /*     R209: HCO + C2H2 = CO + C2H3 */
    rf[209] = rf[209] * c__[12] * c__[21];
    rb[209] = rb[209] * c__[9] * c__[20];
    /*     R210: CH2 + C2H2 = H + C3H3 */
    rf[210] = rf[210] * c__[16] * c__[21];
    rb[210] = rb[210] * c__[1] * c__[35];
    /*     R211: CH2(S) + C2H2 = H + C3H3 */
    rf[211] *= c__[21];
    rb[211] = rb[211] * c__[1] * c__[35];
    /*     R212: C2H2 + HCCO = CO + C3H3 */
    rf[212] = rf[212] * c__[21] * c__[26];
    rb[212] = rb[212] * c__[9] * c__[35];
    /*     R213: H2CC = C2H2 */
    rf[213] *= c__[1];
    rb[213] = rb[213] * c__[1] * c__[21];
    /*     R214: OH + H2CC = H + CH2CO */
    rf[214] *= c__[5];
    rb[214] = rb[214] * c__[1] * c__[25];
    /*     R215: O2 + H2CC = 2HCO */
    rf[215] *= c__[4];
    rb[215] = rb[215] * c__[12] * c__[12];
    /*     R216: HCO + C2H3 = C2H3CHO */
    rf[216] = rf[216] * c__[12] * c__[20];
    rb[216] *= c__[28];
    /*     R217: C3H8 = CH3 + C2H5 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[9] * .5 + c__[10] + c__[14] + c__[
        17] * 2.;
    pr = rklow[29] * ctb / rf[217];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 50.) * .69 + exp(-(*t) / 3e3) * .31 + exp(-9e3 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[217] *= pcor;
    rb[217] *= pcor;
    rf[217] *= c__[29];
    rb[217] = rb[217] * c__[15] * c__[18];
    /*     R218: H + NC3H7 = C3H8 */
    rf[218] *= c__[1];
    rb[218] *= c__[29];
    /*     R219: O2 + C3H8 = HO2 + NC3H7 */
    rf[219] = rf[219] * c__[4] * c__[29];
    rb[219] *= c__[7];
    /*     R220: H + C3H8 = H2 + NC3H7 */
    rf[220] = rf[220] * c__[1] * c__[29];
    rb[220] *= c__[2];
    /*     R221: O + C3H8 = OH + NC3H7 */
    rf[221] = rf[221] * c__[3] * c__[29];
    rb[221] *= c__[5];
    /*     R222: OH + C3H8 = H2O + NC3H7 */
    rf[222] = rf[222] * c__[5] * c__[29];
    rb[222] *= c__[6];
    /*     R223: HO2 + C3H8 = H2O2 + NC3H7 */
    rf[223] = rf[223] * c__[7] * c__[29];
    rb[223] *= c__[8];
    /*     R224: CH3 + C3H8 = CH4 + NC3H7 */
    rf[224] = rf[224] * c__[15] * c__[29];
    rb[224] *= c__[14];
    /*     R225: C2H3 + C3H8 = C2H4 + NC3H7 */
    rf[225] = rf[225] * c__[20] * c__[29];
    rb[225] *= c__[19];
    /*     R226: C2H5 + C3H8 = C2H6 + NC3H7 */
    rf[226] = rf[226] * c__[18] * c__[29];
    rb[226] *= c__[17];
    /*     R227: C3H8 + C3H5-A = NC3H7 + C3H6 */
    rf[227] = rf[227] * c__[29] * c__[31];
    rb[227] *= c__[30];
    /*     R228: CH3O + C3H8 = CH3OH + NC3H7 */
    rf[228] *= c__[29];
    rb[228] *= c__[13];
    /*     R229: CH3 + C2H4 = NC3H7 */
    rf[229] = rf[229] * c__[15] * c__[19];
    /*     R230: H + C3H6 = NC3H7 */
    rf[230] = rf[230] * c__[1] * c__[30];
    /*     R231: O2 + NC3H7 = HO2 + C3H6 */
    rf[231] *= c__[4];
    rb[231] = rb[231] * c__[7] * c__[30];
    /*     R232: CH3 + C2H3 = C3H6 */
    ctb = ctot;
    pr = rklow[30] * ctb / rf[232];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 1341.) * .825 + exp(-(*t) / 6e4) * .175 + exp(-10140.
            / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[232] *= pcor;
    rb[232] *= pcor;
    rf[232] = rf[232] * c__[15] * c__[20];
    rb[232] *= c__[30];
    /*     R233: H + C3H5-A = C3H6 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[31] * ctb / rf[233];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 1097.) * .98 + exp(-(*t) / 1097.) * .02 + exp(-6860. /
            *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[233] *= pcor;
    rb[233] *= pcor;
    rf[233] = rf[233] * c__[1] * c__[31];
    rb[233] *= c__[30];
    /*     R234: C3H6 = H + C3H5-S */
    rf[234] *= c__[30];
    rb[234] = rb[234] * c__[1] * c__[32];
    /*     R235: C3H6 = H + C3H5-T */
    rf[235] *= c__[30];
    rb[235] = rb[235] * c__[1] * c__[33];
    /*     R236: O + C3H6 = HCO + C2H5 */
    rf[236] = rf[236] * c__[3] * c__[30];
    rb[236] = rb[236] * c__[12] * c__[18];
    /*     R237: O + C3H6 = H + CH3 + CH2CO */
    rf[237] = rf[237] * c__[3] * c__[30];
    /*     R238: O + C3H6 = OH + C3H5-A */
    rf[238] = rf[238] * c__[3] * c__[30];
    rb[238] = rb[238] * c__[5] * c__[31];
    /*     R239: O + C3H6 = OH + C3H5-S */
    rf[239] = rf[239] * c__[3] * c__[30];
    rb[239] = rb[239] * c__[5] * c__[32];
    /*     R240: O + C3H6 = OH + C3H5-T */
    rf[240] = rf[240] * c__[3] * c__[30];
    rb[240] = rb[240] * c__[5] * c__[33];
    /*     R241: OH + C3H6 = H2O + C3H5-A */
    rf[241] = rf[241] * c__[5] * c__[30];
    rb[241] = rb[241] * c__[6] * c__[31];
    /*     R242: OH + C3H6 = H2O + C3H5-S */
    rf[242] = rf[242] * c__[5] * c__[30];
    rb[242] = rb[242] * c__[6] * c__[32];
    /*     R243: OH + C3H6 = H2O + C3H5-T */
    rf[243] = rf[243] * c__[5] * c__[30];
    rb[243] = rb[243] * c__[6] * c__[33];
    /*     R244: HO2 + C3H6 = H2O2 + C3H5-A */
    rf[244] = rf[244] * c__[7] * c__[30];
    rb[244] = rb[244] * c__[8] * c__[31];
    /*     R245: HO2 + C3H6 = H2O2 + C3H5-S */
    rf[245] = rf[245] * c__[7] * c__[30];
    rb[245] = rb[245] * c__[8] * c__[32];
    /*     R246: HO2 + C3H6 = H2O2 + C3H5-T */
    rf[246] = rf[246] * c__[7] * c__[30];
    rb[246] = rb[246] * c__[8] * c__[33];
    /*     R247: H + C3H6 = H2 + C3H5-A */
    rf[247] = rf[247] * c__[1] * c__[30];
    rb[247] = rb[247] * c__[2] * c__[31];
    /*     R248: H + C3H6 = H2 + C3H5-S */
    rf[248] = rf[248] * c__[1] * c__[30];
    rb[248] = rb[248] * c__[2] * c__[32];
    /*     R249: H + C3H6 = H2 + C3H5-T */
    rf[249] = rf[249] * c__[1] * c__[30];
    rb[249] = rb[249] * c__[2] * c__[33];
    /*     R250: H + C3H6 = CH3 + C2H4 */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 39.01611320938889 - alogt * 1.05 - 3251.282889621521 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 39.01611320938889 - alogt * 1.05 - 3251.282889621521 / *t;
        r2l = 50.43372849455479 - alogt * 2.39 - 5625.962344214302 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 50.43372849455479 - alogt * 2.39 - 5625.962344214302 / *t;
        r2l = 56.45596470032953 - alogt * 3.04 - 7855.212181859146 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 56.45596470032953 - alogt * 3.04 - 7855.212181859146 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[250] *= pcor;
    rb[250] *= pcor;
    rf[250] = rf[250] * c__[1] * c__[30];
    rb[250] = rb[250] * c__[15] * c__[19];
    /*     R251: O2 + C3H6 = HO2 + C3H5-A */
    rf[251] = rf[251] * c__[4] * c__[30];
    rb[251] = rb[251] * c__[7] * c__[31];
    /*     R252: O2 + C3H6 = HO2 + C3H5-S */
    rf[252] = rf[252] * c__[4] * c__[30];
    rb[252] = rb[252] * c__[7] * c__[32];
    /*     R253: O2 + C3H6 = HO2 + C3H5-T */
    rf[253] = rf[253] * c__[4] * c__[30];
    rb[253] = rb[253] * c__[7] * c__[33];
    /*     R254: CH3 + C3H6 = CH4 + C3H5-A */
    rf[254] = rf[254] * c__[15] * c__[30];
    rb[254] = rb[254] * c__[14] * c__[31];
    /*     R255: CH3 + C3H6 = CH4 + C3H5-S */
    rf[255] = rf[255] * c__[15] * c__[30];
    rb[255] = rb[255] * c__[14] * c__[32];
    /*     R256: CH3 + C3H6 = CH4 + C3H5-T */
    rf[256] = rf[256] * c__[15] * c__[30];
    rb[256] = rb[256] * c__[14] * c__[33];
    /*     R257: C2H5 + C3H6 = C2H6 + C3H5-A */
    rf[257] = rf[257] * c__[18] * c__[30];
    rb[257] = rb[257] * c__[17] * c__[31];
    /*     R258: CH3 + C2H2 = C3H5-A */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 124.1411440829546 - alogt * 13.32 - 16706.79336564534 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 124.1411440829546 - alogt * 13.32 - 16706.79336564534 / *t;
        r2l = 123.0228267232072 - alogt * 12.82 - 17979.93153477433 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = 123.0228267232072 - alogt * 12.82 - 17979.93153477433 / *t;
        r2l = 121.026408517339 - alogt * 12.46 - 18179.70855182738 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = 121.026408517339 - alogt * 12.46 - 18179.70855182738 / *t;
        r2l = 117.4710604558496 - alogt * 11.89 - 18355.33116883371 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = 117.4710604558496 - alogt * 11.89 - 18355.33116883371 / *t;
        r2l = 114.3082740976325 - alogt * 11.4 - 18468.05170238505 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = 114.3082740976325 - alogt * 11.4 - 18468.05170238505 / *t;
        r2l = 102.6487451584704 - alogt * 9.630000000000001 -
            18920.94670326098 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 102.6487451584704 - alogt * 9.630000000000001 -
            18920.94670326098 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[258] *= pcor;
    rb[258] *= pcor;
    rf[258] = rf[258] * c__[15] * c__[21];
    rb[258] *= c__[31];
    /*     R259: O + C3H5-A = H + C2H3CHO */
    rf[259] = rf[259] * c__[3] * c__[31];
    rb[259] = rb[259] * c__[1] * c__[28];
    /*     R260: OH + C3H5-A = 2H + C2H3CHO */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 86.86335526133777 - alogt * 6.71 - 14747.26766185549 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 86.86335526133777 - alogt * 6.71 - 14747.26766185549 / *t;
        r2l = 75.11780750109878 - alogt * 5.16 - 15159.90532932022 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 75.11780750109878 - alogt * 5.16 - 15159.90532932022 / *t;
        r2l = 46.52170548912665 - alogt * 1.56 - 13249.69485895909 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 46.52170548912665 - alogt * 1.56 - 13249.69485895909 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[260] *= pcor;
    rb[260] *= pcor;
    rf[260] = rf[260] * c__[5] * c__[31];
    /*     R261: C2H5 + C3H5-A = C2H4 + C3H6 */
    rf[261] = rf[261] * c__[18] * c__[31];
    rb[261] = rb[261] * c__[19] * c__[30];
    /*     R262: O2 + C3H5-A = CH2O + CH3CO */
    /*     Reaction of PLOG type */
    if (p <= 1032500.) {
        r1l = 34.71272970203413 - alogt * 1.01 - 10128.7450862563 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1.0325e7) {
        r1l = 34.71272970203413 - alogt * 1.01 - 10128.7450862563 / *t;
        r2l = 36.50448917126218 - alogt * 1.21 - 10590.69798714975 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 36.50448917126218 - alogt * 1.21 - 10590.69798714975 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[262] *= pcor;
    rb[262] *= pcor;
    rf[262] = rf[262] * c__[4] * c__[31];
    rb[262] *= c__[11];
    /*     R263: O2 + C3H5-A = OH + C2H3CHO */
    /*     Reaction of PLOG type */
    if (p <= 1032500.) {
        r1l = 30.5324427100113 - alogt * .41 - 11503.02980558092 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1.0325e7) {
        r1l = 30.5324427100113 - alogt * .41 - 11503.02980558092 / *t;
        r2l = 30.83782435956248 - alogt * .45 - 11582.53803906803 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 30.83782435956248 - alogt * .45 - 11582.53803906803 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[263] *= pcor;
    rb[263] *= pcor;
    rf[263] = rf[263] * c__[4] * c__[31];
    rb[263] = rb[263] * c__[5] * c__[28];
    /*     R264: HCO + C3H5-A = CO + C3H6 */
    rf[264] = rf[264] * c__[12] * c__[31];
    rb[264] = rb[264] * c__[9] * c__[30];
    /*     R265: CH3 + C2H3 = H + C3H5-A */
    rf[265] = rf[265] * c__[15] * c__[20];
    rb[265] = rb[265] * c__[1] * c__[31];
    /*     R266: C3H5-A = C3H5-T */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 137.2134970397843 - alogt * 15.42 - 37942.53674004995 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 137.2134970397843 - alogt * 15.42 - 37942.53674004995 / *t;
        r2l = 130.8992102591717 - alogt * 14.08 - 38178.04214050543 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = 130.8992102591717 - alogt * 14.08 - 38178.04214050543 / *t;
        r2l = 128.2107960325864 - alogt * 13.59 - 38218.80269058426 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = 128.2107960325864 - alogt * 13.59 - 38218.80269058426 / *t;
        r2l = 123.6180483665968 - alogt * 12.81 - 38185.59039052003 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = 123.6180483665968 - alogt * 12.81 - 38185.59039052003 / *t;
        r2l = 119.288137733062 - alogt * 12.12 - 38093.50174034192 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = 119.288137733062 - alogt * 12.12 - 38093.50174034192 / *t;
        r2l = 100.0407784159251 - alogt * 9.27 - 37238.03340535406 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 100.0407784159251 - alogt * 9.27 - 37238.03340535406 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[266] *= pcor;
    rb[266] *= pcor;
    rf[266] *= c__[31];
    rb[266] *= c__[33];
    /*     R267: C3H5-A = C3H5-S */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 126.90454437914 - alogt * 14.53 - 37137.39007182607 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 126.90454437914 - alogt * 14.53 - 37137.39007182607 / *t;
        r2l = 119.0412776551304 - alogt * 13.02 - 36885.78173800611 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 119.0412776551304 - alogt * 13.02 - 36885.78173800611 / *t;
        r2l = 112.7962103492235 - alogt * 11.73 - 37087.06840506208 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = 112.7962103492235 - alogt * 11.73 - 37087.06840506208 / *t;
        r2l = 102.8947825296504 - alogt * 9.84 - 36936.1034047701 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 102.8947825296504 - alogt * 9.84 - 36936.1034047701 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[267] *= pcor;
    rb[267] *= pcor;
    rf[267] *= c__[31];
    rb[267] *= c__[32];
    /*     R268: CH3 + C2H2 = C3H5-S */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 74.01919521243067 - alogt * 7.14 - 5032.166676399197 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 74.01919521243067 - alogt * 7.14 - 5032.166676399197 / *t;
        r2l = 81.75362906459728 - alogt * 7.76 - 6692.781679610932 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 81.75362906459728 - alogt * 7.76 - 6692.781679610932 / *t;
        r2l = 88.37370227112764 - alogt * 8.210000000000001 -
            8605.005016642626 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = 88.37370227112764 - alogt * 8.210000000000001 -
            8605.005016642626 / *t;
        r2l = 90.137290863389 - alogt * 8.060000000000001 - 10164.97668632638
            / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 90.137290863389 - alogt * 8.060000000000001 - 10164.97668632638
            / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[268] *= pcor;
    rb[268] *= pcor;
    rf[268] = rf[268] * c__[15] * c__[21];
    rb[268] *= c__[32];
    /*     R269: O2 + C3H5-S = HCO + CH3CHO */
    rf[269] = rf[269] * c__[4] * c__[32];
    rb[269] = rb[269] * c__[12] * c__[22];
    /*     R270: CH3 + C2H2 = C3H5-T */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 47.96862447206298 - alogt * 4.16 - 9057.900017518554 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 47.96862447206298 - alogt * 4.16 - 9057.900017518554 / *t;
        r2l = 52.26430795563243 - alogt * 4.39 - 9485.634185012486 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = 52.26430795563243 - alogt * 4.39 - 9485.634185012486 / *t;
        r2l = 54.75121660809111 - alogt * 4.6 - 9848.453402380868 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = 54.75121660809111 - alogt * 4.6 - 9848.453402380868 / *t;
        r2l = 59.55387059861283 - alogt * 5.06 - 10643.0325205843 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = 59.55387059861283 - alogt * 5.06 - 10643.0325205843 / *t;
        r2l = 64.39981191099844 - alogt * 5.55 - 11523.66168895416 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = 64.39981191099844 - alogt * 5.55 - 11523.66168895416 / *t;
        r2l = 84.22806441451799 - alogt * 7.58 - 15750.68169712949 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 84.22806441451799 - alogt * 7.58 - 15750.68169712949 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[270] *= pcor;
    rb[270] *= pcor;
    rf[270] = rf[270] * c__[15] * c__[21];
    rb[270] *= c__[33];
    /*     R271: H + C3H4-P = C3H5-T */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 102.8398003952331 - alogt * 10.21 - 5132.810009927181 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 102.8398003952331 - alogt * 10.21 - 5132.810009927181 / *t;
        r2l = 108.7283169730886 - alogt * 10.58 - 6889.0361799905 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = 108.7283169730886 - alogt * 10.58 - 6889.0361799905 / *t;
        r2l = 109.8389054528034 - alogt * 10.61 - 7400.807530980299 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = 109.8389054528034 - alogt * 10.61 - 7400.807530980299 / *t;
        r2l = 110.4853436353978 - alogt * 10.55 - 8006.177182151122 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = 110.4853436353978 - alogt * 10.55 - 8006.177182151122 / *t;
        r2l = 110.1674095197755 - alogt * 10.4 - 8353.396682822668 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = 110.1674095197755 - alogt * 10.4 - 8353.396682822668 / *t;
        r2l = 102.4768949015437 - alogt * 9.109999999999999 -
            8755.970016934602 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 102.4768949015437 - alogt * 9.109999999999999 -
            8755.970016934602 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[271] *= pcor;
    rb[271] *= pcor;
    rf[271] = rf[271] * c__[1] * c__[34];
    rb[271] *= c__[33];
    /*     R272: H + C3H5-S = H2 + C3H4-P */
    rf[272] = rf[272] * c__[1] * c__[32];
    rb[272] = rb[272] * c__[2] * c__[34];
    /*     R273: O + C3H5-S = HCO + C2H4 */
    rf[273] = rf[273] * c__[3] * c__[32];
    rb[273] = rb[273] * c__[12] * c__[19];
    /*     R274: OH + C3H5-S = H + HCO + C2H4 */
    rf[274] = rf[274] * c__[5] * c__[32];
    /*     R275: HO2 + C3H5-S = OH + HCO + C2H4 */
    rf[275] = rf[275] * c__[7] * c__[32];
    /*     R276: HCO + C3H5-S = CO + C3H6 */
    rf[276] = rf[276] * c__[12] * c__[32];
    rb[276] = rb[276] * c__[9] * c__[30];
    /*     R277: CH3 + C3H5-S = CH4 + C3H4-P */
    rf[277] = rf[277] * c__[15] * c__[32];
    rb[277] = rb[277] * c__[14] * c__[34];
    /*     R278: C3H5-T = C3H5-S */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 101.7837477209837 - alogt * 12.16 - 26267.91005080381 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 101.7837477209837 - alogt * 12.16 - 26267.91005080381 / *t;
        r2l = 110.9295495718224 - alogt * 12.71 - 27123.37838579167 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 110.9295495718224 - alogt * 12.71 - 27123.37838579167 / *t;
        r2l = 121.3636653754207 - alogt * 13.37 - 28783.9933890034 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = 121.3636653754207 - alogt * 13.37 - 28783.9933890034 / *t;
        r2l = 119.1896976602487 - alogt * 12.43 - 29790.42672428325 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 119.1896976602487 - alogt * 12.43 - 29790.42672428325 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[278] *= pcor;
    rb[278] *= pcor;
    rf[278] *= c__[33];
    rb[278] *= c__[32];
    /*     R279: O2 + C3H5-T = CH2O + CH3CO */
    rf[279] = rf[279] * c__[4] * c__[33];
    rb[279] *= c__[11];
    /*     R280: H + C3H5-T = H2 + C3H4-P */
    rf[280] = rf[280] * c__[1] * c__[33];
    rb[280] = rb[280] * c__[2] * c__[34];
    /*     R281: CH3 + C3H5-T = CH4 + C3H4-P */
    rf[281] = rf[281] * c__[15] * c__[33];
    rb[281] = rb[281] * c__[14] * c__[34];
    /*     R282: O + C3H5-T = CH3 + CH2CO */
    rf[282] = rf[282] * c__[3] * c__[33];
    rb[282] = rb[282] * c__[15] * c__[25];
    /*     R283: OH + C3H5-T = H + CH3 + CH2CO */
    rf[283] = rf[283] * c__[5] * c__[33];
    /*     R284: HO2 + C3H5-T = OH + CH3 + CH2CO */
    rf[284] = rf[284] * c__[7] * c__[33];
    /*     R285: HCO + C3H5-T = CO + C3H6 */
    rf[285] = rf[285] * c__[12] * c__[33];
    rb[285] = rb[285] * c__[9] * c__[30];
    /*     R286: H + C3H4-P = C3H5-S */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 57.56462732485115 - alogt * 5. - 905.7900017518555 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 57.56462732485115 - alogt * 5. - 905.7900017518555 / *t;
        r2l = 66.1771306960717 - alogt * 5.74 - 2163.831670851655 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 1.0325e7) {
        r1l = 66.1771306960717 - alogt * 5.74 - 2163.831670851655 / *t;
        r2l = 78.28789316179756 - alogt * 6.88 - 4478.628341995285 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else if (p < 1.0325e8) {
        r1l = 78.28789316179756 - alogt * 6.88 - 4478.628341995285 / *t;
        r2l = 87.46777432628903 - alogt * 7.63 - 6944.390013430892 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 87.46777432628903 - alogt * 7.63 - 6944.390013430892 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[286] *= pcor;
    rb[286] *= pcor;
    rf[286] = rf[286] * c__[1] * c__[34];
    rb[286] *= c__[32];
    /*     R287: H + C3H4-P = C3H5-A */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = 138.2504157594471 - alogt * 14.56 - 14140.38836068174 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = 138.2504157594471 - alogt * 14.56 - 14140.38836068174 / *t;
        r2l = 139.7463795214492 - alogt * 14.37 - 15923.78823079762 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = 139.7463795214492 - alogt * 14.37 - 15923.78823079762 / *t;
        r2l = 139.2669630950609 - alogt * 14.19 - 16425.99846510226 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = 139.2669630950609 - alogt * 14.19 - 16425.99846510226 / *t;
        r2l = 138.0519648207232 - alogt * 13.89 - 17085.71551637819 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = 138.0519648207232 - alogt * 13.89 - 17085.71551637819 / *t;
        r2l = 136.640977847013 - alogt * 13.61 - 17562.2617006332 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = 136.640977847013 - alogt * 13.61 - 17562.2617006332 / *t;
        r2l = 127.1121837439183 - alogt * 12.07 - 18870.62503649699 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = 127.1121837439183 - alogt * 12.07 - 18870.62503649699 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[287] *= pcor;
    rb[287] *= pcor;
    rf[287] = rf[287] * c__[1] * c__[34];
    rb[287] *= c__[31];
    /*     R288: H + C3H4-P = H2 + C3H3 */
    rf[288] = rf[288] * c__[1] * c__[34];
    rb[288] = rb[288] * c__[2] * c__[35];
    /*     R289: O + C3H4-P = CH3 + HCCO */
    rf[289] = rf[289] * c__[3] * c__[34];
    rb[289] = rb[289] * c__[15] * c__[26];
    /*     R290: O + C3H4-P = CO + C2H4 */
    rf[290] = rf[290] * c__[3] * c__[34];
    rb[290] = rb[290] * c__[9] * c__[19];
    /*     R291: OH + C3H4-P = H2O + C3H3 */
    rf[291] = rf[291] * c__[5] * c__[34];
    rb[291] = rb[291] * c__[6] * c__[35];
    /*     R292: C2H + C3H4-P = C2H2 + C3H3 */
    rf[292] *= c__[34];
    rb[292] = rb[292] * c__[21] * c__[35];
    /*     R293: CH3 + C3H4-P = CH4 + C3H3 */
    rf[293] = rf[293] * c__[15] * c__[34];
    rb[293] = rb[293] * c__[14] * c__[35];
    /*     R294: CH3 + C2H2 = H + C3H4-P */
    /*     Reaction of PLOG type */
    if (p <= 103250.) {
        r1l = alogt * 1.86 + 15.31958795474055 - 5837.313344623069 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1032500.) {
        r1l = alogt * 1.86 + 15.31958795474055 - 5837.313344623069 / *t;
        r2l = alogt * 1.1 + 21.66327309543788 - 6865.888213279064 / *t;
        dpl = 2.302585092994046;
        p1l = 11.54490851082328;
    } else if (p < 2.065e6) {
        r1l = alogt * 1.1 + 21.66327309543788 - 6865.888213279064 / *t;
        r2l = alogt * .85 + 23.75339953721774 - 7253.868264029442 / *t;
        dpl = .6931471805599453;
        p1l = 13.84749360381733;
    } else if (p < 5162500.) {
        r1l = alogt * .85 + 23.75339953721774 - 7253.868264029442 / *t;
        r2l = alogt * .5600000000000001 + 26.2487187760782 -
            7776.207165039679 / *t;
        dpl = .9162907318741551;
        p1l = 14.54064078437727;
    } else if (p < 1.0325e7) {
        r1l = alogt * .5600000000000001 + 26.2487187760782 -
            7776.207165039679 / *t;
        r2l = alogt * .39 + 27.72633129573287 - 8152.110015766699 / *t;
        dpl = .6931471805599453;
        p1l = 15.45693151625143;
    } else if (p < 1.0325e8) {
        r1l = alogt * .39 + 27.72633129573287 - 8152.110015766699 / *t;
        r2l = alogt * .37 + 28.37295846065793 - 9108.221684282546 / *t;
        dpl = 2.302585092994046;
        p1l = 16.15007869681137;
    } else {
        r1l = alogt * .37 + 28.37295846065793 - 9108.221684282546 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[294] *= pcor;
    rb[294] *= pcor;
    rf[294] = rf[294] * c__[15] * c__[21];
    rb[294] = rb[294] * c__[1] * c__[34];
    /*     R295: H + C3H3 = C3H4-P */
    rf[295] = rf[295] * c__[1] * c__[35];
    rb[295] *= c__[34];
    /*     R296: CH3 + C2H = C3H4-P */
    rf[296] *= c__[15];
    rb[296] *= c__[34];
    /*     R297: HO2 + C3H3 = O2 + C3H4-P */
    rf[297] = rf[297] * c__[7] * c__[35];
    rb[297] = rb[297] * c__[4] * c__[34];
    /*     R298: HO2 + C3H4-P = OH + CO + C2H4 */
    rf[298] = rf[298] * c__[7] * c__[34];
    /*     R299: OH + C3H4-P = CH3 + CH2CO */
    rf[299] = rf[299] * c__[5] * c__[34];
    rb[299] = rb[299] * c__[15] * c__[25];
    /*     R300: O + C3H4-P = HCO + C2H3 */
    rf[300] = rf[300] * c__[3] * c__[34];
    rb[300] = rb[300] * c__[12] * c__[20];
    /*     R301: O + C3H4-P = OH + C3H3 */
    rf[301] = rf[301] * c__[3] * c__[34];
    rb[301] = rb[301] * c__[5] * c__[35];
    /*     R302: C2H3 + C3H4-P = C2H4 + C3H3 */
    rf[302] = rf[302] * c__[20] * c__[34];
    rb[302] = rb[302] * c__[19] * c__[35];
    /*     R303: C3H5-A + C3H4-P = C3H6 + C3H3 */
    rf[303] = rf[303] * c__[31] * c__[34];
    rb[303] = rb[303] * c__[30] * c__[35];
    /*     R304: O + C3H3 = CH2O + C2H */
    rf[304] = rf[304] * c__[3] * c__[35];
    rb[304] *= c__[11];
    /*     R305: O2 + C3H3 = HCO + CH2CO */
    rf[305] = rf[305] * c__[4] * c__[35];
    rb[305] = rb[305] * c__[12] * c__[25];
    /*     R306: HO2 + C3H3 = OH + CO + C2H3 */
    rf[306] = rf[306] * c__[7] * c__[35];
    /*     R307: HCO + C3H3 = CO + C3H4-P */
    rf[307] = rf[307] * c__[12] * c__[35];
    rb[307] = rb[307] * c__[9] * c__[34];
    /*     R308: C2H5 + C2H = CH3 + C3H3 */
    rf[308] *= c__[18];
    rb[308] = rb[308] * c__[15] * c__[35];
    /*     R309: C4H6 = H + C4H5-I */
    rf[309] *= c__[36];
    rb[309] = rb[309] * c__[1] * c__[41];
    /*     R310: C4H6 = H2 + C4H4 */
    rf[310] *= c__[36];
    rb[310] = rb[310] * c__[2] * c__[37];
    /*     R311: H + C4H6 = H2 + C4H5-I */
    rf[311] = rf[311] * c__[1] * c__[36];
    rb[311] = rb[311] * c__[2] * c__[41];
    /*     R312: H + C4H6 = C2H4 + C2H3 */
    /*     Reaction of PLOG type */
    if (p <= 1032500.) {
        r1l = 69.45598922554161 - alogt * 4.34 - 10893.13120440134 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 1.0325e7) {
        r1l = 69.45598922554161 - alogt * 4.34 - 10893.13120440134 / *t;
        r2l = 70.77316839849652 - alogt * 4.51 - 11008.87103795852 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 70.77316839849652 - alogt * 4.51 - 11008.87103795852 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[312] *= pcor;
    rb[312] *= pcor;
    rf[312] = rf[312] * c__[1] * c__[36];
    rb[312] = rb[312] * c__[19] * c__[20];
    /*     R313: H + C4H6 = CH3 + C3H4-P */
    rf[313] = rf[313] * c__[1] * c__[36];
    rb[313] = rb[313] * c__[15] * c__[34];
    /*     R314: O + C4H6 = OH + C4H5-I */
    rf[314] = rf[314] * c__[3] * c__[36];
    rb[314] = rb[314] * c__[5] * c__[41];
    /*     R315: O + C4H6 = C2H2 + C2H4O1-2 */
    rf[315] = rf[315] * c__[3] * c__[36];
    rb[315] = rb[315] * c__[21] * c__[27];
    /*     R316: OH + C4H6 = H2O + C4H5-I */
    rf[316] = rf[316] * c__[5] * c__[36];
    rb[316] = rb[316] * c__[6] * c__[41];
    /*     R317: CH3 + C4H6 = CH4 + C4H5-I */
    rf[317] = rf[317] * c__[15] * c__[36];
    rb[317] = rb[317] * c__[14] * c__[41];
    /*     R318: C2H3 + C4H6 = C2H4 + C4H5-I */
    rf[318] = rf[318] * c__[20] * c__[36];
    rb[318] = rb[318] * c__[19] * c__[41];
    /*     R319: C3H5-A + C4H6 = C3H6 + C4H5-I */
    rf[319] = rf[319] * c__[31] * c__[36];
    rb[319] = rb[319] * c__[30] * c__[41];
    /*     R320: C2H3 + C2H2 = H + C4H4 */
    /*     Reaction of PLOG type */
    if (p <= 13629.) {
        r1l = 31.9076872349446 - alogt * .48 - 3069.62167260351 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 27154.75) {
        r1l = 31.9076872349446 - alogt * .48 - 3069.62167260351 / *t;
        r2l = 33.84562921435074 - alogt * .71 - 3371.551673187462 / *t;
        dpl = .6893521095913937;
        p1l = 9.519955154427514;
    } else if (p < 123900.) {
        r1l = 33.84562921435074 - alogt * .71 - 3371.551673187462 / *t;
        r2l = 38.36741779139978 - alogt * 1.25 - 4227.020008175326 / *t;
        dpl = 1.517922803598327;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 38.36741779139978 - alogt * 1.25 - 4227.020008175326 / *t;
        r2l = 42.13967885445277 - alogt * 1.68 - 5334.096676983148 / *t;
        dpl = 2.120263536200091;
        p1l = 11.72723006761723;
    } else if (p < 1.0325e7) {
        r1l = 42.13967885445277 - alogt * 1.68 - 5334.096676983148 / *t;
        r2l = 38.43059669302132 - alogt * 1.13 - 5937.956678151052 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 38.43059669302132 - alogt * 1.13 - 5937.956678151052 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[320] *= pcor;
    rb[320] *= pcor;
    rf[320] = rf[320] * c__[20] * c__[21];
    rb[320] = rb[320] * c__[1] * c__[37];
    /*     R321: C2H3 + C2H2 = C4H5-I */
    /*     Reaction of PLOG type */
    if (p <= 13629.) {
        r1l = 79.89733107423166 - alogt * 8.42 - 3975.411674355365 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 27154.75) {
        r1l = 79.89733107423166 - alogt * 8.42 - 3975.411674355365 / *t;
        r2l = 83.63500069251502 - alogt * 8.779999999999999 -
            4579.271675523269 / *t;
        dpl = .6893521095913937;
        p1l = 9.519955154427514;
    } else if (p < 123900.) {
        r1l = 83.63500069251502 - alogt * 8.779999999999999 -
            4579.271675523269 / *t;
        r2l = 85.19564844077969 - alogt * 8.77 - 4931.523342871213 / *t;
        dpl = 1.517922803598327;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 85.19564844077969 - alogt * 8.77 - 4931.523342871213 / *t;
        r2l = 106.3889179069718 - alogt * 10.98 - 9359.830018102506 / *t;
        dpl = 2.120263536200091;
        p1l = 11.72723006761723;
    } else if (p < 1.0325e7) {
        r1l = 106.3889179069718 - alogt * 10.98 - 9359.830018102506 / *t;
        r2l = 123.6662504684147 - alogt * 12.64 - 14492.64002802969 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 123.6662504684147 - alogt * 12.64 - 14492.64002802969 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[321] *= pcor;
    rb[321] *= pcor;
    rf[321] = rf[321] * c__[20] * c__[21];
    rb[321] *= c__[41];
    /*     R322: 2C2H3 = C4H6 */
    /*     Reaction of PLOG type */
    if (p <= 27154.75) {
        r1l = 133.1932604497159 - alogt * 13.82 - 8871.206633824144 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 123900.) {
        r1l = 133.1932604497159 - alogt * 13.82 - 8871.206633824144 / *t;
        r2l = 120.1398899437985 - alogt * 11.97 - 8079.64681562655 / *t;
        dpl = 1.517922803598327;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 120.1398899437985 - alogt * 11.97 - 8079.64681562655 / *t;
        r2l = 97.11403901385809 - alogt * 8.84 - 6281.653662149118 / *t;
        dpl = 2.120263536200091;
        p1l = 11.72723006761723;
    } else {
        r1l = 97.11403901385809 - alogt * 8.84 - 6281.653662149118 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[322] *= pcor;
    rb[322] *= pcor;
    rf[322] = rf[322] * c__[20] * c__[20];
    rb[322] *= c__[36];
    /*     R323: 2C2H3 = H + C4H5-I */
    /*     Reaction of PLOG type */
    if (p <= 27154.75) {
        r1l = 69.48301789792953 - alogt * 4.95 - 6520.681579278079 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 123900.) {
        r1l = 69.48301789792953 - alogt * 4.95 - 6520.681579278079 / *t;
        r2l = 66.44646362985529 - alogt * 4.49 - 7182.411497224573 / *t;
        dpl = 1.517922803598327;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 66.44646362985529 - alogt * 4.49 - 7182.411497224573 / *t;
        r2l = 50.83919360266296 - alogt * 2.44 - 6870.920379955463 / *t;
        dpl = 2.120263536200091;
        p1l = 11.72723006761723;
    } else {
        r1l = 50.83919360266296 - alogt * 2.44 - 6870.920379955463 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[323] *= pcor;
    rb[323] *= pcor;
    rf[323] = rf[323] * c__[20] * c__[20];
    rb[323] = rb[323] * c__[1] * c__[41];
    /*     R324: H + C4H5-I = H2 + C4H4 */
    rf[324] = rf[324] * c__[1] * c__[41];
    rb[324] = rb[324] * c__[2] * c__[37];
    /*     R325: H + C4H5-I = CH3 + C3H3 */
    rf[325] = rf[325] * c__[1] * c__[41];
    rb[325] = rb[325] * c__[15] * c__[35];
    /*     R326: OH + C4H5-I = H2O + C4H4 */
    rf[326] = rf[326] * c__[5] * c__[41];
    rb[326] = rb[326] * c__[6] * c__[37];
    /*     R327: HCO + C4H5-I = CO + C4H6 */
    rf[327] = rf[327] * c__[12] * c__[41];
    rb[327] = rb[327] * c__[9] * c__[36];
    /*     R328: HO2 + C4H5-I = O2 + C4H6 */
    rf[328] = rf[328] * c__[7] * c__[41];
    rb[328] = rb[328] * c__[4] * c__[36];
    /*     R329: HO2 + C4H5-I = OH + C2H3 + CH2CO */
    rf[329] = rf[329] * c__[7] * c__[41];
    /*     R330: H2O2 + C4H5-I = HO2 + C4H6 */
    rf[330] = rf[330] * c__[8] * c__[41];
    rb[330] = rb[330] * c__[7] * c__[36];
    /*     R331: O2 + C4H5-I = CH2CHO + CH2CO */
    rf[331] = rf[331] * c__[4] * c__[41];
    rb[331] = rb[331] * c__[24] * c__[25];
    /*     R332: C4H5-2 = C4H5-I */
    rf[332] *= c__[42];
    rb[332] *= c__[41];
    /*     R333: C4H5-2 = C4H5-I */
    rf[333] = rf[333] * c__[1] * c__[42];
    rb[333] = rb[333] * c__[1] * c__[41];
    /*     R334: HO2 + C4H5-2 = OH + C2H2 + CH3CO */
    rf[334] = rf[334] * c__[7] * c__[42];
    /*     R335: O2 + C4H5-2 = CH3CO + CH2CO */
    rf[335] = rf[335] * c__[4] * c__[42];
    rb[335] *= c__[25];
    /*     R336: C4H6-2 = C4H6 */
    rf[336] *= c__[43];
    rb[336] *= c__[36];
    /*     R337: H + C4H6-2 = H2 + C4H5-2 */
    rf[337] = rf[337] * c__[1] * c__[43];
    rb[337] = rb[337] * c__[2] * c__[42];
    /*     R338: H + C4H6-2 = CH3 + C3H4-P */
    rf[338] = rf[338] * c__[1] * c__[43];
    rb[338] = rb[338] * c__[15] * c__[34];
    /*     R339: C4H6-2 = H + C4H5-2 */
    rf[339] *= c__[43];
    rb[339] = rb[339] * c__[1] * c__[42];
    /*     R340: CH3 + C4H6-2 = CH4 + C4H5-2 */
    rf[340] = rf[340] * c__[15] * c__[43];
    rb[340] = rb[340] * c__[14] * c__[42];
    /*     R341: H + C4H4 = C4H5-I */
    /*     Reaction of PLOG type */
    if (p <= 13629.) {
        r1l = 123.8452986998637 - alogt * 13.19 - 7145.676680486859 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 27154.75) {
        r1l = 123.8452986998637 - alogt * 13.19 - 7145.676680486859 / *t;
        r2l = 121.9961879341642 - alogt * 12.85 - 7195.998347250851 / *t;
        dpl = .6893521095913937;
        p1l = 9.519955154427514;
    } else if (p < 123900.) {
        r1l = 121.9961879341642 - alogt * 12.85 - 7195.998347250851 / *t;
        r2l = 120.4763621804198 - alogt * 12.44 - 7799.858348418755 / *t;
        dpl = 1.517922803598327;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 120.4763621804198 - alogt * 12.44 - 7799.858348418755 / *t;
        r2l = 119.0210749478129 - alogt * 11.92 - 8906.935017226579 / *t;
        dpl = 2.120263536200091;
        p1l = 11.72723006761723;
    } else if (p < 1.0325e7) {
        r1l = 119.0210749478129 - alogt * 11.92 - 8906.935017226579 / *t;
        r2l = 110.9295495718224 - alogt * 10.58 - 9460.47335163049 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 110.9295495718224 - alogt * 10.58 - 9460.47335163049 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[341] *= pcor;
    rb[341] *= pcor;
    rf[341] = rf[341] * c__[1] * c__[37];
    rb[341] *= c__[41];
    /*     R342: H + C4H4 = H2 + C4H3-N */
    rf[342] = rf[342] * c__[1] * c__[37];
    rb[342] = rb[342] * c__[2] * c__[39];
    /*     R343: H + C4H4 = H2 + C4H3-I */
    rf[343] = rf[343] * c__[1] * c__[37];
    rb[343] = rb[343] * c__[2] * c__[38];
    /*     R344: OH + C4H4 = H2O + C4H3-N */
    rf[344] = rf[344] * c__[5] * c__[37];
    rb[344] = rb[344] * c__[6] * c__[39];
    /*     R345: OH + C4H4 = H2O + C4H3-I */
    rf[345] = rf[345] * c__[5] * c__[37];
    rb[345] = rb[345] * c__[6] * c__[38];
    /*     R346: O + C4H4 = HCO + C3H3 */
    rf[346] = rf[346] * c__[3] * c__[37];
    rb[346] = rb[346] * c__[12] * c__[35];
    /*     R347: HCCO + C3H3 = CO + C4H4 */
    rf[347] = rf[347] * c__[26] * c__[35];
    rb[347] = rb[347] * c__[9] * c__[37];
    /*     R348: CH2 + C3H3 = H + C4H4 */
    rf[348] = rf[348] * c__[16] * c__[35];
    rb[348] = rb[348] * c__[1] * c__[37];
    /*     R349: C2H2 + C2H = C4H3-N */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2. + c__[21] * 1.5 + c__[19] * 1.5;
    pr = rklow[32] * ctb / rf[349];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 100.) * 0. + exp(-(*t) / 5613.) * 1. + exp(-13390. / *
            t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[349] *= pcor;
    rb[349] *= pcor;
    rf[349] *= c__[21];
    rb[349] *= c__[39];
    /*     R350: C4H3-N = C4H3-I */
    rf[350] *= c__[39];
    rb[350] *= c__[38];
    /*     R351: C4H3-N = C4H3-I */
    rf[351] = rf[351] * c__[1] * c__[39];
    rb[351] = rb[351] * c__[1] * c__[38];
    /*     R352: H + C4H3-N = C2H2 + H2CC */
    rf[352] = rf[352] * c__[1] * c__[39];
    rb[352] *= c__[21];
    /*     R353: H + C4H3-N = C4H4 */
    rf[353] = rf[353] * c__[1] * c__[39];
    rb[353] *= c__[37];
    /*     R354: H + C4H3-N = H2 + C4H2 */
    rf[354] = rf[354] * c__[1] * c__[39];
    rb[354] = rb[354] * c__[2] * c__[40];
    /*     R355: OH + C4H3-N = H2O + C4H2 */
    rf[355] = rf[355] * c__[5] * c__[39];
    rb[355] = rb[355] * c__[6] * c__[40];
    /*     R356: C2H2 + C2H = C4H3-I */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2. + c__[21] * 1.5 + c__[19] * 1.5;
    pr = rklow[33] * ctb / rf[356];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 100.) * 0. + exp(-(*t) / 5613.) * 1. + exp(-13390. / *
            t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[356] *= pcor;
    rb[356] *= pcor;
    rf[356] *= c__[21];
    rb[356] *= c__[38];
    /*     R357: H + C4H3-I = C2H2 + H2CC */
    rf[357] = rf[357] * c__[1] * c__[38];
    rb[357] *= c__[21];
    /*     R358: H + C4H3-I = C4H4 */
    rf[358] = rf[358] * c__[1] * c__[38];
    rb[358] *= c__[37];
    /*     R359: H + C4H3-I = H2 + C4H2 */
    rf[359] = rf[359] * c__[1] * c__[38];
    rb[359] = rb[359] * c__[2] * c__[40];
    /*     R360: OH + C4H3-I = H2O + C4H2 */
    rf[360] = rf[360] * c__[5] * c__[38];
    rb[360] = rb[360] * c__[6] * c__[40];
    /*     R361: O2 + C4H3-I = CH2CO + HCCO */
    rf[361] = rf[361] * c__[4] * c__[38];
    rb[361] = rb[361] * c__[25] * c__[26];
    /*     R362: C2H2 + C2H = H + C4H2 */
    rf[362] *= c__[21];
    rb[362] = rb[362] * c__[1] * c__[40];
    /*     R363: H + C4H2 = C4H3-N */
    rf[363] = rf[363] * c__[1] * c__[40];
    rb[363] *= c__[39];
    /*     R364: H + C4H2 = C4H3-I */
    rf[364] = rf[364] * c__[1] * c__[40];
    rb[364] *= c__[38];
    /*     R365: OH + C4H2 = H + H2C4O */
    rf[365] = rf[365] * c__[5] * c__[40];
    rb[365] = rb[365] * c__[1] * c__[44];
    /*     R366: H + H2C4O = C2H2 + HCCO */
    rf[366] = rf[366] * c__[1] * c__[44];
    rb[366] = rb[366] * c__[21] * c__[26];
    /*     R367: OH + H2C4O = CH2CO + HCCO */
    rf[367] = rf[367] * c__[5] * c__[44];
    rb[367] = rb[367] * c__[25] * c__[26];
    /*     R368: C2H2 + H2CC = C4H4 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2. + c__[21] * 2. + c__[19] * 2.;
    pr = rklow[34] * ctb / rf[368];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 56.) * .02000000000000002 + exp(-(*t) / 580.) * .98 +
        exp(-4164. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[368] *= pcor;
    rb[368] *= pcor;
    rf[368] *= c__[21];
    rb[368] *= c__[37];
    /*     R369: C2H4 + H2CC = C4H6 */
    rf[369] *= c__[19];
    rb[369] *= c__[36];
    /*     R370: C2H3 + C3H5-A = 2H + c-C5H6 */
    rf[370] = rf[370] * c__[20] * c__[31];
    /*     R371: C3H5-A + C3H3 = 2H + A1 */
    rf[371] = rf[371] * c__[31] * c__[35];
    /*     R372: 2C3H3 = H + A1- */
    rf[372] = rf[372] * c__[35] * c__[35];
    rb[372] = rb[372] * c__[1] * c__[47];
    /*     R373: 2C3H3 = A1 */
    rf[373] = rf[373] * c__[35] * c__[35];
    rb[373] *= c__[45];
    /*     R374: C2H3 + C4H6 = H + H2 + A1 */
    rf[374] = rf[374] * c__[20] * c__[36];
    /*     R375: C2H2 + C4H5-2 = H + A1 */
    rf[375] = rf[375] * c__[21] * c__[42];
    rb[375] = rb[375] * c__[1] * c__[45];
    /*     R376: C2H4 + C4H5-2 = CH3 + c-C5H6 */
    rf[376] = rf[376] * c__[19] * c__[42];
    rb[376] = rb[376] * c__[15] * c__[46];
    /*     R377: C2H2 + C4H3-N = A1- */
    rf[377] = rf[377] * c__[21] * c__[39];
    rb[377] *= c__[47];
    /*     R378: CH3 + C4H3-I = c-C5H6 */
    rf[378] = rf[378] * c__[15] * c__[38];
    rb[378] *= c__[46];
    /*     R379: H + A1 = H2 + A1- */
    rf[379] = rf[379] * c__[1] * c__[45];
    rb[379] = rb[379] * c__[2] * c__[47];
    /*     R380: HO2 + A1 = H2O2 + A1- */
    rf[380] = rf[380] * c__[7] * c__[45];
    rb[380] = rb[380] * c__[8] * c__[47];
    /*     R381: OH + A1 = H2O + A1- */
    rf[381] = rf[381] * c__[5] * c__[45];
    rb[381] = rb[381] * c__[6] * c__[47];
    /*     R382: A1- = H + C2H2 + C4H2 */
    rf[382] *= c__[47];
    /*     R383: CH2O + A1- = HCO + A1 */
    rf[383] = rf[383] * c__[11] * c__[47];
    rb[383] = rb[383] * c__[12] * c__[45];
    /*     R384: HCO + A1- = CO + A1 */
    rf[384] = rf[384] * c__[12] * c__[47];
    rb[384] = rb[384] * c__[9] * c__[45];
    /*     R385: HO2 + A1- = OH + C6H5O */
    rf[385] = rf[385] * c__[7] * c__[47];
    rb[385] = rb[385] * c__[5] * c__[48];
    /*     R386: O2 + A1- = O + C6H5O */
    rf[386] = rf[386] * c__[4] * c__[47];
    rb[386] = rb[386] * c__[3] * c__[48];
    /*     R387: O + C6H5O = CO2 + c-C5H5 */
    rf[387] = rf[387] * c__[3] * c__[48];
    rb[387] = rb[387] * c__[10] * c__[49];
    /*     R388: O2 + c-C5H6 = HO2 + c-C5H5 */
    rf[388] = rf[388] * c__[4] * c__[46];
    rb[388] = rb[388] * c__[7] * c__[49];
    /*     R389: HO2 + c-C5H6 = H2O2 + c-C5H5 */
    rf[389] = rf[389] * c__[7] * c__[46];
    rb[389] = rb[389] * c__[8] * c__[49];
    /*     R390: H + c-C5H6 = C2H2 + C3H5-A */
    rf[390] = rf[390] * c__[1] * c__[46];
    rb[390] = rb[390] * c__[21] * c__[31];
    /*     R391: CH3 + c-C5H6 = CH4 + c-C5H5 */
    rf[391] = rf[391] * c__[15] * c__[46];
    rb[391] = rb[391] * c__[14] * c__[49];
    /*     R392: O + c-C5H5 = H + C5H4O */
    rf[392] = rf[392] * c__[3] * c__[49];
    rb[392] = rb[392] * c__[1] * c__[50];
    /*     R393: OH + c-C5H5 = CO + C4H6 */
    rf[393] = rf[393] * c__[5] * c__[49];
    rb[393] = rb[393] * c__[9] * c__[36];
    /*     R394: O + C5H4O = CO2 + C4H4 */
    rf[394] = rf[394] * c__[3] * c__[50];
    rb[394] = rb[394] * c__[10] * c__[37];
    /*     R395: HCO + c-C5H6 = CH2O + c-C5H5 */
    rf[395] = rf[395] * c__[12] * c__[46];
    rb[395] = rb[395] * c__[11] * c__[49];
    /*     R396: C2H3 + c-C5H6 = CH3 + A1 */
    rf[396] = rf[396] * c__[20] * c__[46];
    rb[396] = rb[396] * c__[15] * c__[45];
    /*     R397: C4H5-I + c-C5H6 = C4H6 + c-C5H5 */
    rf[397] = rf[397] * c__[41] * c__[46];
    rb[397] = rb[397] * c__[36] * c__[49];
    /*     R398: O2 + c-C5H5 = O + C5H5O(2,4) */
    rf[398] = rf[398] * c__[4] * c__[49];
    rb[398] *= c__[3];
    /*     R399: O + c-C5H5 = C5H5O(2,4) */
    rf[399] = rf[399] * c__[3] * c__[49];
    /*     R400: OH + c-C5H5 = H + C5H5O(2,4) */
    rf[400] = rf[400] * c__[5] * c__[49];
    rb[400] *= c__[1];
    /*     R401: C3H3 + C4H6 = H + C6H5CH3 */
    rf[401] = rf[401] * c__[35] * c__[36];
    rb[401] = rb[401] * c__[1] * c__[52];
    /*     R402: H + C6H5CH2 = C6H5CH3 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[35] * ctb / rf[402];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 383.) * .569 + exp(-(*t) / 152.) * .431 + exp(-4730. /
            *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[402] *= pcor;
    rb[402] *= pcor;
    rf[402] = rf[402] * c__[1] * c__[51];
    rb[402] *= c__[52];
    /*     R403: O + C6H5CH2 = H + A1CHO */
    rf[403] = rf[403] * c__[3] * c__[51];
    rb[403] = rb[403] * c__[1] * c__[54];
    /*     R404: HO2 + C6H5CH2 = H + OH + A1CHO */
    rf[404] = rf[404] * c__[7] * c__[51];
    rb[404] = rb[404] * c__[1] * c__[5] * c__[54];
    /*     R405: H + A1CHO = HCO + A1 */
    rf[405] = rf[405] * c__[1] * c__[54];
    rb[405] = rb[405] * c__[12] * c__[45];
    /*     R406: O + A1 = H + C6H5O */
    rf[406] = rf[406] * c__[3] * c__[45];
    rb[406] = rb[406] * c__[1] * c__[48];
    /*     R407: H + C6H5O = CO + c-C5H6 */
    rf[407] = rf[407] * c__[1] * c__[48];
    rb[407] = rb[407] * c__[9] * c__[46];
    /*     R408: C6H5O = CO + c-C5H5 */
    rf[408] *= c__[48];
    rb[408] = rb[408] * c__[9] * c__[49];
    /*     R409: OH + c-C5H6 = H2O + c-C5H5 */
    rf[409] = rf[409] * c__[5] * c__[46];
    rb[409] = rb[409] * c__[6] * c__[49];
    /*     R410: O + c-C5H6 = OH + c-C5H5 */
    rf[410] = rf[410] * c__[3] * c__[46];
    rb[410] = rb[410] * c__[5] * c__[49];
    /*     R411: C2H2 + C3H3 = c-C5H5 */
    rf[411] = rf[411] * c__[21] * c__[35];
    rb[411] *= c__[49];
    /*     R412: O2 + c-C5H5 = OH + C5H4O */
    rf[412] = rf[412] * c__[4] * c__[49];
    rb[412] = rb[412] * c__[5] * c__[50];
    /*     R413: HO2 + c-C5H5 = OH + C5H5O(2,4) */
    rf[413] = rf[413] * c__[7] * c__[49];
    rb[413] *= c__[5];
    /*     R414: C5H5O(2,4) = H + C5H4O */
    rb[414] = rb[414] * c__[1] * c__[50];
    /*     R415: O2 + C5H5O(2,4) = HO2 + C5H4O */
    rf[415] *= c__[4];
    rb[415] = rb[415] * c__[7] * c__[50];
    /*     R416: C5H4O = CO + 2C2H2 */
    rf[416] *= c__[50];
    rb[416] = rb[416] * c__[9] * c__[21] * c__[21];
    /*     R417: C3H3 + C4H4 = C6H5CH2 */
    rf[417] = rf[417] * c__[35] * c__[37];
    rb[417] *= c__[51];
    /*     R418: H + C6H5CH3 = CH3 + A1 */
    rf[418] = rf[418] * c__[1] * c__[52];
    rb[418] = rb[418] * c__[15] * c__[45];
    /*     R419: HO2 + A1- = O2 + A1 */
    rf[419] = rf[419] * c__[7] * c__[47];
    rb[419] = rb[419] * c__[4] * c__[45];
    /*     R420: H + A1- = A1 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10];
    pr = rklow[36] * ctb / rf[420];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / .1) * 0. + exp(-(*t) / 584.9) * 1. + exp(-6113. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[420] *= pcor;
    rb[420] *= pcor;
    rf[420] = rf[420] * c__[1] * c__[47];
    rb[420] *= c__[45];
    /*     R421: HCO + A1- = A1CHO */
    rf[421] = rf[421] * c__[12] * c__[47];
    rb[421] *= c__[54];
    /*     R422: CH3 + A1- = C6H5CH3 */
    rf[422] = rf[422] * c__[15] * c__[47];
    rb[422] *= c__[52];
    /*     R423: H + C6H5CH3 = H2 + C6H5CH2 */
    rf[423] = rf[423] * c__[1] * c__[52];
    rb[423] = rb[423] * c__[2] * c__[51];
    /*     R424: O2 + C6H5CH3 = HO2 + C6H5CH2 */
    rf[424] = rf[424] * c__[4] * c__[52];
    rb[424] = rb[424] * c__[7] * c__[51];
    /*     R425: OH + C6H5CH3 = H2O + C6H5CH2 */
    rf[425] = rf[425] * c__[5] * c__[52];
    rb[425] = rb[425] * c__[6] * c__[51];
    /*     R426: HO2 + C6H5CH3 = H2O2 + C6H5CH2 */
    rf[426] = rf[426] * c__[7] * c__[52];
    rb[426] = rb[426] * c__[8] * c__[51];
    /*     R427: CH4 + A1- = CH3 + A1 */
    rf[427] = rf[427] * c__[14] * c__[47];
    rb[427] = rb[427] * c__[15] * c__[45];
    /*     R428: H + c-C5H6 = H2 + c-C5H5 */
    rf[428] = rf[428] * c__[1] * c__[46];
    rb[428] = rb[428] * c__[2] * c__[49];
    /*     R429: H + C6H5CH2 = CH3 + A1- */
    rf[429] = rf[429] * c__[1] * c__[51];
    rb[429] = rb[429] * c__[15] * c__[47];
    /*     R430: A1- + C6H5CH3 = A1 + C6H5CH2 */
    rf[430] = rf[430] * c__[47] * c__[52];
    rb[430] = rb[430] * c__[45] * c__[51];
    /*     R431: CH3 + C6H5CH3 = CH4 + C6H5CH2 */
    rf[431] = rf[431] * c__[15] * c__[52];
    rb[431] = rb[431] * c__[14] * c__[51];
    /*     R432: H + c-C5H5 = c-C5H6 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10];
    pr = rklow[37] * ctb / rf[432];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 400.7) * .9319999999999999 + exp(-(*t) / 4135.8) *
        .06800000000000001 + exp(-5501.9 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[432] *= pcor;
    rb[432] *= pcor;
    rf[432] = rf[432] * c__[1] * c__[49];
    rb[432] *= c__[46];
    /*     R433: C2H + A1 = H + A1C2H */
    rf[433] *= c__[45];
    rb[433] = rb[433] * c__[1] * c__[53];
    /*     R434: C2H2 + A1- = H + A1C2H */
    /*     Reaction of PLOG type */
    if (p <= 13629.) {
        r1l = 56.9848088295982 - alogt * 3.38 - 7648.89334812678 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 27154.75) {
        r1l = 56.9848088295982 - alogt * 3.38 - 7648.89334812678 / *t;
        r2l = 61.88211543838745 - alogt * 3.96 - 8605.005016642626 / *t;
        dpl = .6893521095913937;
        p1l = 9.519955154427514;
    } else if (p < 122248.) {
        r1l = 61.88211543838745 - alogt * 3.96 - 8605.005016642626 / *t;
        r2l = 71.37008754696191 - alogt * 5.07 - 10617.87168720231 / *t;
        dpl = 1.504499783266187;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 71.37008754696191 - alogt * 5.07 - 10617.87168720231 / *t;
        r2l = 77.17923053727594 - alogt * 5.7 - 12832.02502481795 / *t;
        dpl = 2.133686556532232;
        p1l = 11.71380704728509;
    } else if (p < 1.0325e7) {
        r1l = 77.17923053727594 - alogt * 5.7 - 12832.02502481795 / *t;
        r2l = 67.69125842870147 - alogt * 4.43 - 13284.92002569388 / *t;
        dpl = 2.302585092994046;
        p1l = 13.84749360381733;
    } else {
        r1l = 67.69125842870147 - alogt * 4.43 - 13284.92002569388 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[434] *= pcor;
    rb[434] *= pcor;
    rf[434] = rf[434] * c__[21] * c__[47];
    rb[434] = rb[434] * c__[1] * c__[53];
    /*     R435: OH + A1C2H = H2O + A1C2H* */
    rf[435] = rf[435] * c__[5] * c__[53];
    rb[435] = rb[435] * c__[6] * c__[55];
    /*     R436: H + A1C2H* = A1C2H */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[38] * ctb / rf[436];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / .1) * 0. + exp(-(*t) / 584.9) * 1. + exp(-6113. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[436] *= pcor;
    rb[436] *= pcor;
    rf[436] = rf[436] * c__[1] * c__[55];
    rb[436] *= c__[53];
    /*     R437: C2H3 + A1 = H + A1C2H3 */
    rf[437] = rf[437] * c__[20] * c__[45];
    rb[437] = rb[437] * c__[1] * c__[56];
    /*     R438: C2H3 + A1- = A1C2H3 */
    /*     Reaction of PLOG type */
    if (p <= 27154.75) {
        r1l = 111.1659383498866 - alogt * 10.52 - 8800.756300354555 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    } else if (p < 122248.) {
        r1l = 111.1659383498866 - alogt * 10.52 - 8800.756300354555 / *t;
        r2l = 88.85921008690934 - alogt * 7.63 - 6475.392079190487 / *t;
        dpl = 1.504499783266187;
        p1l = 10.20930726401891;
    } else if (p < 1032500.) {
        r1l = 88.85921008690934 - alogt * 7.63 - 6475.392079190487 / *t;
        r2l = 62.35211906763319 - alogt * 4.22 - 3640.772590374819 / *t;
        dpl = 2.133686556532232;
        p1l = 11.71380704728509;
    } else {
        r1l = 62.35211906763319 - alogt * 4.22 - 3640.772590374819 / *t;
        r2l = r1l;
        dpl = 1.;
        p1l = pl;
    }
    rl = r1l + (r2l - r1l) / dpl * (pl - p1l);
    pcor = exp(rl);
    rf[438] *= pcor;
    rb[438] *= pcor;
    rf[438] = rf[438] * c__[20] * c__[47];
    rb[438] *= c__[56];
    /*     R439: C2H2 + A1C2H* = A1C2HC2H2 */
    rf[439] = rf[439] * c__[21] * c__[55];
    /*     R440: A1C2HC2H2 = C2H2 + A1C2H* */
    /*     R441: A1C2HC2H2 = A1C2HC2H2u */
    /*     R442: A1C2HC2H2u = A1C2HC2H2 */
    /*     R443: A1C2HC2H2u = A2-1 */
    /*     R444: A2-1 = A1C2HC2H2u */
    rf[444] *= c__[57];
    /*     R445: OH + A2 = H2O + A2-1 */
    rf[445] = rf[445] * c__[5] * c__[58];
    rb[445] = rb[445] * c__[6] * c__[57];
    /*     R446: H + A2-1 = A2 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[39] * ctb / rf[446];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 122.8) * .8 + exp(-(*t) / 478.4) * .2 + exp(-5411.9 /
            *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[446] *= pcor;
    rb[446] *= pcor;
    rf[446] = rf[446] * c__[1] * c__[57];
    rb[446] *= c__[58];
    /*     R447: C2H2 + A2-1 = A2C2H2 */
    rf[447] = rf[447] * c__[21] * c__[57];
    /*     R448: C4H4 + A1- = H + A2 */
    rf[448] = rf[448] * c__[37] * c__[47];
    /*     R449: H + A2R5 = H2 + A2R5- */
    rf[449] = rf[449] * c__[1] * c__[59];
    rb[449] = rb[449] * c__[2] * c__[61];
    /*     R450: OH + A2R5 = H2O + A2R5- */
    rf[450] = rf[450] * c__[5] * c__[59];
    rb[450] = rb[450] * c__[6] * c__[61];
    /*     R451: H + A2R5- = A2R5 */
    ctb = ctot + c__[2] + c__[6] * 5. + c__[14] + c__[9] * .5 + c__[10] + c__[
        17] * 2.;
    pr = rklow[40] * ctb / rf[451];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / .1) * 0. + exp(-(*t) / 584.9) * 1. + exp(-6113. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[451] *= pcor;
    rb[451] *= pcor;
    rf[451] = rf[451] * c__[1] * c__[61];
    rb[451] *= c__[59];
    /*     R452: A2C2H2 = H + A2R5 */
    rb[452] = rb[452] * c__[1] * c__[59];
    /*     R453: C2H2 + A3-4 = H + A4 */
    rf[453] = rf[453] * c__[21] * c__[60];
    /*     R454: OH + A4 = H2O + A4-1 */
    rf[454] = rf[454] * c__[5] * c__[62];
    rb[454] = rb[454] * c__[6] * c__[63];
    /*     R455: H + A4-1 = A4 */
    rf[455] = rf[455] * c__[1] * c__[63];
    rb[455] *= c__[62];
    /*     R456: OH + A4 = H2O + A4-2 */
    rf[456] = rf[456] * c__[5] * c__[62];
    rb[456] = rb[456] * c__[6] * c__[64];
    /*     R457: H + A4-2 = A4 */
    rf[457] = rf[457] * c__[1] * c__[64];
    rb[457] *= c__[62];
    /*     R458: OH + A4 = H2O + A4-4 */
    rf[458] = rf[458] * c__[5] * c__[62];
    rb[458] = rb[458] * c__[6] * c__[65];
    /*     R459: H + A4-4 = A4 */
    rf[459] = rf[459] * c__[1] * c__[65];
    rb[459] *= c__[62];
    /*     R460: O + C6H5O = CO + HCO + 2C2H2 */
    rf[460] = rf[460] * c__[3] * c__[48];
    rb[460] = rb[460] * c__[9] * c__[12] * c__[21] * c__[21];
    /*     R461: OH + A1C2H = CH2CO + A1- */
    rf[461] = rf[461] * c__[5] * c__[53];
    /*     R462: OH + A1C2H = C2H2 + C6H5O */
    rf[462] = rf[462] * c__[5] * c__[53];
    /*     R463: OH + A1C2H3 = C2H4 + C6H5O */
    rf[463] = rf[463] * c__[5] * c__[56];
    /*     R464: OH + A2 = H + CH2CO + A1C2H */
    rf[464] = rf[464] * c__[5] * c__[58];
    /*     R465: OH + A4 = CH2CO + A3-4 */
    rf[465] = rf[465] * c__[5] * c__[62];
    /*     R466: O + A1C2H = HCCO + A1- */
    rf[466] = rf[466] * c__[3] * c__[53];
    /*     R467: O + A1C2H3 = CO + CH3 + A1- */
    rf[467] = rf[467] * c__[3] * c__[56];
    /*     R468: O + A1C2H = C2H + C6H5O */
    rf[468] = rf[468] * c__[3] * c__[53];
    /*     R469: O + A1C2H3 = C2H3 + C6H5O */
    rf[469] = rf[469] * c__[3] * c__[56];
    /*     R470: O + A2 = CH2CO + A1C2H */
    rf[470] = rf[470] * c__[3] * c__[58];
    /*     R471: O + A4 = HCCO + A3-4 */
    rf[471] = rf[471] * c__[3] * c__[62];
    /*     R472: C2H3 + A1C2H = H + A2 */
    rf[472] = rf[472] * c__[20] * c__[53];
    /*     R473: A1C2H = H + A1C2H* */
    rf[473] *= c__[53];
    rb[473] = rb[473] * c__[1] * c__[55];
    /*     R474: CH3 + A1C2H4 = CH4 + A1C2H3 */
    rf[474] *= c__[15];
    rb[474] = rb[474] * c__[14] * c__[56];
    /*     R475: H + A1C2H4 = H2 + A1C2H3 */
    rf[475] *= c__[1];
    rb[475] = rb[475] * c__[2] * c__[56];
    /*     R476: H + A1C2H4 = A1C2H5 */
    ctb = ctot + c__[2] + c__[6] * 11. + c__[10] * 2.6 + c__[9] * .75 + c__[
        14] + c__[17] * 2.;
    pr = rklow[41] * ctb / rf[476];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 369.) * .6850000000000001 + exp(-(*t) / 3285.) * .315
        + exp(-6667. / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[476] *= pcor;
    rb[476] *= pcor;
    rf[476] *= c__[1];
    rb[476] *= c__[66];
    /*     R477: O2 + A1C2H4 = HO2 + A1C2H3 */
    rf[477] *= c__[4];
    rb[477] = rb[477] * c__[7] * c__[56];
    /*     R478: OH + A1C2H4 = H2O + A1C2H3 */
    rf[478] *= c__[5];
    rb[478] = rb[478] * c__[6] * c__[56];
    /*     R479: A1C2H4 = H + A1C2H3 */
    rb[479] = rb[479] * c__[1] * c__[56];
    /*     R480: CH3 + A1C2H5 = CH4 + A1C2H4 */
    rf[480] = rf[480] * c__[15] * c__[66];
    rb[480] *= c__[14];
    /*     R481: H + A1C2H5 = H2 + A1C2H4 */
    rf[481] = rf[481] * c__[1] * c__[66];
    rb[481] *= c__[2];
    /*     R482: H + A1C2H5 = C2H5 + A1 */
    rf[482] = rf[482] * c__[1] * c__[66];
    rb[482] = rb[482] * c__[18] * c__[45];
    /*     R483: HO2 + A1C2H5 = H2O2 + A1C2H4 */
    rf[483] = rf[483] * c__[7] * c__[66];
    rb[483] *= c__[8];
    /*     R484: O + A1C2H5 = OH + A1C2H4 */
    rf[484] = rf[484] * c__[3] * c__[66];
    rb[484] *= c__[5];
    /*     R485: OH + A1C2H5 = H2O + A1C2H4 */
    rf[485] = rf[485] * c__[5] * c__[66];
    rb[485] *= c__[6];
    /*     R486: A1C2H5 = C2H5 + A1- */
    rf[486] *= c__[66];
    rb[486] = rb[486] * c__[18] * c__[47];
    /*     R487: C2H4 + A1C2H* = H + A2 */
    rf[487] = rf[487] * c__[19] * c__[55];
    /*     R488: O + A1 = OH + A1- */
    rf[488] = rf[488] * c__[3] * c__[45];
    rb[488] = rb[488] * c__[5] * c__[47];
    /*     R489: C2H4 + A1- = A1C2H4 */
    rf[489] = rf[489] * c__[19] * c__[47];
    /*     R490: C2H4 + A1- = C2H3 + A1 */
    rf[490] = rf[490] * c__[19] * c__[47];
    rb[490] = rb[490] * c__[20] * c__[45];
    /*     R491: CH2O + A1- = H + A1CHO */
    rf[491] = rf[491] * c__[11] * c__[47];
    rb[491] = rb[491] * c__[1] * c__[54];
    /*     R492: CH3 + A1CHO = CO + CH4 + A1- */
    rf[492] = rf[492] * c__[15] * c__[54];
    /*     R493: H + A1CHO = H2 + CO + A1- */
    rf[493] = rf[493] * c__[1] * c__[54];
    /*     R494: HO2 + A1CHO = H2O2 + CO + A1- */
    rf[494] = rf[494] * c__[7] * c__[54];
    /*     R495: O + A1CHO = OH + CO + A1- */
    rf[495] = rf[495] * c__[3] * c__[54];
    /*     R496: O2 + A1CHO = HO2 + CO + A1- */
    rf[496] = rf[496] * c__[4] * c__[54];
    /*     R497: OH + A1CHO = H2O + CO + A1- */
    rf[497] = rf[497] * c__[5] * c__[54];
    /*     R498: O + A1- = C6H5O */
    rf[498] = rf[498] * c__[3] * c__[47];
    rb[498] *= c__[48];
    /*     R499: C2H3 + A2 = H2 + A2C2H2 */
    rf[499] = rf[499] * c__[20] * c__[58];
    rb[499] *= c__[2];
    /*     R500: H + A2CH3 = CH3 + A2 */
    rf[500] = rf[500] * c__[1] * c__[68];
    rb[500] = rb[500] * c__[15] * c__[58];
    /*     R501: A2CH3 = CH3 + A2-1 */
    rf[501] *= c__[68];
    rb[501] = rb[501] * c__[15] * c__[57];
    /*     R502: A2R5 = H + A2R5- */
    rf[502] *= c__[59];
    rb[502] = rb[502] * c__[1] * c__[61];
    /*     R503: C2H3 + A2-1 = H + A2C2H2 */
    rf[503] = rf[503] * c__[20] * c__[57];
    rb[503] *= c__[1];
    /*     R504: C2H4 + A2-1 = H2 + A2C2H2 */
    rf[504] = rf[504] * c__[19] * c__[57];
    rb[504] *= c__[2];
    /*     R505: A3-4 = C2H2 + A2R5- */
    rf[505] *= c__[60];
    /*     R506: C2H2 + A2R5- = A3-4 */
    rf[506] = rf[506] * c__[21] * c__[61];
    /*     R507: C2H2 + A4-2 = H + A4R5 */
    rf[507] = rf[507] * c__[21] * c__[64];
    rb[507] = rb[507] * c__[1] * c__[71];
    /*     R508: C3H3 + C9H7 = H2 + A2R5 */
    rf[508] = rf[508] * c__[35] * c__[67];
    rb[508] = rb[508] * c__[2] * c__[59];
    /*     R509: HO2 + C9H7 = H2O + C9H6O */
    rf[509] = rf[509] * c__[7] * c__[67];
    rb[509] = rb[509] * c__[6] * c__[69];
    /*     R510: HO2 + C9H7 = H + OH + C9H6O */
    rf[510] = rf[510] * c__[7] * c__[67];
    /*     R511: O + C9H7 = H + C9H6O */
    rf[511] = rf[511] * c__[3] * c__[67];
    rb[511] = rb[511] * c__[1] * c__[69];
    /*     R512: O2 + C9H7 = OH + C9H6O */
    rf[512] = rf[512] * c__[4] * c__[67];
    rb[512] = rb[512] * c__[5] * c__[69];
    /*     R513: O2 + C9H7 = H + O + C9H6O */
    rf[513] = rf[513] * c__[4] * c__[67];
    /*     R514: HCO + C9H8 = CH2O + C9H7 */
    rf[514] = rf[514] * c__[12] * c__[70];
    rb[514] = rb[514] * c__[11] * c__[67];
    /*     R515: HO2 + C9H8 = H2O2 + C9H7 */
    rf[515] = rf[515] * c__[7] * c__[70];
    rb[515] = rb[515] * c__[8] * c__[67];
    /*     R516: O + C9H8 = OH + C9H7 */
    rf[516] = rf[516] * c__[3] * c__[70];
    rb[516] = rb[516] * c__[5] * c__[67];
    /*     R517: O + C9H8 = 2H + C9H6O */
    rf[517] = rf[517] * c__[3] * c__[70];
    /*     R518: O2 + C9H8 = HO2 + C9H7 */
    rf[518] = rf[518] * c__[4] * c__[70];
    rb[518] = rb[518] * c__[7] * c__[67];
    /*     R519: OH + C9H8 = H2O + C9H7 */
    rf[519] = rf[519] * c__[5] * c__[70];
    rb[519] = rb[519] * c__[6] * c__[67];
    /*     R520: C9H8 = H + C9H7 */
    rf[520] *= c__[70];
    rb[520] = rb[520] * c__[1] * c__[67];
    /*     R521: C2H2 + C6H5CH2 = H + C9H8 */
    rf[521] = rf[521] * c__[21] * c__[51];
    rb[521] = rb[521] * c__[1] * c__[70];
    /*     R522: C6H5CH2 = C2H2 + c-C5H5 */
    rf[522] *= c__[51];
    rb[522] = rb[522] * c__[21] * c__[49];
    /*     R523: O + A1C2H4 = CH2O + C6H5CH2 */
    rf[523] *= c__[3];
    rb[523] = rb[523] * c__[11] * c__[51];
    /*     R524: A1C2H5 = CH3 + C6H5CH2 */
    ctb = ctot + c__[2] + c__[6] * 11. + c__[10] * 2.6 + c__[9] * .75 + c__[
        14] + c__[17] * 2.;
    pr = rklow[42] * ctb / rf[524];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,1e-200);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 50.) * .69 + exp(-(*t) / 3e3) * .31 + exp(-9e3 / *t);
    d__1 = max(fcent,1e-200);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
    /* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b8, flog);
    pcor = fc * pcor;
    rf[524] *= pcor;
    rb[524] *= pcor;
    rf[524] *= c__[66];
    rb[524] = rb[524] * c__[15] * c__[51];
    /*     R525: 2c-C5H5 = 2H + A2 */
    rf[525] = rf[525] * c__[49] * c__[49];
    /*     R526: OH + c-C5H6 = HCO + C4H6 */
    rf[526] = rf[526] * c__[5] * c__[46];
    rb[526] = rb[526] * c__[12] * c__[36];
    /*     R527: HO2 + c-C5H5 = H2O + C5H4O */
    rf[527] = rf[527] * c__[7] * c__[49];
    rb[527] = rb[527] * c__[6] * c__[50];
    /*     R528: c-C5H5 + C6H5CH3 = c-C5H6 + C6H5CH2 */
    rf[528] = rf[528] * c__[49] * c__[52];
    rb[528] = rb[528] * c__[46] * c__[51];
    /*     R529: c-C5H6 + A1- = A1 + c-C5H5 */
    rf[529] = rf[529] * c__[46] * c__[47];
    rb[529] = rb[529] * c__[45] * c__[49];
    /*     R530: C2H3 + c-C5H6 = C2H4 + c-C5H5 */
    rf[530] = rf[530] * c__[20] * c__[46];
    rb[530] = rb[530] * c__[19] * c__[49];
    /*     R531: c-C5H5 = C2H2 + C3H3 */
    ctb = ctot;
    pr = rklow[43] * ctb / rf[531];
    pcor = pr / (pr + 1.f);
    rf[531] *= pcor;
    rb[531] *= pcor;
    rf[531] *= c__[49];
    rb[531] = rb[531] * c__[21] * c__[35];
    /*     R532: C2H3 + C4H6 = CH3 + c-C5H6 */
    rf[532] = rf[532] * c__[20] * c__[36];
    rb[532] = rb[532] * c__[15] * c__[46];
    /*     R533: CH2 + C6H5CH2 = H + A1C2H3 */
    rf[533] = rf[533] * c__[16] * c__[51];
    rb[533] = rb[533] * c__[1] * c__[56];
    /*     R534: CH3 + C6H5CH2 = A1C2H5 */
    rf[534] = rf[534] * c__[15] * c__[51];
    rb[534] *= c__[66];
    /*     R535: C4H4 + A1- = C2H + A1C2H3 */
    rf[535] = rf[535] * c__[37] * c__[47];
    rb[535] *= c__[56];
    /*     R536: C4H6 + A1- = C2H3 + A1C2H3 */
    rf[536] = rf[536] * c__[36] * c__[47];
    rb[536] = rb[536] * c__[20] * c__[56];
    /*     R537: 2C4H4 = A1C2H3 */
    rf[537] = rf[537] * c__[37] * c__[37];
    rb[537] *= c__[56];
    /*     R538: A1C2H5 = H2 + A1C2H3 */
    rf[538] = rf[538] * c__[1] * c__[66];
    rb[538] = rb[538] * c__[1] * c__[2] * c__[56];
    /*     R539: A1C2H5 = H2 + A1C2H3 */
    rf[539] *= c__[66];
    rb[539] = rb[539] * c__[2] * c__[56];
    /*     R540: C2H6 + A1- = H + A1C2H5 */
    rf[540] = rf[540] * c__[17] * c__[47];
    rb[540] = rb[540] * c__[1] * c__[66];
    /*     R541: C4H3-N + A1- = A2 */
    rf[541] = rf[541] * c__[39] * c__[47];
    rb[541] *= c__[58];
    /*     R542: C4H3-N + A1- = H + A2-1 */
    rf[542] = rf[542] * c__[39] * c__[47];
    rb[542] = rb[542] * c__[1] * c__[57];
    /*     R543: C3H3 + A1- = C9H8 */
    rf[543] = rf[543] * c__[35] * c__[47];
    rb[543] *= c__[70];
    /*     R544: C3H3 + A1 = H + C9H8 */
    rf[544] = rf[544] * c__[35] * c__[45];
    rb[544] = rb[544] * c__[1] * c__[70];
    /*     R545: CH2 + C6H5CH3 = H2 + A1C2H3 */
    rf[545] = rf[545] * c__[16] * c__[52];
    rb[545] = rb[545] * c__[2] * c__[56];
    /*     R546: H + C6H5CH3 = CH4 + A1- */
    rf[546] = rf[546] * c__[1] * c__[52];
    rb[546] = rb[546] * c__[14] * c__[47];
    /*     R547: O + C6H5CH3 = OH + C6H5CH2 */
    rf[547] = rf[547] * c__[3] * c__[52];
    rb[547] = rb[547] * c__[5] * c__[51];
    /*     R548: C6H5CH2 + C9H7 = 2H2 + A4 */
    rf[548] = rf[548] * c__[51] * c__[67];
    /*     R549: 2C9H7 = H2 + C2H2 + A4 */
    rf[549] = rf[549] * c__[67] * c__[67];
    /*     R550: A1C2H + A1C2H* = H + A4 */
    rf[550] = rf[550] * c__[53] * c__[55];
    /*     R551: CH3 + A2 = CH4 + A2-1 */
    rf[551] = rf[551] * c__[15] * c__[58];
    /*     R552: CH4 + A2-1 = CH3 + A2 */
    rf[552] = rf[552] * c__[14] * c__[57];
    /*     R553: CH3 + A4 = CH4 + A4-1 */
    rf[553] = rf[553] * c__[15] * c__[62];
    /*     R554: CH4 + A4-1 = CH3 + A4 */
    rf[554] = rf[554] * c__[14] * c__[63];
    /*     R555: CH3 + A4 = CH4 + A4-2 */
    rf[555] = rf[555] * c__[15] * c__[62];
    /*     R556: CH4 + A4-2 = CH3 + A4 */
    rf[556] = rf[556] * c__[14] * c__[64];
    /*     R557: CH3 + A4 = CH4 + A4-4 */
    rf[557] = rf[557] * c__[15] * c__[62];
    /*     R558: CH4 + A4-4 = CH3 + A4 */
    rf[558] = rf[558] * c__[14] * c__[65];
    /*     R559: H + A2 = H2 + A2-1 */
    rf[559] = rf[559] * c__[1] * c__[58];
    /*     R560: H2 + A2-1 = H + A2 */
    rf[560] = rf[560] * c__[2] * c__[57];
    /*     R561: H + A4 = H2 + A4-1 */
    rf[561] = rf[561] * c__[1] * c__[62];
    /*     R562: H2 + A4-1 = H + A4 */
    rf[562] = rf[562] * c__[2] * c__[63];
    /*     R563: H + A4 = H2 + A4-2 */
    rf[563] = rf[563] * c__[1] * c__[62];
    /*     R564: H2 + A4-2 = H + A4 */
    rf[564] = rf[564] * c__[2] * c__[64];
    /*     R565: H + A4 = H2 + A4-4 */
    rf[565] = rf[565] * c__[1] * c__[62];
    /*     R566: H2 + A4-4 = H + A4 */
    rf[566] = rf[566] * c__[2] * c__[65];
    /*     R567: H + C9H8 = H2 + C9H7 */
    rf[567] = rf[567] * c__[1] * c__[70];
    /*     R568: H2 + C9H7 = H + C9H8 */
    rf[568] = rf[568] * c__[2] * c__[67];
    /*     R569: CH3 + C9H8 = CH4 + C9H7 */
    rf[569] = rf[569] * c__[15] * c__[70];
    /*     R570: CH4 + C9H7 = CH3 + C9H8 */
    rf[570] = rf[570] * c__[14] * c__[67];
    /*     R571: H + A1C2H = H2 + A1C2H* */
    rf[571] = rf[571] * c__[1] * c__[53];
    /*     R572: H2 + A1C2H* = H + A1C2H */
    rf[572] = rf[572] * c__[2] * c__[55];
    /*     R573: CH3 + A1C2H = CH4 + A1C2H* */
    rf[573] = rf[573] * c__[15] * c__[53];
    /*     R574: CH4 + A1C2H* = CH3 + A1C2H */
    rf[574] = rf[574] * c__[14] * c__[55];
    /*     R575: C4H4 + A4-1 = H + BAPYR */
    rf[575] = rf[575] * c__[37] * c__[63];
    /*     R576: C4H4 + A4-2 = H + BAPYR */
    rf[576] = rf[576] * c__[37] * c__[64];
    /*     R577: C4H4 + A4-4 = H + BEPYREN */
    rf[577] = rf[577] * c__[37] * c__[65];
    /*     R578: C2H2 + A4-1 = H + PYC2H-1 */
    rf[578] = rf[578] * c__[21] * c__[63];
    rb[578] = rb[578] * c__[1] * c__[74];
    /*     R579: C2H2 + A4-2 = H + PYC2H-2 */
    rf[579] = rf[579] * c__[21] * c__[64];
    rb[579] = rb[579] * c__[1] * c__[75];
    /*     R580: C2H2 + A4-4 = H + PYC2H-4 */
    rf[580] = rf[580] * c__[21] * c__[65];
    rb[580] = rb[580] * c__[1] * c__[76];
    /*     R581: H + PYC2H-1 = H2 + PYC2H-1JP */
    rf[581] = rf[581] * c__[1] * c__[74];
    /*     R582: H + PYC2H-2 = H2 + PYC2H-2JS */
    rf[582] = rf[582] * c__[1] * c__[75];
    /*     R583: H2 + PYC2H-1JP = H + PYC2H-1 */
    rf[583] = rf[583] * c__[2] * c__[77];
    /*     R584: H2 + PYC2H-2JS = H + PYC2H-2 */
    rf[584] = rf[584] * c__[2] * c__[78];
    /*     R585: CH3 + PYC2H-1 = CH4 + PYC2H-1JP */
    rf[585] = rf[585] * c__[15] * c__[74];
    /*     R586: CH3 + PYC2H-2 = CH4 + PYC2H-2JS */
    rf[586] = rf[586] * c__[15] * c__[75];
    /*     R587: CH4 + PYC2H-1JP = CH3 + PYC2H-1 */
    rf[587] = rf[587] * c__[14] * c__[77];
    /*     R588: CH4 + PYC2H-2JS = CH3 + PYC2H-2 */
    rf[588] = rf[588] * c__[14] * c__[78];
    /*     R589: OH + PYC2H-1 = H2O + PYC2H-1JP */
    rf[589] = rf[589] * c__[5] * c__[74];
    rb[589] = rb[589] * c__[6] * c__[77];
    /*     R590: OH + PYC2H-2 = H2O + PYC2H-2JS */
    rf[590] = rf[590] * c__[5] * c__[75];
    rb[590] = rb[590] * c__[6] * c__[78];
    /*     R591: C2H2 + PYC2H-1JP = BAPYRJS */
    rf[591] = rf[591] * c__[21] * c__[77];
    /*     R592: C2H2 + PYC2H-2JS = BAPYRJS */
    rf[592] = rf[592] * c__[21] * c__[78];
    /*     R593: H + BAPYRJS = BAPYR */
    rf[593] = rf[593] * c__[1] * c__[79];
    rb[593] *= c__[72];
    /*     R594: OH + BAPYR = H + CH2CO + PYC2H-2 */
    rf[594] = rf[594] * c__[5] * c__[72];
    /*     R595: O + BAPYR = CH2CO + PYC2H-2 */
    rf[595] = rf[595] * c__[3] * c__[72];
    rb[595] = rb[595] * c__[25] * c__[75];
    /*     R596: O2 + BAPYRJS = CO + HCO + PYC2H-2 */
    rf[596] = rf[596] * c__[4] * c__[79];
    /*     R597: OH + PYC2H-2 = CH2CO + A4-2 */
    rf[597] = rf[597] * c__[5] * c__[75];
    rb[597] = rb[597] * c__[25] * c__[64];
    /*     R598: O + PYC2H-2 = HCCO + A4-2 */
    rf[598] = rf[598] * c__[3] * c__[75];
    rb[598] = rb[598] * c__[26] * c__[64];
    /*     R599: H + PYC2H-4 = H2 + PYC2H-4JS */
    rf[599] = rf[599] * c__[1] * c__[76];
    /*     R600: H2 + PYC2H-4JS = H + PYC2H-4 */
    rf[600] = rf[600] * c__[2] * c__[80];
    /*     R601: CH3 + PYC2H-4 = CH4 + PYC2H-4JS */
    rf[601] = rf[601] * c__[15] * c__[76];
    /*     R602: CH4 + PYC2H-4JS = CH3 + PYC2H-4 */
    rf[602] = rf[602] * c__[14] * c__[80];
    /*     R603: OH + PYC2H-4 = H2O + PYC2H-4JS */
    rf[603] = rf[603] * c__[5] * c__[76];
    rb[603] = rb[603] * c__[6] * c__[80];
    /*     R604: C2H2 + PYC2H-4JS = BEPYRENJS */
    rf[604] = rf[604] * c__[21] * c__[80];
    /*     R605: H + BEPYRENJS = BEPYREN */
    rf[605] = rf[605] * c__[1] * c__[81];
    rb[605] *= c__[73];
    /*     R606: H + BEPYREN = H2 + BEPYRENJS */
    rf[606] = rf[606] * c__[1] * c__[73];
    /*     R607: H2 + BEPYRENJS = H + BEPYREN */
    rf[607] = rf[607] * c__[2] * c__[81];
    /*     R608: CH3 + BEPYREN = CH4 + BEPYRENJS */
    rf[608] = rf[608] * c__[15] * c__[73];
    /*     R609: CH4 + BEPYRENJS = CH3 + BEPYREN */
    rf[609] = rf[609] * c__[14] * c__[81];
    /*     R610: OH + BEPYREN = H2O + BEPYRENJS */
    rf[610] = rf[610] * c__[5] * c__[73];
    rb[610] = rb[610] * c__[6] * c__[81];
    /*     R611: C2H2 + BEPYRENJS = H + BGHIPER */
    rf[611] = rf[611] * c__[21] * c__[81];
    /*     R612: H + BAPYR = H2 + BAPYRJS */
    rf[612] = rf[612] * c__[1] * c__[72];
    /*     R613: H2 + BAPYRJS = H + BAPYR */
    rf[613] = rf[613] * c__[2] * c__[79];
    /*     R614: CH3 + BAPYR = CH4 + BAPYRJS */
    rf[614] = rf[614] * c__[15] * c__[72];
    /*     R615: CH4 + BAPYRJS = CH3 + BAPYR */
    rf[615] = rf[615] * c__[14] * c__[79];
    /*     R616: OH + BAPYR = H2O + BAPYRJS */
    rf[616] = rf[616] * c__[5] * c__[72];
    rb[616] = rb[616] * c__[6] * c__[79];
    /*     R617: C2H2 + BAPYRJS = H + ANTHAN */
    rf[617] = rf[617] * c__[21] * c__[79];
    /*     R618: OH + ANTHAN = CH2CO + BAPYRJS */
    rf[618] = rf[618] * c__[5] * c__[83];
    /*     R619: O + ANTHAN = HCCO + BAPYRJS */
    rf[619] = rf[619] * c__[3] * c__[83];
    /*     R620: H + BGHIPER = H2 + BGHIPEJS1 */
    rf[620] = rf[620] * c__[1] * c__[82];
    /*     R621: H2 + BGHIPEJS1 = H + BGHIPER */
    rf[621] = rf[621] * c__[2] * c__[84];
    /*     R622: CH3 + BGHIPER = CH4 + BGHIPEJS1 */
    rf[622] = rf[622] * c__[15] * c__[82];
    /*     R623: CH4 + BGHIPEJS1 = CH3 + BGHIPER */
    rf[623] = rf[623] * c__[14] * c__[84];
    /*     R624: OH + BGHIPER = H2O + BGHIPEJS1 */
    rf[624] = rf[624] * c__[5] * c__[82];
    rb[624] = rb[624] * c__[6] * c__[84];
    /*     R625: C2H2 + BGHIPEJS1 = H + CORONEN */
    rf[625] = rf[625] * c__[21] * c__[84];

    return 0;
} /* ratx_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int rdot_(double *rf, double *rb, double *wdot)
{
    static long int k;
    static double rop;


    /* Parameter adjustments */
    --wdot;
    --rb;
    --rf;

    /* Function Body */
    for (k = 1; k <= 86; ++k) {
        wdot[k] = 0.;
    }

    rop = rf[1] - rb[1];
    wdot[1] -= rop;
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[5] += rop;
    rop = rf[2] - rb[2];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[3] -= rop;
    wdot[5] += rop;
    rop = rf[3] - rb[3];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[5] -= rop;
    wdot[6] += rop;
    rop = rf[4] - rb[4];
    wdot[3] -= rop;
    wdot[5] = wdot[5] + rop + rop;
    wdot[6] -= rop;
    rop = rf[5] - rb[5];
    wdot[1] = wdot[1] + rop + rop;
    wdot[2] -= rop;
    rop = rf[6] - rb[6];
    wdot[3] = wdot[3] - rop - rop;
    wdot[4] += rop;
    rop = rf[7] - rb[7];
    wdot[1] -= rop;
    wdot[3] -= rop;
    wdot[5] += rop;
    rop = rf[8] - rb[8];
    wdot[1] -= rop;
    wdot[5] -= rop;
    wdot[6] += rop;
    rop = rf[9] - rb[9];
    wdot[1] -= rop;
    wdot[4] -= rop;
    wdot[7] += rop;
    rop = rf[10] - rb[10];
    wdot[1] -= rop;
    wdot[5] = wdot[5] + rop + rop;
    wdot[7] -= rop;
    rop = rf[11] - rb[11];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[4] -= rop;
    wdot[7] += rop;
    rop = rf[12] - rb[12];
    wdot[3] -= rop;
    wdot[4] += rop;
    wdot[5] += rop;
    wdot[7] -= rop;
    rop = rf[13] - rb[13];
    wdot[4] += rop;
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[7] -= rop;
    rop = rf[14] - rb[14];
    wdot[4] += rop;
    wdot[7] = wdot[7] - rop - rop;
    wdot[8] += rop;
    rop = rf[15] - rb[15];
    wdot[4] += rop;
    wdot[7] = wdot[7] - rop - rop;
    wdot[8] += rop;
    rop = rf[16] - rb[16];
    wdot[5] = wdot[5] + rop + rop;
    wdot[8] -= rop;
    rop = rf[17] - rb[17];
    wdot[5] = wdot[5] + rop + rop;
    wdot[8] -= rop;
    rop = rf[18] - rb[18];
    wdot[1] -= rop;
    wdot[5] += rop;
    wdot[6] += rop;
    wdot[8] -= rop;
    rop = rf[19] - rb[19];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[7] += rop;
    wdot[8] -= rop;
    rop = rf[20] - rb[20];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[7] += rop;
    wdot[8] -= rop;
    rop = rf[21] - rb[21];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[7] += rop;
    wdot[8] -= rop;
    rop = rf[22] - rb[22];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[7] += rop;
    wdot[8] -= rop;
    rop = rf[23] - rb[23];
    wdot[3] -= rop;
    wdot[9] -= rop;
    wdot[10] += rop;
    rop = rf[24] - rb[24];
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[9] -= rop;
    wdot[10] += rop;
    rop = rf[25] - rb[25];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[9] -= rop;
    wdot[10] += rop;
    rop = rf[26] - rb[26];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[9] -= rop;
    wdot[10] += rop;
    rop = rf[27] - rb[27];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[9] -= rop;
    wdot[10] += rop;
    rop = rf[28] - rb[28];
    wdot[1] += rop;
    wdot[9] += rop;
    wdot[12] -= rop;
    rop = rf[29] - rb[29];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[9] += rop;
    wdot[12] -= rop;
    rop = rf[30] - rb[30];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[9] += rop;
    wdot[12] -= rop;
    rop = rf[31] - rb[31];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[9] += rop;
    wdot[12] -= rop;
    rop = rf[32] - rb[32];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[10] += rop;
    wdot[12] -= rop;
    rop = rf[33] - rb[33];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[9] += rop;
    wdot[12] -= rop;
    rop = rf[34];
    wdot[1] += rop;
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[10] += rop;
    wdot[12] -= rop;
    rop = rf[35];
    wdot[2] += rop;
    wdot[9] = wdot[9] + rop + rop;
    wdot[12] = wdot[12] - rop - rop;
    rop = rf[36] - rb[36];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[37] - rb[37];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[11] -= rop;
    wdot[12] += rop;
    rop = rf[38] - rb[38];
    wdot[9] += rop;
    wdot[11] += rop;
    wdot[12] = wdot[12] - rop - rop;
    rop = rf[39] - rb[39];
    wdot[1] -= rop;
    wdot[11] += rop;
    wdot[12] -= rop;
    rop = rf[40] - rb[40];
    wdot[2] -= rop;
    wdot[9] -= rop;
    wdot[11] += rop;
    rop = rf[41] - rb[41];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[11] -= rop;
    wdot[12] += rop;
    rop = rf[42] - rb[42];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[11] -= rop;
    wdot[12] += rop;
    rop = rf[43] - rb[43];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[11] -= rop;
    wdot[12] += rop;
    rop = rf[44] - rb[44];
    wdot[11] -= rop;
    wdot[12] += rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[45] - rb[45];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[11] -= rop;
    wdot[12] += rop;
    rop = rf[46] - rb[46];
    wdot[1] += rop;
    wdot[11] += rop;
    rop = rf[47] - rb[47];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[11] += rop;
    rop = rf[48] - rb[48];
    wdot[11] -= rop;
    wdot[12] += rop;
    wdot[13] += rop;
    rop = rf[49] - rb[49];
    wdot[13] -= rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[50] - rb[50];
    wdot[11] += rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[51] - rb[51];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[11] += rop;
    rop = rf[52] - rb[52];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[11] += rop;
    rop = rf[53] - rb[53];
    wdot[1] -= rop;
    wdot[11] -= rop;
    rop = rf[54] - rb[54];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[11] += rop;
    rop = rf[55] - rb[55];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[11] += rop;
    rop = rf[56] - rb[56];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[11] += rop;
    rop = rf[57] - rb[57];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[11] += rop;
    rop = rf[58] - rb[58];
    wdot[11] = wdot[11] + rop + rop;
    wdot[12] -= rop;
    rop = rf[59] - rb[59];
    wdot[11] += rop;
    wdot[13] += rop;
    rop = rf[60] - rb[60];
    wdot[11] += rop;
    wdot[12] -= rop;
    wdot[13] -= rop;
    rop = rf[61] - rb[61];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[11] += rop;
    rop = rf[62] - rb[62];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[11] += rop;
    rop = rf[63] - rb[63];
    wdot[11] += rop;
    wdot[13] += rop;
    rop = rf[64] - rb[64];
    wdot[5] += rop;
    wdot[13] -= rop;
    wdot[15] += rop;
    rop = rf[65] - rb[65];
    wdot[6] += rop;
    wdot[13] -= rop;
    rop = rf[66] - rb[66];
    wdot[1] += rop;
    wdot[13] -= rop;
    rop = rf[67] - rb[67];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[13] -= rop;
    rop = rf[68] - rb[68];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[13] -= rop;
    rop = rf[69] - rb[69];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[13] -= rop;
    rop = rf[70] - rb[70];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[13] -= rop;
    rop = rf[71] - rb[71];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[13] -= rop;
    rop = rf[72] - rb[72];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[13] -= rop;
    rop = rf[73] - rb[73];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[13] -= rop;
    rop = rf[74] - rb[74];
    wdot[13] -= rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[75] - rb[75];
    rop = rf[76] - rb[76];
    wdot[11] += rop;
    wdot[13] += rop;
    rop = rf[77] - rb[77];
    wdot[1] -= rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[78] - rb[78];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[14] -= rop;
    wdot[15] += rop;
    rop = rf[79] - rb[79];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[14] -= rop;
    wdot[15] += rop;
    rop = rf[80] - rb[80];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[14] -= rop;
    wdot[15] += rop;
    rop = rf[81] - rb[81];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[14] -= rop;
    wdot[15] += rop;
    rop = rf[82] - rb[82];
    wdot[14] -= rop;
    wdot[15] = wdot[15] + rop + rop;
    wdot[16] -= rop;
    rop = rf[83] - rb[83];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[15] -= rop;
    rop = rf[84] - rb[84];
    wdot[2] += rop;
    wdot[5] -= rop;
    wdot[11] += rop;
    wdot[15] -= rop;
    rop = rf[85] - rb[85];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[15] -= rop;
    rop = rf[86] - rb[86];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[15] -= rop;
    rop = rf[87] - rb[87];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[15] -= rop;
    rop = rf[88] - rb[88];
    wdot[4] += rop;
    wdot[7] -= rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    rop = rf[89] - rb[89];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[11] += rop;
    wdot[15] -= rop;
    rop = rf[90] - rb[90];
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[15] -= rop;
    rop = rf[91] - rb[91];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[11] += rop;
    wdot[15] -= rop;
    rop = rf[92] - rb[92];
    wdot[16] += rop;
    rop = rf[93] - rb[93];
    wdot[2] += rop;
    wdot[3] -= rop;
    wdot[9] += rop;
    rop = rf[94] - rb[94];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[12] += rop;
    rop = rf[95] - rb[95];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[11] += rop;
    rop = rf[96] - rb[96];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[15] += rop;
    rop = rf[97];
    wdot[1] += rop;
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[9] += rop;
    rop = rf[98] - rb[98];
    wdot[4] -= rop;
    wdot[6] += rop;
    wdot[9] += rop;
    rop = rf[99] - rb[99];
    wdot[16] += rop;
    rop = rf[100] - rb[100];
    wdot[16] += rop;
    rop = rf[101] - rb[101];
    wdot[16] += rop;
    rop = rf[102] - rb[102];
    wdot[9] += rop;
    wdot[10] -= rop;
    wdot[11] += rop;
    rop = rf[103] - rb[103];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[16] -= rop;
    rop = rf[104] - rb[104];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[12] += rop;
    wdot[16] -= rop;
    rop = rf[105];
    wdot[1] = wdot[1] + rop + rop;
    wdot[4] -= rop;
    wdot[10] += rop;
    wdot[16] -= rop;
    rop = rf[106];
    wdot[1] = wdot[1] + rop + rop;
    wdot[3] -= rop;
    wdot[9] += rop;
    wdot[16] -= rop;
    rop = rf[107] - rb[107];
    wdot[15] = wdot[15] - rop - rop;
    wdot[17] += rop;
    rop = rf[108] - rb[108];
    wdot[1] -= rop;
    wdot[17] += rop;
    wdot[18] -= rop;
    rop = rf[109] - rb[109];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[110] - rb[110];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[111] - rb[111];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[112] - rb[112];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[113] - rb[113];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[114] - rb[114];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[115] - rb[115];
    wdot[13] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[116] - rb[116];
    wdot[15] += rop;
    wdot[17] -= rop;
    wdot[18] += rop;
    rop = rf[117] - rb[117];
    wdot[1] -= rop;
    wdot[18] += rop;
    wdot[19] -= rop;
    rop = rf[118] - rb[118];
    wdot[18] += rop;
    wdot[19] = wdot[19] - rop - rop;
    wdot[20] += rop;
    rop = rf[119] - rb[119];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[18] -= rop;
    wdot[19] += rop;
    rop = rf[120] - rb[120];
    wdot[1] += rop;
    wdot[15] = wdot[15] - rop - rop;
    wdot[18] += rop;
    rop = rf[121] - rb[121];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[18] -= rop;
    wdot[19] += rop;
    rop = rf[122] - rb[122];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[18] -= rop;
    wdot[22] += rop;
    rop = rf[123] - rb[123];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[18] -= rop;
    wdot[19] += rop;
    rop = rf[124] - rb[124];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[18] -= rop;
    wdot[19] += rop;
    rop = rf[125] - rb[125];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[18] -= rop;
    wdot[27] += rop;
    rop = rf[126] - rb[126];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[18] -= rop;
    wdot[22] += rop;
    rop = rf[127] - rb[127];
    wdot[12] += rop;
    wdot[15] += rop;
    wdot[27] -= rop;
    rop = rf[128] - rb[128];
    wdot[22] += rop;
    wdot[27] -= rop;
    rop = rf[129] - rb[129];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[27] -= rop;
    rop = rf[130] - rb[130];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[27] -= rop;
    rop = rf[131] - rb[131];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[27] -= rop;
    rop = rf[132] - rb[132];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[27] -= rop;
    rop = rf[133] - rb[133];
    wdot[13] += rop;
    wdot[27] -= rop;
    rop = rf[134] - rb[134];
    rop = rf[135] - rb[135];
    wdot[24] += rop;
    rop = rf[136] - rb[136];
    wdot[12] += rop;
    wdot[15] += rop;
    wdot[22] -= rop;
    rop = rf[137] - rb[137];
    wdot[9] += rop;
    wdot[14] += rop;
    wdot[22] -= rop;
    rop = rf[138] - rb[138];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[22] -= rop;
    rop = rf[139] - rb[139];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[22] -= rop;
    wdot[24] += rop;
    rop = rf[140] - rb[140];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[22] -= rop;
    rop = rf[141] - rb[141];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[22] -= rop;
    rop = rf[142] - rb[142];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[22] -= rop;
    rop = rf[143] - rb[143];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[22] -= rop;
    rop = rf[144] - rb[144];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[22] -= rop;
    rop = rf[145] - rb[145];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[22] -= rop;
    wdot[24] += rop;
    rop = rf[146] - rb[146];
    wdot[9] += rop;
    wdot[15] += rop;
    rop = rf[147] - rb[147];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[25] += rop;
    rop = rf[148] - rb[148];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[25] += rop;
    rop = rf[149] - rb[149];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[25] += rop;
    rop = rf[150] - rb[150];
    wdot[1] += rop;
    wdot[24] -= rop;
    wdot[25] += rop;
    rop = rf[151] - rb[151];
    wdot[9] += rop;
    wdot[15] += rop;
    wdot[24] -= rop;
    rop = rf[152] - rb[152];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[24] -= rop;
    wdot[25] += rop;
    rop = rf[153];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[9] += rop;
    wdot[11] += rop;
    wdot[24] -= rop;
    rop = rf[154] - rb[154];
    wdot[9] -= rop;
    wdot[16] -= rop;
    wdot[25] += rop;
    rop = rf[155] - rb[155];
    wdot[1] += rop;
    wdot[25] += rop;
    rop = rf[156] - rb[156];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[25] -= rop;
    wdot[26] += rop;
    rop = rf[157] - rb[157];
    wdot[1] -= rop;
    wdot[9] += rop;
    wdot[15] += rop;
    wdot[25] -= rop;
    rop = rf[158] - rb[158];
    wdot[3] -= rop;
    wdot[10] += rop;
    wdot[16] += rop;
    wdot[25] -= rop;
    rop = rf[159] - rb[159];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[25] -= rop;
    wdot[26] += rop;
    rop = rf[160] - rb[160];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[25] -= rop;
    wdot[26] += rop;
    rop = rf[161] - rb[161];
    wdot[5] -= rop;
    wdot[9] += rop;
    wdot[25] -= rop;
    rop = rf[162] - rb[162];
    wdot[9] += rop;
    wdot[15] -= rop;
    wdot[18] += rop;
    wdot[25] -= rop;
    rop = rf[163] - rb[163];
    wdot[9] += rop;
    wdot[19] += rop;
    wdot[25] -= rop;
    rop = rf[164];
    wdot[2] += rop;
    wdot[5] -= rop;
    wdot[9] = wdot[9] + rop + rop;
    wdot[26] -= rop;
    rop = rf[165];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[9] = wdot[9] + rop + rop;
    wdot[26] -= rop;
    rop = rf[166] - rb[166];
    wdot[1] -= rop;
    wdot[9] += rop;
    wdot[26] -= rop;
    rop = rf[167];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[9] = wdot[9] + rop + rop;
    wdot[26] -= rop;
    rop = rf[168];
    wdot[1] += rop;
    wdot[4] -= rop;
    wdot[9] += rop;
    wdot[10] += rop;
    wdot[26] -= rop;
    rop = rf[169] - rb[169];
    wdot[1] -= rop;
    wdot[19] += rop;
    wdot[20] -= rop;
    rop = rf[170] - rb[170];
    wdot[2] += rop;
    wdot[19] -= rop;
    rop = rf[171] - rb[171];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[19] -= rop;
    wdot[20] += rop;
    rop = rf[172] - rb[172];
    wdot[3] -= rop;
    wdot[12] += rop;
    wdot[15] += rop;
    wdot[19] -= rop;
    rop = rf[173] - rb[173];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[19] -= rop;
    wdot[24] += rop;
    rop = rf[174] - rb[174];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[19] -= rop;
    wdot[20] += rop;
    rop = rf[175] - rb[175];
    wdot[5] -= rop;
    wdot[11] += rop;
    wdot[15] += rop;
    wdot[19] -= rop;
    rop = rf[176] - rb[176];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[19] -= rop;
    wdot[22] += rop;
    rop = rf[177] - rb[177];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[19] -= rop;
    wdot[23] += rop;
    rop = rf[178] - rb[178];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[23] -= rop;
    wdot[24] += rop;
    rop = rf[179] - rb[179];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[23] -= rop;
    wdot[24] += rop;
    rop = rf[180] - rb[180];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[23] -= rop;
    wdot[24] += rop;
    rop = rf[181] - rb[181];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[23] -= rop;
    wdot[24] += rop;
    rop = rf[182] - rb[182];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[23] -= rop;
    wdot[24] += rop;
    rop = rf[183] - rb[183];
    wdot[22] += rop;
    wdot[23] -= rop;
    rop = rf[184] - rb[184];
    wdot[22] += rop;
    wdot[23] -= rop;
    rop = rf[185] - rb[185];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[19] -= rop;
    wdot[20] += rop;
    rop = rf[186] - rb[186];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[19] -= rop;
    wdot[20] += rop;
    rop = rf[187] - rb[187];
    wdot[13] += rop;
    wdot[19] -= rop;
    wdot[20] += rop;
    rop = rf[188] - rb[188];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[19] -= rop;
    wdot[27] += rop;
    rop = rf[189] - rb[189];
    wdot[1] += rop;
    wdot[15] -= rop;
    wdot[19] += rop;
    rop = rf[190] - rb[190];
    wdot[1] -= rop;
    wdot[20] += rop;
    wdot[21] -= rop;
    rop = rf[191] - rb[191];
    wdot[4] -= rop;
    wdot[11] += rop;
    wdot[12] += rop;
    wdot[20] -= rop;
    rop = rf[192] - rb[192];
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[20] -= rop;
    wdot[24] += rop;
    rop = rf[193];
    wdot[1] += rop;
    wdot[4] -= rop;
    wdot[9] += rop;
    wdot[11] += rop;
    wdot[20] -= rop;
    rop = rf[194] - rb[194];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[20] -= rop;
    wdot[21] += rop;
    rop = rf[195] - rb[195];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[20] -= rop;
    wdot[21] += rop;
    rop = rf[196] - rb[196];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[20] -= rop;
    rop = rf[197] - rb[197];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[20] -= rop;
    wdot[21] += rop;
    rop = rf[198] - rb[198];
    wdot[19] += rop;
    wdot[20] = wdot[20] - rop - rop;
    wdot[21] += rop;
    rop = rf[199] - rb[199];
    wdot[1] -= rop;
    wdot[21] += rop;
    rop = rf[200] - rb[200];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[26] += rop;
    rop = rf[201] - rb[201];
    wdot[4] -= rop;
    wdot[9] += rop;
    wdot[12] += rop;
    rop = rf[202] - rb[202];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[21] += rop;
    rop = rf[203] - rb[203];
    wdot[21] -= rop;
    rop = rf[204] - rb[204];
    wdot[3] -= rop;
    wdot[9] += rop;
    wdot[16] += rop;
    wdot[21] -= rop;
    rop = rf[205] - rb[205];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[21] -= rop;
    wdot[26] += rop;
    rop = rf[206] - rb[206];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[21] -= rop;
    rop = rf[207] - rb[207];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[21] -= rop;
    wdot[25] += rop;
    rop = rf[208] - rb[208];
    wdot[5] -= rop;
    wdot[9] += rop;
    wdot[15] += rop;
    wdot[21] -= rop;
    rop = rf[209] - rb[209];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[20] += rop;
    wdot[21] -= rop;
    rop = rf[210] - rb[210];
    wdot[1] += rop;
    wdot[16] -= rop;
    wdot[21] -= rop;
    wdot[35] += rop;
    rop = rf[211] - rb[211];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[35] += rop;
    rop = rf[212] - rb[212];
    wdot[9] += rop;
    wdot[21] -= rop;
    wdot[26] -= rop;
    wdot[35] += rop;
    rop = rf[213] - rb[213];
    wdot[21] += rop;
    rop = rf[214] - rb[214];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[25] += rop;
    rop = rf[215] - rb[215];
    wdot[4] -= rop;
    wdot[12] = wdot[12] + rop + rop;
    rop = rf[216] - rb[216];
    wdot[12] -= rop;
    wdot[20] -= rop;
    wdot[28] += rop;
    rop = rf[217] - rb[217];
    wdot[15] += rop;
    wdot[18] += rop;
    wdot[29] -= rop;
    rop = rf[218] - rb[218];
    wdot[1] -= rop;
    wdot[29] += rop;
    rop = rf[219] - rb[219];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[29] -= rop;
    rop = rf[220] - rb[220];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[29] -= rop;
    rop = rf[221] - rb[221];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[29] -= rop;
    rop = rf[222] - rb[222];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[29] -= rop;
    rop = rf[223] - rb[223];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[29] -= rop;
    rop = rf[224] - rb[224];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[29] -= rop;
    rop = rf[225] - rb[225];
    wdot[19] += rop;
    wdot[20] -= rop;
    wdot[29] -= rop;
    rop = rf[226] - rb[226];
    wdot[17] += rop;
    wdot[18] -= rop;
    wdot[29] -= rop;
    rop = rf[227] - rb[227];
    wdot[29] -= rop;
    wdot[30] += rop;
    wdot[31] -= rop;
    rop = rf[228] - rb[228];
    wdot[13] += rop;
    wdot[29] -= rop;
    rop = rf[229] - rb[229];
    wdot[15] -= rop;
    wdot[19] -= rop;
    rop = rf[230] - rb[230];
    wdot[1] -= rop;
    wdot[30] -= rop;
    rop = rf[231] - rb[231];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[30] += rop;
    rop = rf[232] - rb[232];
    wdot[15] -= rop;
    wdot[20] -= rop;
    wdot[30] += rop;
    rop = rf[233] - rb[233];
    wdot[1] -= rop;
    wdot[30] += rop;
    wdot[31] -= rop;
    rop = rf[234] - rb[234];
    wdot[1] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[235] - rb[235];
    wdot[1] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[236] - rb[236];
    wdot[3] -= rop;
    wdot[12] += rop;
    wdot[18] += rop;
    wdot[30] -= rop;
    rop = rf[237];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[15] += rop;
    wdot[25] += rop;
    wdot[30] -= rop;
    rop = rf[238] - rb[238];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[239] - rb[239];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[240] - rb[240];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[241] - rb[241];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[242] - rb[242];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[243] - rb[243];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[244] - rb[244];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[245] - rb[245];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[246] - rb[246];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[247] - rb[247];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[248] - rb[248];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[249] - rb[249];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[250] - rb[250];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[19] += rop;
    wdot[30] -= rop;
    rop = rf[251] - rb[251];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[252] - rb[252];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[253] - rb[253];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[254] - rb[254];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[255] - rb[255];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[30] -= rop;
    wdot[32] += rop;
    rop = rf[256] - rb[256];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[30] -= rop;
    wdot[33] += rop;
    rop = rf[257] - rb[257];
    wdot[17] += rop;
    wdot[18] -= rop;
    wdot[30] -= rop;
    wdot[31] += rop;
    rop = rf[258] - rb[258];
    wdot[15] -= rop;
    wdot[21] -= rop;
    wdot[31] += rop;
    rop = rf[259] - rb[259];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[28] += rop;
    wdot[31] -= rop;
    rop = rf[260];
    wdot[1] = wdot[1] + rop + rop;
    wdot[5] -= rop;
    wdot[28] += rop;
    wdot[31] -= rop;
    rop = rf[261] - rb[261];
    wdot[18] -= rop;
    wdot[19] += rop;
    wdot[30] += rop;
    wdot[31] -= rop;
    rop = rf[262] - rb[262];
    wdot[4] -= rop;
    wdot[11] += rop;
    wdot[31] -= rop;
    rop = rf[263] - rb[263];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[28] += rop;
    wdot[31] -= rop;
    rop = rf[264] - rb[264];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[30] += rop;
    wdot[31] -= rop;
    rop = rf[265] - rb[265];
    wdot[1] += rop;
    wdot[15] -= rop;
    wdot[20] -= rop;
    wdot[31] += rop;
    rop = rf[266] - rb[266];
    wdot[31] -= rop;
    wdot[33] += rop;
    rop = rf[267] - rb[267];
    wdot[31] -= rop;
    wdot[32] += rop;
    rop = rf[268] - rb[268];
    wdot[15] -= rop;
    wdot[21] -= rop;
    wdot[32] += rop;
    rop = rf[269] - rb[269];
    wdot[4] -= rop;
    wdot[12] += rop;
    wdot[22] += rop;
    wdot[32] -= rop;
    rop = rf[270] - rb[270];
    wdot[15] -= rop;
    wdot[21] -= rop;
    wdot[33] += rop;
    rop = rf[271] - rb[271];
    wdot[1] -= rop;
    wdot[33] += rop;
    wdot[34] -= rop;
    rop = rf[272] - rb[272];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[32] -= rop;
    wdot[34] += rop;
    rop = rf[273] - rb[273];
    wdot[3] -= rop;
    wdot[12] += rop;
    wdot[19] += rop;
    wdot[32] -= rop;
    rop = rf[274];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[12] += rop;
    wdot[19] += rop;
    wdot[32] -= rop;
    rop = rf[275];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[12] += rop;
    wdot[19] += rop;
    wdot[32] -= rop;
    rop = rf[276] - rb[276];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[30] += rop;
    wdot[32] -= rop;
    rop = rf[277] - rb[277];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[32] -= rop;
    wdot[34] += rop;
    rop = rf[278] - rb[278];
    wdot[32] += rop;
    wdot[33] -= rop;
    rop = rf[279] - rb[279];
    wdot[4] -= rop;
    wdot[11] += rop;
    wdot[33] -= rop;
    rop = rf[280] - rb[280];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[33] -= rop;
    wdot[34] += rop;
    rop = rf[281] - rb[281];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[33] -= rop;
    wdot[34] += rop;
    rop = rf[282] - rb[282];
    wdot[3] -= rop;
    wdot[15] += rop;
    wdot[25] += rop;
    wdot[33] -= rop;
    rop = rf[283];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[15] += rop;
    wdot[25] += rop;
    wdot[33] -= rop;
    rop = rf[284];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[15] += rop;
    wdot[25] += rop;
    wdot[33] -= rop;
    rop = rf[285] - rb[285];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[30] += rop;
    wdot[33] -= rop;
    rop = rf[286] - rb[286];
    wdot[1] -= rop;
    wdot[32] += rop;
    wdot[34] -= rop;
    rop = rf[287] - rb[287];
    wdot[1] -= rop;
    wdot[31] += rop;
    wdot[34] -= rop;
    rop = rf[288] - rb[288];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[289] - rb[289];
    wdot[3] -= rop;
    wdot[15] += rop;
    wdot[26] += rop;
    wdot[34] -= rop;
    rop = rf[290] - rb[290];
    wdot[3] -= rop;
    wdot[9] += rop;
    wdot[19] += rop;
    wdot[34] -= rop;
    rop = rf[291] - rb[291];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[292] - rb[292];
    wdot[21] += rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[293] - rb[293];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[294] - rb[294];
    wdot[1] += rop;
    wdot[15] -= rop;
    wdot[21] -= rop;
    wdot[34] += rop;
    rop = rf[295] - rb[295];
    wdot[1] -= rop;
    wdot[34] += rop;
    wdot[35] -= rop;
    rop = rf[296] - rb[296];
    wdot[15] -= rop;
    wdot[34] += rop;
    rop = rf[297] - rb[297];
    wdot[4] += rop;
    wdot[7] -= rop;
    wdot[34] += rop;
    wdot[35] -= rop;
    rop = rf[298];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[9] += rop;
    wdot[19] += rop;
    wdot[34] -= rop;
    rop = rf[299] - rb[299];
    wdot[5] -= rop;
    wdot[15] += rop;
    wdot[25] += rop;
    wdot[34] -= rop;
    rop = rf[300] - rb[300];
    wdot[3] -= rop;
    wdot[12] += rop;
    wdot[20] += rop;
    wdot[34] -= rop;
    rop = rf[301] - rb[301];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[302] - rb[302];
    wdot[19] += rop;
    wdot[20] -= rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[303] - rb[303];
    wdot[30] += rop;
    wdot[31] -= rop;
    wdot[34] -= rop;
    wdot[35] += rop;
    rop = rf[304] - rb[304];
    wdot[3] -= rop;
    wdot[11] += rop;
    wdot[35] -= rop;
    rop = rf[305] - rb[305];
    wdot[4] -= rop;
    wdot[12] += rop;
    wdot[25] += rop;
    wdot[35] -= rop;
    rop = rf[306];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[9] += rop;
    wdot[20] += rop;
    wdot[35] -= rop;
    rop = rf[307] - rb[307];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[34] += rop;
    wdot[35] -= rop;
    rop = rf[308] - rb[308];
    wdot[15] += rop;
    wdot[18] -= rop;
    wdot[35] += rop;
    rop = rf[309] - rb[309];
    wdot[1] += rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[310] - rb[310];
    wdot[2] += rop;
    wdot[36] -= rop;
    wdot[37] += rop;
    rop = rf[311] - rb[311];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[312] - rb[312];
    wdot[1] -= rop;
    wdot[19] += rop;
    wdot[20] += rop;
    wdot[36] -= rop;
    rop = rf[313] - rb[313];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[34] += rop;
    wdot[36] -= rop;
    rop = rf[314] - rb[314];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[315] - rb[315];
    wdot[3] -= rop;
    wdot[21] += rop;
    wdot[27] += rop;
    wdot[36] -= rop;
    rop = rf[316] - rb[316];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[317] - rb[317];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[318] - rb[318];
    wdot[19] += rop;
    wdot[20] -= rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[319] - rb[319];
    wdot[30] += rop;
    wdot[31] -= rop;
    wdot[36] -= rop;
    wdot[41] += rop;
    rop = rf[320] - rb[320];
    wdot[1] += rop;
    wdot[20] -= rop;
    wdot[21] -= rop;
    wdot[37] += rop;
    rop = rf[321] - rb[321];
    wdot[20] -= rop;
    wdot[21] -= rop;
    wdot[41] += rop;
    rop = rf[322] - rb[322];
    wdot[20] = wdot[20] - rop - rop;
    wdot[36] += rop;
    rop = rf[323] - rb[323];
    wdot[1] += rop;
    wdot[20] = wdot[20] - rop - rop;
    wdot[41] += rop;
    rop = rf[324] - rb[324];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[37] += rop;
    wdot[41] -= rop;
    rop = rf[325] - rb[325];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[35] += rop;
    wdot[41] -= rop;
    rop = rf[326] - rb[326];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[37] += rop;
    wdot[41] -= rop;
    rop = rf[327] - rb[327];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[36] += rop;
    wdot[41] -= rop;
    rop = rf[328] - rb[328];
    wdot[4] += rop;
    wdot[7] -= rop;
    wdot[36] += rop;
    wdot[41] -= rop;
    rop = rf[329];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[20] += rop;
    wdot[25] += rop;
    wdot[41] -= rop;
    rop = rf[330] - rb[330];
    wdot[7] += rop;
    wdot[8] -= rop;
    wdot[36] += rop;
    wdot[41] -= rop;
    rop = rf[331] - rb[331];
    wdot[4] -= rop;
    wdot[24] += rop;
    wdot[25] += rop;
    wdot[41] -= rop;
    rop = rf[332] - rb[332];
    wdot[41] += rop;
    wdot[42] -= rop;
    rop = rf[333] - rb[333];
    wdot[41] += rop;
    wdot[42] -= rop;
    rop = rf[334];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[21] += rop;
    wdot[42] -= rop;
    rop = rf[335] - rb[335];
    wdot[4] -= rop;
    wdot[25] += rop;
    wdot[42] -= rop;
    rop = rf[336] - rb[336];
    wdot[36] += rop;
    wdot[43] -= rop;
    rop = rf[337] - rb[337];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[42] += rop;
    wdot[43] -= rop;
    rop = rf[338] - rb[338];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[34] += rop;
    wdot[43] -= rop;
    rop = rf[339] - rb[339];
    wdot[1] += rop;
    wdot[42] += rop;
    wdot[43] -= rop;
    rop = rf[340] - rb[340];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[42] += rop;
    wdot[43] -= rop;
    rop = rf[341] - rb[341];
    wdot[1] -= rop;
    wdot[37] -= rop;
    wdot[41] += rop;
    rop = rf[342] - rb[342];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[37] -= rop;
    wdot[39] += rop;
    rop = rf[343] - rb[343];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[37] -= rop;
    wdot[38] += rop;
    rop = rf[344] - rb[344];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[37] -= rop;
    wdot[39] += rop;
    rop = rf[345] - rb[345];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[37] -= rop;
    wdot[38] += rop;
    rop = rf[346] - rb[346];
    wdot[3] -= rop;
    wdot[12] += rop;
    wdot[35] += rop;
    wdot[37] -= rop;
    rop = rf[347] - rb[347];
    wdot[9] += rop;
    wdot[26] -= rop;
    wdot[35] -= rop;
    wdot[37] += rop;
    rop = rf[348] - rb[348];
    wdot[1] += rop;
    wdot[16] -= rop;
    wdot[35] -= rop;
    wdot[37] += rop;
    rop = rf[349] - rb[349];
    wdot[21] -= rop;
    wdot[39] += rop;
    rop = rf[350] - rb[350];
    wdot[38] += rop;
    wdot[39] -= rop;
    rop = rf[351] - rb[351];
    wdot[38] += rop;
    wdot[39] -= rop;
    rop = rf[352] - rb[352];
    wdot[1] -= rop;
    wdot[21] += rop;
    wdot[39] -= rop;
    rop = rf[353] - rb[353];
    wdot[1] -= rop;
    wdot[37] += rop;
    wdot[39] -= rop;
    rop = rf[354] - rb[354];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[39] -= rop;
    wdot[40] += rop;
    rop = rf[355] - rb[355];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[39] -= rop;
    wdot[40] += rop;
    rop = rf[356] - rb[356];
    wdot[21] -= rop;
    wdot[38] += rop;
    rop = rf[357] - rb[357];
    wdot[1] -= rop;
    wdot[21] += rop;
    wdot[38] -= rop;
    rop = rf[358] - rb[358];
    wdot[1] -= rop;
    wdot[37] += rop;
    wdot[38] -= rop;
    rop = rf[359] - rb[359];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[38] -= rop;
    wdot[40] += rop;
    rop = rf[360] - rb[360];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[38] -= rop;
    wdot[40] += rop;
    rop = rf[361] - rb[361];
    wdot[4] -= rop;
    wdot[25] += rop;
    wdot[26] += rop;
    wdot[38] -= rop;
    rop = rf[362] - rb[362];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[40] += rop;
    rop = rf[363] - rb[363];
    wdot[1] -= rop;
    wdot[39] += rop;
    wdot[40] -= rop;
    rop = rf[364] - rb[364];
    wdot[1] -= rop;
    wdot[38] += rop;
    wdot[40] -= rop;
    rop = rf[365] - rb[365];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[40] -= rop;
    wdot[44] += rop;
    rop = rf[366] - rb[366];
    wdot[1] -= rop;
    wdot[21] += rop;
    wdot[26] += rop;
    wdot[44] -= rop;
    rop = rf[367] - rb[367];
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[26] += rop;
    wdot[44] -= rop;
    rop = rf[368] - rb[368];
    wdot[21] -= rop;
    wdot[37] += rop;
    rop = rf[369] - rb[369];
    wdot[19] -= rop;
    wdot[36] += rop;
    rop = rf[370];
    wdot[1] = wdot[1] + rop + rop;
    wdot[20] -= rop;
    wdot[31] -= rop;
    wdot[46] += rop;
    rop = rf[371];
    wdot[1] = wdot[1] + rop + rop;
    wdot[31] -= rop;
    wdot[35] -= rop;
    wdot[45] += rop;
    rop = rf[372] - rb[372];
    wdot[1] += rop;
    wdot[35] = wdot[35] - rop - rop;
    wdot[47] += rop;
    rop = rf[373] - rb[373];
    wdot[35] = wdot[35] - rop - rop;
    wdot[45] += rop;
    rop = rf[374];
    wdot[1] += rop;
    wdot[2] += rop;
    wdot[20] -= rop;
    wdot[36] -= rop;
    wdot[45] += rop;
    rop = rf[375] - rb[375];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[42] -= rop;
    wdot[45] += rop;
    rop = rf[376] - rb[376];
    wdot[15] += rop;
    wdot[19] -= rop;
    wdot[42] -= rop;
    wdot[46] += rop;
    rop = rf[377] - rb[377];
    wdot[21] -= rop;
    wdot[39] -= rop;
    wdot[47] += rop;
    rop = rf[378] - rb[378];
    wdot[15] -= rop;
    wdot[38] -= rop;
    wdot[46] += rop;
    rop = rf[379] - rb[379];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[45] -= rop;
    wdot[47] += rop;
    rop = rf[380] - rb[380];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[45] -= rop;
    wdot[47] += rop;
    rop = rf[381] - rb[381];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[45] -= rop;
    wdot[47] += rop;
    rop = rf[382];
    wdot[1] += rop;
    wdot[21] += rop;
    wdot[40] += rop;
    wdot[47] -= rop;
    rop = rf[383] - rb[383];
    wdot[11] -= rop;
    wdot[12] += rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[384] - rb[384];
    wdot[9] += rop;
    wdot[12] -= rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[385] - rb[385];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[47] -= rop;
    wdot[48] += rop;
    rop = rf[386] - rb[386];
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[47] -= rop;
    wdot[48] += rop;
    rop = rf[387] - rb[387];
    wdot[3] -= rop;
    wdot[10] += rop;
    wdot[48] -= rop;
    wdot[49] += rop;
    rop = rf[388] - rb[388];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[389] - rb[389];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[390] - rb[390];
    wdot[1] -= rop;
    wdot[21] += rop;
    wdot[31] += rop;
    wdot[46] -= rop;
    rop = rf[391] - rb[391];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[392] - rb[392];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[49] -= rop;
    wdot[50] += rop;
    rop = rf[393] - rb[393];
    wdot[5] -= rop;
    wdot[9] += rop;
    wdot[36] += rop;
    wdot[49] -= rop;
    rop = rf[394] - rb[394];
    wdot[3] -= rop;
    wdot[10] += rop;
    wdot[37] += rop;
    wdot[50] -= rop;
    rop = rf[395] - rb[395];
    wdot[11] += rop;
    wdot[12] -= rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[396] - rb[396];
    wdot[15] += rop;
    wdot[20] -= rop;
    wdot[45] += rop;
    wdot[46] -= rop;
    rop = rf[397] - rb[397];
    wdot[36] += rop;
    wdot[41] -= rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[398] - rb[398];
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[49] -= rop;
    rop = rf[399] - rb[399];
    wdot[3] -= rop;
    wdot[49] -= rop;
    rop = rf[400] - rb[400];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[49] -= rop;
    rop = rf[401] - rb[401];
    wdot[1] += rop;
    wdot[35] -= rop;
    wdot[36] -= rop;
    wdot[52] += rop;
    rop = rf[402] - rb[402];
    wdot[1] -= rop;
    wdot[51] -= rop;
    wdot[52] += rop;
    rop = rf[403] - rb[403];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[51] -= rop;
    wdot[54] += rop;
    rop = rf[404] - rb[404];
    wdot[1] += rop;
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[51] -= rop;
    wdot[54] += rop;
    rop = rf[405] - rb[405];
    wdot[1] -= rop;
    wdot[12] += rop;
    wdot[45] += rop;
    wdot[54] -= rop;
    rop = rf[406] - rb[406];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[45] -= rop;
    wdot[48] += rop;
    rop = rf[407] - rb[407];
    wdot[1] -= rop;
    wdot[9] += rop;
    wdot[46] += rop;
    wdot[48] -= rop;
    rop = rf[408] - rb[408];
    wdot[9] += rop;
    wdot[48] -= rop;
    wdot[49] += rop;
    rop = rf[409] - rb[409];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[410] - rb[410];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[411] - rb[411];
    wdot[21] -= rop;
    wdot[35] -= rop;
    wdot[49] += rop;
    rop = rf[412] - rb[412];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[49] -= rop;
    wdot[50] += rop;
    rop = rf[413] - rb[413];
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[49] -= rop;
    rop = rf[414] - rb[414];
    wdot[1] += rop;
    wdot[50] += rop;
    rop = rf[415] - rb[415];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[50] += rop;
    rop = rf[416] - rb[416];
    wdot[9] += rop;
    wdot[21] = wdot[21] + rop + rop;
    wdot[50] -= rop;
    rop = rf[417] - rb[417];
    wdot[35] -= rop;
    wdot[37] -= rop;
    wdot[51] += rop;
    rop = rf[418] - rb[418];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[45] += rop;
    wdot[52] -= rop;
    rop = rf[419] - rb[419];
    wdot[4] += rop;
    wdot[7] -= rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[420] - rb[420];
    wdot[1] -= rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[421] - rb[421];
    wdot[12] -= rop;
    wdot[47] -= rop;
    wdot[54] += rop;
    rop = rf[422] - rb[422];
    wdot[15] -= rop;
    wdot[47] -= rop;
    wdot[52] += rop;
    rop = rf[423] - rb[423];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[424] - rb[424];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[425] - rb[425];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[426] - rb[426];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[427] - rb[427];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[428] - rb[428];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[429] - rb[429];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[47] += rop;
    wdot[51] -= rop;
    rop = rf[430] - rb[430];
    wdot[45] += rop;
    wdot[47] -= rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[431] - rb[431];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[432] - rb[432];
    wdot[1] -= rop;
    wdot[46] += rop;
    wdot[49] -= rop;
    rop = rf[433] - rb[433];
    wdot[1] += rop;
    wdot[45] -= rop;
    wdot[53] += rop;
    rop = rf[434] - rb[434];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[47] -= rop;
    wdot[53] += rop;
    rop = rf[435] - rb[435];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[53] -= rop;
    wdot[55] += rop;
    rop = rf[436] - rb[436];
    wdot[1] -= rop;
    wdot[53] += rop;
    wdot[55] -= rop;
    rop = rf[437] - rb[437];
    wdot[1] += rop;
    wdot[20] -= rop;
    wdot[45] -= rop;
    wdot[56] += rop;
    rop = rf[438] - rb[438];
    wdot[20] -= rop;
    wdot[47] -= rop;
    wdot[56] += rop;
    rop = rf[439];
    wdot[21] -= rop;
    wdot[55] -= rop;
    rop = rf[440];
    wdot[21] += rop;
    wdot[55] += rop;
    rop = rf[441];
    rop = rf[442];
    rop = rf[443];
    wdot[57] += rop;
    rop = rf[444];
    wdot[57] -= rop;
    rop = rf[445] - rb[445];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[57] += rop;
    wdot[58] -= rop;
    rop = rf[446] - rb[446];
    wdot[1] -= rop;
    wdot[57] -= rop;
    wdot[58] += rop;
    rop = rf[447] - rb[447];
    wdot[21] -= rop;
    wdot[57] -= rop;
    rop = rf[448];
    wdot[1] += rop;
    wdot[37] -= rop;
    wdot[47] -= rop;
    wdot[58] += rop;
    rop = rf[449] - rb[449];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[59] -= rop;
    wdot[61] += rop;
    rop = rf[450] - rb[450];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[59] -= rop;
    wdot[61] += rop;
    rop = rf[451] - rb[451];
    wdot[1] -= rop;
    wdot[59] += rop;
    wdot[61] -= rop;
    rop = rf[452] - rb[452];
    wdot[1] += rop;
    wdot[59] += rop;
    rop = rf[453];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[60] -= rop;
    wdot[62] += rop;
    rop = rf[454] - rb[454];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[62] -= rop;
    wdot[63] += rop;
    rop = rf[455] - rb[455];
    wdot[1] -= rop;
    wdot[62] += rop;
    wdot[63] -= rop;
    rop = rf[456] - rb[456];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[62] -= rop;
    wdot[64] += rop;
    rop = rf[457] - rb[457];
    wdot[1] -= rop;
    wdot[62] += rop;
    wdot[64] -= rop;
    rop = rf[458] - rb[458];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[62] -= rop;
    wdot[65] += rop;
    rop = rf[459] - rb[459];
    wdot[1] -= rop;
    wdot[62] += rop;
    wdot[65] -= rop;
    rop = rf[460] - rb[460];
    wdot[3] -= rop;
    wdot[9] += rop;
    wdot[12] += rop;
    wdot[21] = wdot[21] + rop + rop;
    wdot[48] -= rop;
    rop = rf[461];
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[47] += rop;
    wdot[53] -= rop;
    rop = rf[462];
    wdot[5] -= rop;
    wdot[21] += rop;
    wdot[48] += rop;
    wdot[53] -= rop;
    rop = rf[463];
    wdot[5] -= rop;
    wdot[19] += rop;
    wdot[48] += rop;
    wdot[56] -= rop;
    rop = rf[464];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[53] += rop;
    wdot[58] -= rop;
    rop = rf[465];
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[60] += rop;
    wdot[62] -= rop;
    rop = rf[466];
    wdot[3] -= rop;
    wdot[26] += rop;
    wdot[47] += rop;
    wdot[53] -= rop;
    rop = rf[467];
    wdot[3] -= rop;
    wdot[9] += rop;
    wdot[15] += rop;
    wdot[47] += rop;
    wdot[56] -= rop;
    rop = rf[468];
    wdot[3] -= rop;
    wdot[48] += rop;
    wdot[53] -= rop;
    rop = rf[469];
    wdot[3] -= rop;
    wdot[20] += rop;
    wdot[48] += rop;
    wdot[56] -= rop;
    rop = rf[470];
    wdot[3] -= rop;
    wdot[25] += rop;
    wdot[53] += rop;
    wdot[58] -= rop;
    rop = rf[471];
    wdot[3] -= rop;
    wdot[26] += rop;
    wdot[60] += rop;
    wdot[62] -= rop;
    rop = rf[472];
    wdot[1] += rop;
    wdot[20] -= rop;
    wdot[53] -= rop;
    wdot[58] += rop;
    rop = rf[473] - rb[473];
    wdot[1] += rop;
    wdot[53] -= rop;
    wdot[55] += rop;
    rop = rf[474] - rb[474];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[56] += rop;
    rop = rf[475] - rb[475];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[56] += rop;
    rop = rf[476] - rb[476];
    wdot[1] -= rop;
    wdot[66] += rop;
    rop = rf[477] - rb[477];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[56] += rop;
    rop = rf[478] - rb[478];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[56] += rop;
    rop = rf[479] - rb[479];
    wdot[1] += rop;
    wdot[56] += rop;
    rop = rf[480] - rb[480];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[66] -= rop;
    rop = rf[481] - rb[481];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[66] -= rop;
    rop = rf[482] - rb[482];
    wdot[1] -= rop;
    wdot[18] += rop;
    wdot[45] += rop;
    wdot[66] -= rop;
    rop = rf[483] - rb[483];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[66] -= rop;
    rop = rf[484] - rb[484];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[66] -= rop;
    rop = rf[485] - rb[485];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[66] -= rop;
    rop = rf[486] - rb[486];
    wdot[18] += rop;
    wdot[47] += rop;
    wdot[66] -= rop;
    rop = rf[487];
    wdot[1] += rop;
    wdot[19] -= rop;
    wdot[55] -= rop;
    wdot[58] += rop;
    rop = rf[488] - rb[488];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[45] -= rop;
    wdot[47] += rop;
    rop = rf[489] - rb[489];
    wdot[19] -= rop;
    wdot[47] -= rop;
    rop = rf[490] - rb[490];
    wdot[19] -= rop;
    wdot[20] += rop;
    wdot[45] += rop;
    wdot[47] -= rop;
    rop = rf[491] - rb[491];
    wdot[1] += rop;
    wdot[11] -= rop;
    wdot[47] -= rop;
    wdot[54] += rop;
    rop = rf[492];
    wdot[9] += rop;
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[493];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[9] += rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[494];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[9] += rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[495];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[9] += rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[496];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[9] += rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[497];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[9] += rop;
    wdot[47] += rop;
    wdot[54] -= rop;
    rop = rf[498] - rb[498];
    wdot[3] -= rop;
    wdot[47] -= rop;
    wdot[48] += rop;
    rop = rf[499] - rb[499];
    wdot[2] += rop;
    wdot[20] -= rop;
    wdot[58] -= rop;
    rop = rf[500] - rb[500];
    wdot[1] -= rop;
    wdot[15] += rop;
    wdot[58] += rop;
    wdot[68] -= rop;
    rop = rf[501] - rb[501];
    wdot[15] += rop;
    wdot[57] += rop;
    wdot[68] -= rop;
    rop = rf[502] - rb[502];
    wdot[1] += rop;
    wdot[59] -= rop;
    wdot[61] += rop;
    rop = rf[503] - rb[503];
    wdot[1] += rop;
    wdot[20] -= rop;
    wdot[57] -= rop;
    rop = rf[504] - rb[504];
    wdot[2] += rop;
    wdot[19] -= rop;
    wdot[57] -= rop;
    rop = rf[505];
    wdot[21] += rop;
    wdot[60] -= rop;
    wdot[61] += rop;
    rop = rf[506];
    wdot[21] -= rop;
    wdot[60] += rop;
    wdot[61] -= rop;
    rop = rf[507] - rb[507];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[64] -= rop;
    wdot[71] += rop;
    rop = rf[508] - rb[508];
    wdot[2] += rop;
    wdot[35] -= rop;
    wdot[59] += rop;
    wdot[67] -= rop;
    rop = rf[509] - rb[509];
    wdot[6] += rop;
    wdot[7] -= rop;
    wdot[67] -= rop;
    wdot[69] += rop;
    rop = rf[510];
    wdot[1] += rop;
    wdot[5] += rop;
    wdot[7] -= rop;
    wdot[67] -= rop;
    wdot[69] += rop;
    rop = rf[511] - rb[511];
    wdot[1] += rop;
    wdot[3] -= rop;
    wdot[67] -= rop;
    wdot[69] += rop;
    rop = rf[512] - rb[512];
    wdot[4] -= rop;
    wdot[5] += rop;
    wdot[67] -= rop;
    wdot[69] += rop;
    rop = rf[513];
    wdot[1] += rop;
    wdot[3] += rop;
    wdot[4] -= rop;
    wdot[67] -= rop;
    wdot[69] += rop;
    rop = rf[514] - rb[514];
    wdot[11] += rop;
    wdot[12] -= rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[515] - rb[515];
    wdot[7] -= rop;
    wdot[8] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[516] - rb[516];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[517];
    wdot[1] = wdot[1] + rop + rop;
    wdot[3] -= rop;
    wdot[69] += rop;
    wdot[70] -= rop;
    rop = rf[518] - rb[518];
    wdot[4] -= rop;
    wdot[7] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[519] - rb[519];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[520] - rb[520];
    wdot[1] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[521] - rb[521];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[51] -= rop;
    wdot[70] += rop;
    rop = rf[522] - rb[522];
    wdot[21] += rop;
    wdot[49] += rop;
    wdot[51] -= rop;
    rop = rf[523] - rb[523];
    wdot[3] -= rop;
    wdot[11] += rop;
    wdot[51] += rop;
    rop = rf[524] - rb[524];
    wdot[15] += rop;
    wdot[51] += rop;
    wdot[66] -= rop;
    rop = rf[525];
    wdot[1] = wdot[1] + rop + rop;
    wdot[49] = wdot[49] - rop - rop;
    wdot[58] += rop;
    rop = rf[526] - rb[526];
    wdot[5] -= rop;
    wdot[12] += rop;
    wdot[36] += rop;
    wdot[46] -= rop;
    rop = rf[527] - rb[527];
    wdot[6] += rop;
    wdot[7] -= rop;
    wdot[49] -= rop;
    wdot[50] += rop;
    rop = rf[528] - rb[528];
    wdot[46] += rop;
    wdot[49] -= rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[529] - rb[529];
    wdot[45] += rop;
    wdot[46] -= rop;
    wdot[47] -= rop;
    wdot[49] += rop;
    rop = rf[530] - rb[530];
    wdot[19] += rop;
    wdot[20] -= rop;
    wdot[46] -= rop;
    wdot[49] += rop;
    rop = rf[531] - rb[531];
    wdot[21] += rop;
    wdot[35] += rop;
    wdot[49] -= rop;
    rop = rf[532] - rb[532];
    wdot[15] += rop;
    wdot[20] -= rop;
    wdot[36] -= rop;
    wdot[46] += rop;
    rop = rf[533] - rb[533];
    wdot[1] += rop;
    wdot[16] -= rop;
    wdot[51] -= rop;
    wdot[56] += rop;
    rop = rf[534] - rb[534];
    wdot[15] -= rop;
    wdot[51] -= rop;
    wdot[66] += rop;
    rop = rf[535] - rb[535];
    wdot[37] -= rop;
    wdot[47] -= rop;
    wdot[56] += rop;
    rop = rf[536] - rb[536];
    wdot[20] += rop;
    wdot[36] -= rop;
    wdot[47] -= rop;
    wdot[56] += rop;
    rop = rf[537] - rb[537];
    wdot[37] = wdot[37] - rop - rop;
    wdot[56] += rop;
    rop = rf[538] - rb[538];
    wdot[2] += rop;
    wdot[56] += rop;
    wdot[66] -= rop;
    rop = rf[539] - rb[539];
    wdot[2] += rop;
    wdot[56] += rop;
    wdot[66] -= rop;
    rop = rf[540] - rb[540];
    wdot[1] += rop;
    wdot[17] -= rop;
    wdot[47] -= rop;
    wdot[66] += rop;
    rop = rf[541] - rb[541];
    wdot[39] -= rop;
    wdot[47] -= rop;
    wdot[58] += rop;
    rop = rf[542] - rb[542];
    wdot[1] += rop;
    wdot[39] -= rop;
    wdot[47] -= rop;
    wdot[57] += rop;
    rop = rf[543] - rb[543];
    wdot[35] -= rop;
    wdot[47] -= rop;
    wdot[70] += rop;
    rop = rf[544] - rb[544];
    wdot[1] += rop;
    wdot[35] -= rop;
    wdot[45] -= rop;
    wdot[70] += rop;
    rop = rf[545] - rb[545];
    wdot[2] += rop;
    wdot[16] -= rop;
    wdot[52] -= rop;
    wdot[56] += rop;
    rop = rf[546] - rb[546];
    wdot[1] -= rop;
    wdot[14] += rop;
    wdot[47] += rop;
    wdot[52] -= rop;
    rop = rf[547] - rb[547];
    wdot[3] -= rop;
    wdot[5] += rop;
    wdot[51] += rop;
    wdot[52] -= rop;
    rop = rf[548];
    wdot[2] = wdot[2] + rop + rop;
    wdot[51] -= rop;
    wdot[62] += rop;
    wdot[67] -= rop;
    rop = rf[549];
    wdot[2] += rop;
    wdot[21] += rop;
    wdot[62] += rop;
    wdot[67] = wdot[67] - rop - rop;
    rop = rf[550];
    wdot[1] += rop;
    wdot[53] -= rop;
    wdot[55] -= rop;
    wdot[62] += rop;
    rop = rf[551];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[57] += rop;
    wdot[58] -= rop;
    rop = rf[552];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[57] -= rop;
    wdot[58] += rop;
    rop = rf[553];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[62] -= rop;
    wdot[63] += rop;
    rop = rf[554];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[62] += rop;
    wdot[63] -= rop;
    rop = rf[555];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[62] -= rop;
    wdot[64] += rop;
    rop = rf[556];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[62] += rop;
    wdot[64] -= rop;
    rop = rf[557];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[62] -= rop;
    wdot[65] += rop;
    rop = rf[558];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[62] += rop;
    wdot[65] -= rop;
    rop = rf[559];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[57] += rop;
    wdot[58] -= rop;
    rop = rf[560];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[57] -= rop;
    wdot[58] += rop;
    rop = rf[561];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[62] -= rop;
    wdot[63] += rop;
    rop = rf[562];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[62] += rop;
    wdot[63] -= rop;
    rop = rf[563];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[62] -= rop;
    wdot[64] += rop;
    rop = rf[564];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[62] += rop;
    wdot[64] -= rop;
    rop = rf[565];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[62] -= rop;
    wdot[65] += rop;
    rop = rf[566];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[62] += rop;
    wdot[65] -= rop;
    rop = rf[567];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[568];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[67] -= rop;
    wdot[70] += rop;
    rop = rf[569];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[67] += rop;
    wdot[70] -= rop;
    rop = rf[570];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[67] -= rop;
    wdot[70] += rop;
    rop = rf[571];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[53] -= rop;
    wdot[55] += rop;
    rop = rf[572];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[53] += rop;
    wdot[55] -= rop;
    rop = rf[573];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[53] -= rop;
    wdot[55] += rop;
    rop = rf[574];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[53] += rop;
    wdot[55] -= rop;
    rop = rf[575];
    wdot[1] += rop;
    wdot[37] -= rop;
    wdot[63] -= rop;
    wdot[72] += rop;
    rop = rf[576];
    wdot[1] += rop;
    wdot[37] -= rop;
    wdot[64] -= rop;
    wdot[72] += rop;
    rop = rf[577];
    wdot[1] += rop;
    wdot[37] -= rop;
    wdot[65] -= rop;
    wdot[73] += rop;
    rop = rf[578] - rb[578];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[63] -= rop;
    wdot[74] += rop;
    rop = rf[579] - rb[579];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[64] -= rop;
    wdot[75] += rop;
    rop = rf[580] - rb[580];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[65] -= rop;
    wdot[76] += rop;
    rop = rf[581];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[74] -= rop;
    wdot[77] += rop;
    rop = rf[582];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[75] -= rop;
    wdot[78] += rop;
    rop = rf[583];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[74] += rop;
    wdot[77] -= rop;
    rop = rf[584];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[75] += rop;
    wdot[78] -= rop;
    rop = rf[585];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[74] -= rop;
    wdot[77] += rop;
    rop = rf[586];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[75] -= rop;
    wdot[78] += rop;
    rop = rf[587];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[74] += rop;
    wdot[77] -= rop;
    rop = rf[588];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[75] += rop;
    wdot[78] -= rop;
    rop = rf[589] - rb[589];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[74] -= rop;
    wdot[77] += rop;
    rop = rf[590] - rb[590];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[75] -= rop;
    wdot[78] += rop;
    rop = rf[591];
    wdot[21] -= rop;
    wdot[77] -= rop;
    wdot[79] += rop;
    rop = rf[592];
    wdot[21] -= rop;
    wdot[78] -= rop;
    wdot[79] += rop;
    rop = rf[593] - rb[593];
    wdot[1] -= rop;
    wdot[72] += rop;
    wdot[79] -= rop;
    rop = rf[594];
    wdot[1] += rop;
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[72] -= rop;
    wdot[75] += rop;
    rop = rf[595] - rb[595];
    wdot[3] -= rop;
    wdot[25] += rop;
    wdot[72] -= rop;
    wdot[75] += rop;
    rop = rf[596];
    wdot[4] -= rop;
    wdot[9] += rop;
    wdot[12] += rop;
    wdot[75] += rop;
    wdot[79] -= rop;
    rop = rf[597] - rb[597];
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[64] += rop;
    wdot[75] -= rop;
    rop = rf[598] - rb[598];
    wdot[3] -= rop;
    wdot[26] += rop;
    wdot[64] += rop;
    wdot[75] -= rop;
    rop = rf[599];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[76] -= rop;
    wdot[80] += rop;
    rop = rf[600];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[76] += rop;
    wdot[80] -= rop;
    rop = rf[601];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[76] -= rop;
    wdot[80] += rop;
    rop = rf[602];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[76] += rop;
    wdot[80] -= rop;
    rop = rf[603] - rb[603];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[76] -= rop;
    wdot[80] += rop;
    rop = rf[604];
    wdot[21] -= rop;
    wdot[80] -= rop;
    wdot[81] += rop;
    rop = rf[605] - rb[605];
    wdot[1] -= rop;
    wdot[73] += rop;
    wdot[81] -= rop;
    rop = rf[606];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[73] -= rop;
    wdot[81] += rop;
    rop = rf[607];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[73] += rop;
    wdot[81] -= rop;
    rop = rf[608];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[73] -= rop;
    wdot[81] += rop;
    rop = rf[609];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[73] += rop;
    wdot[81] -= rop;
    rop = rf[610] - rb[610];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[73] -= rop;
    wdot[81] += rop;
    rop = rf[611];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[81] -= rop;
    wdot[82] += rop;
    rop = rf[612];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[72] -= rop;
    wdot[79] += rop;
    rop = rf[613];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[72] += rop;
    wdot[79] -= rop;
    rop = rf[614];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[72] -= rop;
    wdot[79] += rop;
    rop = rf[615];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[72] += rop;
    wdot[79] -= rop;
    rop = rf[616] - rb[616];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[72] -= rop;
    wdot[79] += rop;
    rop = rf[617];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[79] -= rop;
    wdot[83] += rop;
    rop = rf[618];
    wdot[5] -= rop;
    wdot[25] += rop;
    wdot[79] += rop;
    wdot[83] -= rop;
    rop = rf[619];
    wdot[3] -= rop;
    wdot[26] += rop;
    wdot[79] += rop;
    wdot[83] -= rop;
    rop = rf[620];
    wdot[1] -= rop;
    wdot[2] += rop;
    wdot[82] -= rop;
    wdot[84] += rop;
    rop = rf[621];
    wdot[1] += rop;
    wdot[2] -= rop;
    wdot[82] += rop;
    wdot[84] -= rop;
    rop = rf[622];
    wdot[14] += rop;
    wdot[15] -= rop;
    wdot[82] -= rop;
    wdot[84] += rop;
    rop = rf[623];
    wdot[14] -= rop;
    wdot[15] += rop;
    wdot[82] += rop;
    wdot[84] -= rop;
    rop = rf[624] - rb[624];
    wdot[5] -= rop;
    wdot[6] += rop;
    wdot[82] -= rop;
    wdot[84] += rop;
    rop = rf[625];
    wdot[1] += rop;
    wdot[21] -= rop;
    wdot[84] -= rop;
    wdot[85] += rop;

    return 0;
} /* rdot_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int qssa_(double *rf, double *rb, double *xq)
{
    static double a1_0__, a2_0__, a1_2__, a2_1__, a3_0__, a4_0__, a5_0__,
                  a2_6__, a2_7__, a6_0__, a5_6__, a6_2__, a6_5__, a7_0__, a7_2__,
                  a8_0__, a9_0__, den, a10_0__, a11_0__, a12_0__, a13_0__, a10_11__,
                  a11_10__;


    /* Parameter adjustments */
    --xq;
    --rb;
    --rf;

    /* Function Body */
    rf[59] = 0.;
    rf[63] = 0.;
    rf[76] = 0.;

    /*     CH2OH */
    den = rf[54] + rf[55] + rf[56] + rf[57] + rf[58] + rf[61] + rf[62] + rb[
        53] + rb[60] + rb[66] + rb[67] + rb[69] + rb[70] + rb[72] + rb[73]
            + rb[74] + rb[75] + rb[85] + rb[161];
    a1_0__ = (rf[53] + rb[54] + rb[55] + rb[56] + rb[57] + rb[58] + rb[59] +
            rf[60] + rb[61] + rb[62] + rb[63] + rb[63] + rf[66] + rf[67] + rf[
            69] + rf[70] + rf[72] + rf[73] + rf[74] + rf[85] + rf[161]) / den;
    a1_2__ = rf[75] / den;
    /*     CH3O */
    den = rf[46] + rf[47] + rf[48] + rf[50] + rf[51] + rf[52] + rf[75] + rf[
        115] + rf[133] + rf[187] + rf[228] + rb[49] + rb[68] + rb[71] +
            rb[86] + rb[87] + rb[90];
    a2_0__ = (rb[46] + rb[47] + rb[48] + rf[49] + rb[50] + rb[51] + rb[52] +
            rb[59] + rf[68] + rf[71] + rb[76] + rb[76] + rf[86] + rf[87] + rf[
            90] + rb[115] + rb[187]) / den;
    a2_1__ = rb[75] / den;
    a2_6__ = rb[133] / den;
    a2_7__ = rb[228] / den;
    /*     CH2(S) */
    den = rf[92] + rf[93] + rf[94] + rf[95] + rf[96] + rf[97] + rf[98] + rf[
        99] + rf[100] + rf[101] + rf[102] + rf[116] + rf[163] + rf[189] +
            rf[211] + rb[65] + rb[83] + rb[166];
    a3_0__ = (rf[65] + rf[83] + rb[92] + rb[93] + rb[94] + rb[95] + rb[96] +
            rb[97] + rb[98] + rb[99] + rb[100] + rb[101] + rb[102] + rb[116]
            + rb[163] + rf[166] + rb[189] + rb[211]) / den;
    /*     C2H */
    den = rf[199] + rf[200] + rf[201] + rf[202] + rf[292] + rf[296] + rf[308]
        + rf[349] + rf[356] + rf[362] + rf[433] + rb[206] + rb[304] + rb[
        468] + rb[535];
    a4_0__ = (rb[199] + rb[200] + rb[201] + rb[202] + rf[206] + rb[292] + rb[
            296] + rf[304] + rb[308] + rb[349] + rb[356] + rb[362] + rb[433]
            + rf[468] + rf[535]) / den;
    /*     CH3CO */
    den = rf[146] + rf[147] + rf[148] + rf[149] + rf[155] + rb[134] + rb[138]
        + rb[140] + rb[141] + rb[142] + rb[143] + rb[144] + rb[262] + rb[
        279] + rb[334] + rb[335];
    a5_0__ = (rf[138] + rf[140] + rf[141] + rf[142] + rf[143] + rf[144] + rb[
            146] + rb[147] + rb[148] + rb[149] + rb[155] + rf[262] + rf[279]
            + rf[334] + rf[335]) / den;
    a5_6__ = rf[134] / den;
    /*     C2H3O1-2 */
    den = rf[134] + rf[135] + rb[129] + rb[130] + rb[131] + rb[132] + rb[133];
    a6_0__ = (rf[129] + rf[130] + rf[131] + rf[132] + rb[135]) / den;
    a6_2__ = rf[133] / den;
    a6_5__ = rb[134] / den;
    /*     NC3H7 */
    den = rf[218] + rf[231] + rb[219] + rb[220] + rb[221] + rb[222] + rb[223]
        + rb[224] + rb[225] + rb[226] + rb[227] + rb[228] + rb[229] + rb[
        230];
    a7_0__ = (rb[218] + rf[219] + rf[220] + rf[221] + rf[222] + rf[223] + rf[
            224] + rf[225] + rf[226] + rf[227] + rf[229] + rf[230] + rb[231])
        / den;
    a7_2__ = rf[228] / den;
    /*     H2CC */
    den = rf[213] + rf[214] + rf[215] + rf[368] + rf[369] + rb[170] + rb[196]
        + rb[203] + rb[352] + rb[357];
    a8_0__ = (rf[170] + rf[196] + rf[203] + rb[213] + rb[214] + rb[215] + rf[
            352] + rf[357] + rb[368] + rb[369]) / den;
    /*     C5H5O(2,4) */
    den = rf[414] + rf[415] + rb[398] + rb[399] + rb[400] + rb[413];
    a9_0__ = (rf[398] + rf[399] + rf[400] + rf[413] + rb[414] + rb[415]) /
        den;
    /*     A1C2HC2H2 */
    den = rf[440] + rf[441] + rb[439] + rb[442];
    a10_0__ = (rf[439] + rb[440]) / den;
    a10_11__ = (rb[441] + rf[442]) / den;
    /*     A1C2HC2H2u */
    den = rf[442] + rf[443] + rb[441] + rb[444];
    a11_0__ = (rb[443] + rf[444]) / den;
    a11_10__ = (rf[441] + rb[442]) / den;
    /*     A1C2H4 */
    den = rf[474] + rf[475] + rf[476] + rf[477] + rf[478] + rf[479] + rf[523]
        + rb[480] + rb[481] + rb[483] + rb[484] + rb[485] + rb[489];
    a12_0__ = (rb[474] + rb[475] + rb[476] + rb[477] + rb[478] + rb[479] + rf[
            480] + rf[481] + rf[483] + rf[484] + rf[485] + rf[489] + rb[523])
        / den;
    /*     A2C2H2 */
    den = rf[452] + rb[447] + rb[499] + rb[503] + rb[504];
    a13_0__ = (rf[447] + rb[452] + rf[499] + rf[503] + rf[504]) / den;

    a6_0__ += a6_5__ * a5_0__;
    den = 1 - a6_5__ * a5_6__;
    a6_0__ /= den;
    a6_2__ /= den;
    a2_0__ += a2_7__ * a7_0__;
    den = 1 - a2_7__ * a7_2__;
    a2_0__ /= den;
    a2_6__ /= den;
    a2_1__ /= den;
    a2_0__ += a2_1__ * a1_0__;
    den = 1 - a2_1__ * a1_2__;
    a2_0__ /= den;
    a2_6__ /= den;
    a2_0__ += a2_6__ * a6_0__;
    den = 1 - a2_6__ * a6_2__;
    a2_0__ /= den;
    xq[2] = a2_0__;
    xq[6] = a6_0__ + a6_2__ * xq[2];
    xq[1] = a1_0__ + a1_2__ * xq[2];
    xq[7] = a7_0__ + a7_2__ * xq[2];
    xq[5] = a5_0__ + a5_6__ * xq[6];
    xq[3] = a3_0__;
    xq[4] = a4_0__;
    xq[8] = a8_0__;
    xq[9] = a9_0__;
    a10_0__ += a10_11__ * a11_0__;
    den = 1 - a10_11__ * a11_10__;
    a10_0__ /= den;
    xq[10] = a10_0__;
    xq[11] = a11_0__ + a11_10__ * xq[10];
    xq[12] = a12_0__;
    xq[13] = a13_0__;

    rf[46] *= xq[2];
    rf[47] *= xq[2];
    rf[48] *= xq[2];
    rb[49] *= xq[2];
    rf[50] *= xq[2];
    rf[51] *= xq[2];
    rf[52] *= xq[2];
    rb[53] *= xq[1];
    rf[54] *= xq[1];
    rf[55] *= xq[1];
    rf[56] *= xq[1];
    rf[57] *= xq[1];
    rf[58] *= xq[1];
    rb[60] *= xq[1];
    rf[61] *= xq[1];
    rf[62] *= xq[1];
    rb[65] *= xq[3];
    rb[66] *= xq[1];
    rb[67] *= xq[1];
    rb[68] *= xq[2];
    rb[69] *= xq[1];
    rb[70] *= xq[1];
    rb[71] *= xq[2];
    rb[72] *= xq[1];
    rb[73] *= xq[1];
    rb[74] *= xq[1];
    rf[75] *= xq[2];
    rb[75] *= xq[1];
    rb[83] *= xq[3];
    rb[85] *= xq[1];
    rb[86] *= xq[2];
    rb[87] *= xq[2];
    rb[90] *= xq[2];
    rf[92] *= xq[3];
    rf[93] *= xq[3];
    rf[94] *= xq[3];
    rf[95] *= xq[3];
    rf[96] *= xq[3];
    rf[97] *= xq[3];
    rf[98] *= xq[3];
    rf[99] *= xq[3];
    rf[100] *= xq[3];
    rf[101] *= xq[3];
    rf[102] *= xq[3];
    rf[115] *= xq[2];
    rf[116] *= xq[3];
    rb[129] *= xq[6];
    rb[130] *= xq[6];
    rb[131] *= xq[6];
    rb[132] *= xq[6];
    rf[133] *= xq[2];
    rb[133] *= xq[6];
    rf[134] *= xq[6];
    rb[134] *= xq[5];
    rf[135] *= xq[6];
    rb[138] *= xq[5];
    rb[140] *= xq[5];
    rb[141] *= xq[5];
    rb[142] *= xq[5];
    rb[143] *= xq[5];
    rb[144] *= xq[5];
    rf[146] *= xq[5];
    rf[147] *= xq[5];
    rf[148] *= xq[5];
    rf[149] *= xq[5];
    rf[155] *= xq[5];
    rb[161] *= xq[1];
    rf[163] *= xq[3];
    rb[166] *= xq[3];
    rb[170] *= xq[8];
    rf[187] *= xq[2];
    rf[189] *= xq[3];
    rb[196] *= xq[8];
    rf[199] *= xq[4];
    rf[200] *= xq[4];
    rf[201] *= xq[4];
    rf[202] *= xq[4];
    rb[203] *= xq[8];
    rb[206] *= xq[4];
    rf[211] *= xq[3];
    rf[213] *= xq[8];
    rf[214] *= xq[8];
    rf[215] *= xq[8];
    rf[218] *= xq[7];
    rb[219] *= xq[7];
    rb[220] *= xq[7];
    rb[221] *= xq[7];
    rb[222] *= xq[7];
    rb[223] *= xq[7];
    rb[224] *= xq[7];
    rb[225] *= xq[7];
    rb[226] *= xq[7];
    rb[227] *= xq[7];
    rf[228] *= xq[2];
    rb[228] *= xq[7];
    rb[229] *= xq[7];
    rb[230] *= xq[7];
    rf[231] *= xq[7];
    rb[262] *= xq[5];
    rb[279] *= xq[5];
    rf[292] *= xq[4];
    rf[296] *= xq[4];
    rb[304] *= xq[4];
    rf[308] *= xq[4];
    rb[334] *= xq[5];
    rb[335] *= xq[5];
    rf[349] *= xq[4];
    rb[352] *= xq[8];
    rf[356] *= xq[4];
    rb[357] *= xq[8];
    rf[362] *= xq[4];
    rf[368] *= xq[8];
    rf[369] *= xq[8];
    rb[398] *= xq[9];
    rb[399] *= xq[9];
    rb[400] *= xq[9];
    rb[413] *= xq[9];
    rf[414] *= xq[9];
    rf[415] *= xq[9];
    rf[433] *= xq[4];
    rb[439] *= xq[10];
    rf[440] *= xq[10];
    rf[441] *= xq[10];
    rb[441] *= xq[11];
    rf[442] *= xq[11];
    rb[442] *= xq[10];
    rf[443] *= xq[11];
    rb[444] *= xq[11];
    rb[447] *= xq[13];
    rf[452] *= xq[13];
    rb[468] *= xq[4];
    rf[474] *= xq[12];
    rf[475] *= xq[12];
    rf[476] *= xq[12];
    rf[477] *= xq[12];
    rf[478] *= xq[12];
    rf[479] *= xq[12];
    rb[480] *= xq[12];
    rb[481] *= xq[12];
    rb[483] *= xq[12];
    rb[484] *= xq[12];
    rb[485] *= xq[12];
    rb[489] *= xq[12];
    rb[499] *= xq[13];
    rb[503] *= xq[13];
    rb[504] *= xq[13];
    rf[523] *= xq[12];
    rb[535] *= xq[4];

    return 0;
} /* qssa_ */

