/** getrates.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "c2h4red.h"

#include <cmath>
#include <algorithm>

using namespace std;

/* Table of constant values */

static double c_b5 = 10.;


/*     15-step Reduced Mechanism for C2H4/Air */

/*     August 25, 2005 */

/*     Developed by Tianfeng Lu */
/*     MAE Dept., Princeton Univ. */
/*     D103 E-Quad, Olden St */
/*     Princeton, NJ 08544 */
/*     Email: tlu@princeton.edu */

/*     Applicable parameter range: */
/*         Equivalence ratio: 0.5-1.7 (for premixed mixtures) */
/*         Atmospheric or slightly higer pressure */
/*         Ingition, extinction, and flames */
/*                                                                      C */
/* ----------------------------------------------------------------------C */

///* Subroutine */ int getrates_(double *p, double *t, double *y,
//    int *ickwrk, double *rckwrk, double *rrsp)

/* Subroutine */

// int getrates_(double *p, double *t, double *y, double *rrsp)

void c2h4red::getProblemSpecificRR(double rho, double temp, double pres,
                          double *yi, double *rrsp) const {

    // rho not used but fits the generic interface: units kg/m3
    // temp units are K
    // rrsp is output: units are kmol/m^3*s
    // pres input units are Pa

    pres *= 10.0;      // convert Pa to dynes/cm^2


    /* Initialized data */

    static double small = 1e-50;
    static double rf[167] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
    static double rb[167] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
    static double xq[10] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

    /* System generated locals */
    double d__1, d__2, d__3;

    /* Builtin functions */
    //dol double log(double);

    /* Local variables */
    static double c__[19];
    static int i__, k;
    static double a1, a2, a3, a4, a5, a6, a9, b10, a12, a21, a23, a24,
        a25, a32, a34, a35, a29, a42, a43, a45, a46, a49, a52, a53, a54,
        a64, a92, a94, b12, b13, b14, b20, b22, b23, b24, b30, b32, b33,
        b34, c12, c13, c10, c22, c23, c20, den[10], abv[10], rop[167],
        sum;
    extern /* Subroutine */ int ratt_(double *, double *, double *
        , double *, double *), ratx_(double *, double *,
        double *, double *, double *, double *);
    static double alogt, rklow[18];


    /* Parameter adjustments */
    --rrsp;
//    --rckwrk;
//    --ickwrk;
    --yi;

    /* Function Body */

/*     Convert mass fraction to mole concentration */

    c__[0] = yi[1] * .496046521;
    c__[1] = yi[2] * .992093043;
    c__[2] = yi[3] * .0625023433;
    c__[3] = yi[4] * .0312511716;
    c__[4] = yi[5] * .0587980383;
    c__[5] = yi[6] * .0555082499;
    c__[6] = yi[7] * .0302968146;
    c__[7] = yi[8] * .0293990192;
    c__[8] = yi[9] * .0665112065;
    c__[9] = yi[10] * .0623323639;
    c__[10] = yi[11] * .0357008335;
    c__[11] = yi[12] * .0227221341;
    c__[12] = yi[13] * .0333039255;
    c__[13] = yi[14] * .0384050525;
    c__[14] = yi[15] * .0356453112;
    c__[15] = yi[16] * .0332556033;
    c__[16] = yi[17] * .0237882046;
    c__[17] = yi[18] * .0237635408;
    c__[18] = yi[19] * .0356972032;

    sum = 0.f;
    for (k = 1; k <= 19; ++k) {
    sum += c__[k - 1];
    }
    sum = pres / (sum * temp * 83145100.);

    for (k = 1; k <= 19; ++k) {
/* Computing MAX */
    d__1 = c__[k - 1];
    c__[k - 1] = max(d__1,small);
    c__[k - 1] *= sum;
    }

/*     Compute T- and P- dependent elementary reaction rates */

    alogt = log(temp);
    ratt_(&temp, &alogt, rf, rb, rklow);
    ratx_(&temp, &alogt, c__, rf, rb, rklow);

/*     Compute QSS species concetrations */

/*     C */
    abv[0] = rb[72] + rb[97] + rb[98];
    den[0] = rb[42] + rf[72] + rf[97] + rf[98];
/*     CH */
    abv[1] = rb[5] + rb[101] + rb[102] + rb[104] + rb[107] + rb[108] + rb[149]
        ;
    den[1] = rf[5] + rf[42] + rb[44] + rf[73] + rb[75] + rf[99] + rf[100] +
        rf[101] + rf[104] + rf[105] + rf[106] + rf[107] + rf[149];
/*     CH2 */
    abv[2] = rf[17] + rf[24] + rb[43] + rb[74] + rf[77] + rb[92] + rb[102] +
        rb[109] + rb[110] + rb[111] + rb[111] + rb[112] + rb[113] + rb[
        114] + rb[150] + rb[151] + rb[152] + rb[152];
    den[2] = rf[6] + rb[17] + rb[24] + rf[43] + rf[74] + rf[75] + rb[77] + rf[
        92] + rb[100] + rf[109] + rf[110] + rf[112] + rf[113] + rf[114] +
        rb[116] + rb[120] + rb[123] + rb[124] + rf[150] + rf[151];
/*     CH2 */
    abv[3] = rb[7] + rb[76] + rf[78] + rb[117] + rb[118] + rb[119] + rb[121]
        + rb[122] + rb[125] + rb[153];
    den[3] = rf[7] + rf[8] + rf[44] + rb[53] + rb[62] + rf[76] + rb[78] + rf[
        116] + rf[117] + rf[118] + rf[119] + rf[120] + rf[121] + rf[122]
        + rf[123] + rf[124] + rf[125] + rf[153];
/*     HCO */
    abv[4] = rb[12] + rb[13] + rf[14] + rf[19] + rf[26] + rb[47] + rb[48] +
        rf[50] + rb[81] + rf[82] + rf[96] + rb[132] + rf[133] + rb[136] +
        rb[137] + rb[138] + rf[165];
    den[4] = rb[6] + rb[8] + rf[12] + rf[13] + rb[14] + rb[19] + rb[26] + rf[
        47] + rf[48] + rb[50] + rb[73] + rf[81] + rb[82] + rb[96] + rb[99]
         + rb[106] + rf[132] + rb[133] + rf[136] + rf[137] + rf[138] + rb[
        160];
/*     CH3O */
    abv[5] = rb[15] + rf[49] + rb[51] + rb[52] + rb[83] + rf[94] + rf[127] +
        rb[139];
    den[5] = rf[15] + rb[49] + rf[51] + rf[52] + rf[53] + rf[83] + rb[94] +
        rb[127] + rf[139];
/*     C2H3 */
    abv[6] = rb[18] + rf[54] + rb[55] + rb[56] + rf[58] + rb[86] + rf[87] +
        rf[134] + rb[155] + rb[166];
    den[6] = rf[18] + rb[54] + rf[55] + rf[56] + rb[58] + rf[86] + rb[87] +
        rb[134] + rf[140] + rf[154] + rf[155] + rf[166];
/*     C2H5 */
    abv[7] = rb[20] + rf[21] + rf[57] + rb[59] + rb[60] + rf[61] + rf[88] +
        rf[131] + rf[135] + rb[142] + rf[165];
    den[7] = rf[20] + rb[21] + rb[57] + rf[59] + rf[60] + rb[61] + rb[88] +
        rb[131] + rb[135] + rf[142];
/*     HCCO */
    abv[8] = rf[16] + rb[22] + rf[23] + rf[63] + rf[89] + rb[108] + rb[143] +
        rb[144] + rb[144];
    den[8] = rb[16] + rf[22] + rb[23] + rf[62] + rb[63] + rb[89] + rb[105] +
        rf[143];
/*     CH2CHO */
    abv[9] = rf[146] + rf[156] + rb[158] + rb[161] + rb[162];
    den[9] = rb[146] + rb[156] + rf[158] + rf[161] + rf[162];

/*     C2H3 */
    xq[6] = abv[6] / den[6];
/*     C2H5 */
    xq[7] = abv[7] / den[7];

    a1 = abv[0] / den[0];
    a12 = rf[42] / den[0];
    a2 = abv[1] / den[1];
    a21 = rb[42] / den[1];
    a23 = (rf[75] + rb[100]) / den[1];
    a24 = rf[44] / den[1];
    a25 = (rb[73] + rb[99] + rb[106]) / den[1];
    a29 = rb[105] / den[1];
    a3 = abv[2] / den[2];
    a32 = (rb[75] + rf[100]) / den[2];
    a34 = (rf[116] + rf[120] + rf[123] + rf[124]) / den[2];
    a35 = rb[6] / den[2];
    a4 = abv[3] / den[3];
    a42 = rb[44] / den[3];
    a43 = (rb[116] + rb[120] + rb[123] + rb[124]) / den[3];
    a45 = rb[8] / den[3];
    a46 = rf[53] / den[3];
    a49 = rf[62] / den[3];
    a5 = (abv[4] + rf[140] * xq[6]) / den[4];
    a52 = (rf[73] + rf[99] + rf[106]) / den[4];
    a53 = rf[6] / den[4];
    a54 = rf[8] / den[4];
    a6 = abv[5] / den[5];
    a64 = rb[53] / den[5];
    a9 = abv[8] / den[8];
    a92 = rf[105] / den[8];
    a94 = rb[62] / den[8];

    b10 = a2 + a21 * a1 + a5 * a25 + a9 * a29;
    b12 = a21 * a12 + a25 * a52 + a29 * a92 - 1;
    b13 = a23 + a25 * a53;
    b14 = a24 + a25 * a54 + a29 * a94;
    b20 = a3 + a35 * a5;
    b22 = a32 + a35 * a52;
    b23 = a35 * a53 - 1;
    b24 = a34 + a35 * a54;
    b30 = a4 + a45 * a5 + a46 * a6 + a49 * a9;
    b32 = a42 + a45 * a52 + a49 * a92;
    b33 = a43 + a45 * a53;
    b34 = a45 * a54 + a46 * a64 + a49 * a94 - 1;
    c12 = b34 * b12 - b14 * b32;
    c13 = b34 * b13 - b14 * b33;
    c10 = b34 * b10 - b14 * b30;
    c22 = b34 * b22 - b24 * b32;
    c23 = b34 * b23 - b24 * b33;
    c20 = b34 * b20 - b24 * b30;

/*     CH */
/* Computing MAX */
    d__3 = (d__2 = c23 * c12 - c13 * c22, abs(d__2));
    xq[1] = (d__1 = c13 * c20 - c23 * c10, abs(d__1)) / max(d__3,small);
/*     CH2 */
/* Computing MAX */
    d__2 = abs(c23);
    xq[2] = (d__1 = c22 * xq[1] + c20, abs(d__1)) / max(d__2,small);
/*     CH2* */
/* Computing MAX */
    d__2 = abs(b34);
    xq[3] = (d__1 = b32 * xq[1] + b33 * xq[2] + b30, abs(d__1)) / max(d__2,
        small);
/*     C */
    xq[0] = a1 + a12 * xq[1];
/*     HCO */
    xq[4] = a5 + a52 * xq[1] + a53 * xq[2] + a54 * xq[3];
/*     CH3O */
    xq[5] = a6 + a64 * xq[3];
/*     HCCO */
    xq[8] = a9 + a92 * xq[1] + a94 * xq[3];
/*     CH2CHO */
    xq[9] = (abv[9] + rb[160] * xq[4] + rf[154] * xq[6]) / den[9];

/*     Multiply QSS species concetratios to involved reactions */

    rf[5] *= xq[1];
    rf[6] *= xq[2];
    rb[6] *= xq[4];
    rf[7] *= xq[3];
    rf[8] *= xq[3];
    rb[8] *= xq[4];
    rf[12] *= xq[4];
    rf[13] *= xq[4];
    rb[14] *= xq[4];
    rf[15] *= xq[5];
    rb[16] *= xq[8];
    rb[17] *= xq[2];
    rf[18] *= xq[6];
    rb[19] *= xq[4];
    rf[20] *= xq[7];
    rb[21] *= xq[7];
    rf[22] *= xq[8];
    rb[23] *= xq[8];
    rb[24] *= xq[2];
    rb[26] *= xq[4];
    rf[42] *= xq[1];
    rb[42] *= xq[0];
    rf[43] *= xq[2];
    rf[44] *= xq[3];
    rb[44] *= xq[1];
    rf[47] *= xq[4];
    rf[48] *= xq[4];
    rb[49] *= xq[5];
    rb[50] *= xq[4];
    rf[51] *= xq[5];
    rf[52] *= xq[5];
    rf[53] *= xq[5];
    rb[53] *= xq[3];
    rb[54] *= xq[6];
    rf[55] *= xq[6];
    rf[56] *= xq[6];
    rb[57] *= xq[7];
    rb[58] *= xq[6];
    rf[59] *= xq[7];
    rf[60] *= xq[7];
    rb[61] *= xq[7];
    rf[62] *= xq[8];
    rb[62] *= xq[3];
    rb[63] *= xq[8];
    rf[72] *= xq[0];
    rf[73] *= xq[1];
    rb[73] *= xq[4];
    rf[74] *= xq[2];
    rf[75] *= xq[2];
    rb[75] *= xq[1];
    rf[76] *= xq[3];
    rb[77] *= xq[2];
    rb[78] *= xq[3];
    rf[81] *= xq[4];
    rb[82] *= xq[4];
    rf[83] *= xq[5];
    rf[86] *= xq[6];
    rb[87] *= xq[6];
    rb[88] *= xq[7];
    rb[89] *= xq[8];
    rf[92] *= xq[2];
    rb[94] *= xq[5];
    rb[96] *= xq[4];
    rf[97] *= xq[0];
    rf[98] *= xq[0];
    rf[99] *= xq[1];
    rb[99] *= xq[4];
    rf[100] *= xq[1];
    rb[100] *= xq[2];
    rf[101] *= xq[1];
    rf[104] *= xq[1];
    rf[105] *= xq[1];
    rb[105] *= xq[8];
    rf[106] *= xq[1];
    rb[106] *= xq[4];
    rf[107] *= xq[1];
    rf[109] *= xq[2];
    rf[110] *= xq[2];
    rf[112] *= xq[2];
    rf[113] *= xq[2];
    rf[114] *= xq[2];
    rf[116] *= xq[3];
    rb[116] *= xq[2];
    rf[117] *= xq[3];
    rf[118] *= xq[3];
    rf[119] *= xq[3];
    rf[120] *= xq[3];
    rb[120] *= xq[2];
    rf[121] *= xq[3];
    rf[122] *= xq[3];
    rf[123] *= xq[3];
    rb[123] *= xq[2];
    rf[124] *= xq[3];
    rb[124] *= xq[2];
    rf[125] *= xq[3];
    rb[127] *= xq[5];
    rb[131] *= xq[7];
    rf[132] *= xq[4];
    rb[133] *= xq[4];
    rb[134] *= xq[6];
    rb[135] *= xq[7];
    rf[136] *= xq[4];
    rf[137] *= xq[4];
    rf[138] *= xq[4];
    rf[139] *= xq[5];
    rf[140] *= xq[6];
    rf[142] *= xq[7];
    rf[143] *= xq[8];
    rb[146] *= xq[9];
    rf[149] *= xq[1];
    rf[150] *= xq[2];
    rf[151] *= xq[2];
    rf[153] *= xq[3];
    rf[154] *= xq[6];
    rf[155] *= xq[6];
    rb[156] *= xq[9];
    rf[158] *= xq[9];
    rb[160] *= xq[4];
    rf[161] *= xq[9];
    rf[162] *= xq[9];
    rf[166] *= xq[6];

    for (i__ = 1; i__ <= 167; ++i__) {
    rop[i__ - 1] = rf[i__ - 1] - rb[i__ - 1];
    }

    rop[27] += rop[28];
    rop[27] += rop[29];
    rop[27] += rop[30];
    rop[32] += rop[33];
    rop[32] += rop[34];
    rop[32] += rop[35];
    rop[69] += rop[147];
    rop[70] += rop[71];
    rop[90] += rop[91];
    rop[116] += rop[120];
    rop[116] += rop[123];
    rop[116] += rop[124];
    rop[136] += rop[137];

    rrsp[1] = -rop[2] + rop[7] + rop[32] + rop[38] + rop[40] + rop[42] + rop[
        44] + rop[46] + rop[48] + rop[50] + rop[51] + rop[56] + rop[58] +
        rop[60] + rop[61] + rop[63] - rop[65] - rop[66] - rop[100] - rop[
        110] + rop[111] - rop[119] + rop[141] + rop[145] + rop[148] - rop[
        149] + rop[153] + rop[161];
    rrsp[2] = -rop[1] + rop[2] + rop[5] + rop[6] + rop[8] + rop[9] + rop[13]
        + rop[16] + rop[18] + rop[22] - rop[27] - rop[31] - rop[32] - rop[
        32] - rop[36] - rop[37] - rop[38] - rop[39] - rop[40] - rop[41] -
        rop[42] - rop[43] - rop[44] - rop[45] - rop[46] - rop[47] - rop[
        48] - rop[49] - rop[50] - rop[51] - rop[52] - rop[53] - rop[54] -
        rop[55] - rop[56] - rop[57] - rop[58] - rop[59] - rop[60] - rop[
        61] - rop[62] - rop[63] - rop[64] + rop[66] + rop[72] + rop[73] +
        rop[74] + rop[76] + rop[80] + rop[84] + rop[98] + rop[100] + rop[
        101] + rop[102] + rop[104] + rop[107] + rop[109] + rop[110] + rop[
        112] + rop[117] + rop[119] + rop[121] + rop[131] + rop[136] + rop[
        145] + rop[146] + rop[150] + rop[150] + rop[152] + rop[152] - rop[
        156] - rop[160] - rop[161] - rop[163] + rop[164];
    rrsp[3] = -rop[0] - rop[0] - rop[1] - rop[2] - rop[3] - rop[4] - rop[5] -
        rop[6] - rop[7] - rop[8] - rop[9] - rop[10] - rop[11] - rop[12] -
        rop[13] - rop[14] - rop[15] - rop[16] - rop[17] - rop[18] - rop[
        19] - rop[20] - rop[21] - rop[22] - rop[23] - rop[24] + rop[25] +
        rop[31] + rop[37] + rop[68] + rop[97] + rop[99] + rop[127] - rop[
        145] - rop[146] + rop[151] + rop[154] - rop[164] - rop[165];
    rrsp[4] = rop[0] + rop[3] - rop[25] - rop[26] - rop[27] - rop[31] + rop[
        38] + rop[69] + rop[90] + rop[93] - rop[97] - rop[99] - rop[109]
        - rop[117] - rop[118] - rop[127] - rop[128] - rop[138] - rop[139]
        - rop[140] - rop[142] - rop[143] - rop[150] - rop[151] - rop[154]
        - rop[155] - rop[158];
    rrsp[5] = rop[1] + rop[2] + rop[3] + rop[4] + rop[10] + rop[12] + rop[14]
        + rop[15] + rop[21] + rop[23] + rop[31] - rop[36] + rop[39] + rop[
        39] + rop[41] + rop[52] - rop[66] - rop[67] - rop[67] - rop[68] -
        rop[68] - rop[69] - rop[70] - rop[72] - rop[73] - rop[74] - rop[
        75] - rop[76] - rop[77] - rop[78] - rop[79] - rop[80] - rop[81] -
        rop[82] - rop[83] - rop[84] - rop[85] - rop[86] - rop[87] - rop[
        88] - rop[89] + rop[92] + rop[94] + rop[95] + rop[109] + rop[117]
        + rop[128] + rop[143] - rop[148] + rop[158] - rop[162];
    rrsp[6] = rop[36] + rop[37] + rop[41] + rop[53] + rop[66] + rop[68] + rop[
        69] + rop[70] + rop[75] + rop[77] + rop[78] + rop[79] + rop[81] +
        rop[82] + rop[83] + rop[86] + rop[87] + rop[88] + rop[89] - rop[
        101] + rop[118] - rop[153] + rop[162];
    rrsp[7] = -rop[3] + rop[4] + rop[26] + rop[27] - rop[37] - rop[38] - rop[
        39] + rop[40] - rop[69] + rop[70] - rop[90] - rop[90] - rop[92] -
        rop[93] - rop[94] - rop[95] - rop[96] + rop[129] + rop[138] + rop[
        139] + rop[142] + rop[155];
    rrsp[8] = -rop[4] - rop[40] - rop[41] + rop[67] - rop[70] + rop[90] + rop[
        96] - rop[129];
    rrsp[9] = -rop[9] + rop[10] + rop[19] + rop[20] + rop[43] - rop[45] + rop[
        46] + rop[52] + rop[64] - rop[77] - rop[78] + rop[79] + rop[85] -
        rop[93] - rop[94] - rop[98] + rop[110] - rop[112] + rop[113] +
        rop[113] + rop[119] - rop[121] + rop[122] + rop[122] - rop[127] -
        rop[128] - rop[129] - rop[130] - rop[130] - rop[131] - rop[131] -
        rop[132] - rop[133] - rop[134] - rop[135] - rop[145] - rop[148] +
        rop[149] + rop[160] + rop[163] + rop[164] - rop[166];
    rrsp[10] = -rop[10] + rop[45] - rop[46] - rop[79] + rop[93] - rop[104] -
        rop[113] - rop[122] + rop[129] + rop[132] + rop[133] + rop[134] +
        rop[135];
    rrsp[11] = rop[5] + rop[7] - rop[11] + rop[12] + rop[17] + rop[22] + rop[
        22] - rop[25] + rop[48] + rop[62] + rop[64] - rop[65] + rop[72] -
        rop[80] + rop[81] + rop[85] - rop[95] + rop[97] - rop[105] + rop[
        106] + rop[108] + rop[109] - rop[114] + rop[117] + rop[118] + rop[
        125] + rop[132] + rop[136] + rop[138] + rop[143] + rop[143] + rop[
        144] + rop[144] + rop[145] + rop[158];
    rrsp[12] = rop[11] + rop[13] + rop[24] + rop[25] + rop[80] + rop[95] -
        rop[106] - rop[125] + rop[150];
    rrsp[13] = rop[9] - rop[14] + rop[15] + rop[20] - rop[26] + rop[47] - rop[
        49] - rop[50] + rop[51] + rop[65] + rop[74] + rop[76] - rop[82] +
        rop[83] + rop[92] - rop[96] + rop[101] - rop[107] + rop[125] +
        rop[128] - rop[133] + rop[139] + rop[140] + rop[148] + rop[151] +
        rop[153] + rop[158];
    rrsp[14] = -rop[16] - rop[17] - rop[54] + rop[56] - rop[84] - rop[85] +
        rop[86] + rop[98] + rop[102] + rop[108] + rop[111] + rop[141] +
        rop[144] + rop[152] + rop[155];
    rrsp[15] = -rop[19] + rop[55] - rop[57] - rop[58] + rop[60] - rop[87] +
        rop[104] + rop[112] + rop[121] - rop[134] - rop[141] + rop[142] -
        rop[146] + rop[163];
    rrsp[16] = -rop[21] + rop[59] - rop[61] - rop[88] + rop[130] - rop[135];
    rrsp[17] = rop[18] - rop[23] - rop[24] - rop[63] - rop[64] + rop[84] -
        rop[89] + rop[107] + rop[114] - rop[156] + rop[161] + rop[162] +
        rop[164];
    rrsp[18] = -rop[163] - rop[164] - rop[165] + rop[166];
    rrsp[19] = 0.f;


    // convert units from mol/cm3*s to kmol/m3*s

    for(int k=1; k<=19; k++)
        rrsp[k] *= 1000.0;

}


/* ----------------------------------------------------------------------C */

/* Subroutine */ int ratt_(double *t, double *alogt, double *rf,
    double *rb, double *rklow)
{
    /* Initialized data */

    static double small = 1e-50;

    /* Builtin functions */
    //dol double exp(double);

    /* Local variables */
    static int n;
    static double eg[29], ti, ti2, eqk[167], smh[29], tmp, pfac, pfac2,
        pfac3;
    extern /* Subroutine */ int rdsmh_(double *, double *, double
        *);


/* *****precision > double */
/* *****END precision > double */

    /* Parameter adjustments */
    --rklow;
    --rb;
    --rf;

    /* Function Body */

    ti = 1. / *t;
    ti2 = ti * ti;

    rf[1] = ti * 1.2e17;
    rf[2] = ti * 5e17;
    rf[3] = exp(*alogt * 2.7 + 10.5635949 - ti * 3150.13634);
    rf[4] = 2e13;
    rf[5] = exp(*alogt * 2. + 16.0803938 - ti * 2012.86667);
    rf[6] = 5.7e13;
    rf[7] = 8e13;
    rf[8] = 1.5e13;
    rf[9] = 1.5e13;
    rf[10] = 5.06e13;
    rf[11] = exp(*alogt * 1.5 + 20.7430685 - ti * 4327.66334);
    rf[12] = exp(23.6136376 - ti * 1200.17175);
    rf[13] = 3e13;
    rf[14] = 3e13;
    rf[15] = exp(31.2945828 - ti * 1781.387);
    rf[16] = 1e13;
    tmp = exp(*alogt * 2. - ti * 956.111669);
    rf[17] = tmp * 1.35e7;
    rf[18] = tmp * 6.94e6;
    rf[19] = 3e13;
    tmp = exp(*alogt * 1.83 - ti * 110.707667);
    rf[20] = tmp * 2.5e7;
    rf[147] = tmp * 3.35e6;
    rf[21] = 2.24e13;
    rf[22] = exp(*alogt * 1.92 + 18.3130955 - ti * 2863.30284);
    rf[23] = 1e14;
    tmp = exp(ti * -4025.73334);
    rf[24] = tmp * 1e13;
    rf[64] = tmp * 5e13;
    rf[25] = exp(28.1906369 - ti * 679.342501);
    rf[26] = exp(28.5473118 - ti * 24053.7567);
    rf[27] = exp(32.2361913 - ti * 20128.6667);
    rf[28] = exp(42.4761511 - *alogt * .86);
    tmp = exp(*alogt * -1.24);
    rf[29] = tmp * 2.08e19;
    rf[31] = tmp * 2.6e19;
    rf[30] = exp(43.8677883 - *alogt * .76);
    rf[32] = exp(37.8159211 - *alogt * .6707 - ti * 8575.31523);
    rf[33] = ti * 1e18;
    rf[34] = exp(39.0385861 - *alogt * .6);
    rf[35] = exp(45.5408762 - *alogt * 1.25);
    rf[36] = ti2 * 5.5e20;
    rf[37] = ti2 * 2.2e22;
    rf[38] = exp(29.0097872 - ti * 337.658384);
    rf[39] = exp(31.9812991 - ti * 537.435401);
    rf[40] = exp(31.3686907 - ti * 319.542584);
    rf[41] = exp(*alogt * 2. + 16.308716 - ti * 2616.72667);
    rf[42] = exp(29.9336062 - ti * 1811.58);
    rf[43] = 1.65e14;
    rf[44] = 6e14;
    rf[45] = 3e13;
    rf[46] = exp(37.1706652 - *alogt * .534 - ti * 269.724134);
    rf[47] = exp(*alogt * 1.62 + 20.3077504 - ti * 5454.86868);
    rf[48] = exp(*alogt * .48 + 27.7171988 + ti * 130.836334);
    rf[49] = 7.34e13;
    rf[50] = exp(*alogt * .454 + 27.014835 - ti * 1308.36334);
    rf[51] = exp(*alogt * 1.9 + 17.8655549 - ti * 1379.8201);
    rf[52] = 2e13;
    rf[53] = exp(*alogt * .5 + 28.0364862 + ti * 55.3538334);
    rf[54] = exp(33.1993656 - *alogt * .23 - ti * 538.441834);
    rf[55] = exp(29.3537877 - ti * 1207.72);
    rf[56] = exp(*alogt * .27 + 29.4360258 - ti * 140.900667);
    rf[57] = 3.03e13;
    rf[58] = exp(*alogt * .454 + 27.3863985 - ti * 915.854335);
    rf[59] = exp(*alogt * 2.53 + 14.096923 - ti * 6159.37201);
    rf[60] = exp(40.7945264 - *alogt * .99 - ti * 795.082335);
    rf[61] = 2e12;
    rf[62] = exp(*alogt * 1.9 + 18.5604427 - ti * 3789.22151);
    rf[63] = 1e14;
    rf[65] = exp(30.0558238 - ti * 1725.02674);
    rf[66] = exp(*alogt * 1.5 + 17.5767107 - ti * 40056.0467);
    rf[67] = exp(*alogt * 1.51 + 19.190789 - ti * 1726.03317);
    rf[68] = exp(31.9350862 - *alogt * .37);
    rf[69] = exp(*alogt * 2.4 + 10.482906 + ti * 1061.78717);
    rf[70] = exp(ti * 251.608334 + 30.3051698);
    rf[71] = exp(28.3241683 - ti * 214.873517);
    rf[72] = exp(41.9771599 - ti * 14799.6022);
    rf[73] = 5e13;
    rf[74] = 3e13;
    rf[75] = 2e13;
    rf[76] = exp(*alogt * 2. + 16.2403133 - ti * 1509.65);
    rf[77] = 3e13;
    rf[78] = exp(*alogt * 1.6 + 17.8408622 - ti * 2727.43434);
    rf[79] = exp(41.0064751 - *alogt * 1.34 - ti * 713.058018);
    rf[80] = exp(*alogt * 1.6 + 18.4206807 - ti * 1570.036);
    rf[81] = exp(*alogt * 1.228 + 17.6783433 - ti * 35.2251667);
    rf[82] = 5e13;
    rf[83] = exp(*alogt * 1.18 + 21.9558261 + ti * 224.93785);
    rf[84] = 5e12;
    rf[85] = exp(*alogt * 4.5 - 8.4310155 + ti * 503.216668);
    rf[86] = exp(*alogt * 4. - 7.6354939 + ti * 1006.43334);
    rf[87] = 5e12;
    rf[88] = exp(*alogt * 2. + 14.4032972 - ti * 1258.04167);
    rf[89] = exp(*alogt * 2.12 + 15.0796373 - ti * 437.798501);
    rf[90] = exp(29.6459241 - ti * 1006.43334);
    rf[91] = exp(ti * 820.243168 + 25.5908003);
    rf[92] = exp(33.6712758 - ti * 6038.60001);
    rf[93] = 2e13;
    rf[94] = 1e12;
    rf[95] = 2.87e13;
    rf[96] = exp(32.6416564 - ti * 11875.9134);
    rf[97] = exp(*alogt * 2. + 15.5382772 - ti * 6038.60001);
    rf[98] = exp(31.6914641 - ti * 289.852801);
    rf[99] = 5e13;
    rf[100] = 6.71e13;
    rf[101] = exp(32.3131523 - ti * 1565.00384);
    rf[102] = exp(ti * 379.928584 + 29.3732401);
    rf[103] = 4e13;
    rf[105] = 6e13;
    rf[106] = 5e13;
    rf[107] = exp(32.8780452 - ti * 7946.79762);
    rf[108] = exp(ti * 259.156584 + 32.1806786);
    rf[109] = 5e13;
    tmp = exp(ti * -754.825001);
    rf[110] = tmp * 5e12;
    rf[151] = tmp * 5.8e12;
    rf[152] = tmp * 2.4e12;
    rf[111] = exp(*alogt * 2. + 13.1223634 - ti * 3638.25651);
    rf[112] = exp(35.00878 - ti * 6010.41988);
    rf[113] = 4e13;
    rf[114] = exp(*alogt * 2. + 14.7156719 - ti * 4161.60184);
    rf[115] = exp(*alogt * .5 + 27.4203001 - ti * 2269.50717);
    rf[117] = exp(30.3390713 - ti * 301.930001);
    rf[118] = 2.8e13;
    rf[119] = 1.2e13;
    rf[120] = 7e13;
    rf[121] = 3e13;
    tmp = exp(ti * 286.833501);
    rf[122] = tmp * 1.2e13;
    rf[123] = tmp * 1.6e13;
    rf[124] = 9e12;
    rf[125] = 7e12;
    rf[126] = 1.4e13;
    rf[128] = exp(31.2033668 - ti * 15338.044);
    rf[129] = exp(28.4682686 - ti * 10222.8466);
    rf[130] = exp(*alogt * 2.47 + 10.1064284 - ti * 2606.66234);
    rf[131] = exp(38.7538626 - *alogt * 1.18 - ti * 329.103701);
    rf[132] = exp(*alogt * .1 + 29.5538088 - ti * 5334.09668);
    rf[133] = 2.648e13;
    rf[134] = exp(*alogt * 2.81 + 8.10772006 - ti * 2948.84967);
    rf[135] = exp(*alogt * 2. + 12.3327053 - ti * 4629.59334);
    rf[136] = exp(*alogt * 1.74 + 15.6303353 - ti * 5258.61418);
    tmp = exp(*alogt * -1. - ti * 8554.68335);
    rf[137] = tmp * 1.5e18;
    rf[138] = tmp * 1.87e17;
    rf[139] = exp(30.2300002 - ti * 201.286667);
    rf[140] = exp(*alogt * 7.6 - 28.4796532 + ti * 1776.35484);
    rf[141] = exp(38.3630605 - *alogt * 1.39 - ti * 510.764918);
    rf[142] = exp(*alogt * .44 + 29.7104627 - ti * 43664.1103);
    rf[143] = exp(27.4566677 - ti * 1949.96459);
    rf[144] = exp(28.7941719 - ti * 429.747034);
    rf[145] = 1e13;
    rf[146] = 3.37e13;
    rf[148] = exp(36.1482143 - ti * 8720.74485);
    rf[149] = exp(*alogt * .5 + 22.8027074 + ti * 883.145252);
    rf[150] = exp(*alogt * .43 + 28.3090547 + ti * 186.190167);
    rf[154] = exp(*alogt * .25 + 24.9457104 + ti * 470.507584);
    rf[155] = exp(*alogt * .29 + 25.5207079 - ti * 5.53538334);
    rf[156] = exp(*alogt * 1.61 + 14.1059389 + ti * 193.2352);
    rf[157] = exp(*alogt * .422 + 26.9105027 + ti * 883.145252);
    rf[159] = 1.81e10;
    rf[161] = 2.2e13;
    rf[162] = 1.1e13;
    rf[163] = 1.2e13;
    rf[164] = exp(51.6031099 - *alogt * 2.39 - ti * 5625.96234);
    rf[165] = exp(*alogt * 1.65 + 18.6030023 - ti * 164.55185);
    rf[166] = exp(*alogt * 1.65 + 17.3708586 + ti * 489.126601);
    rf[167] = 2.5e13;

    rdsmh_(t, alogt, smh);
    for (n = 1; n <= 28; ++n) {
    eg[n - 1] = exp(smh[n - 1]);
    }

    pfac = 1013250. / (*t * 83145100.);
    pfac2 = pfac * pfac;
    pfac3 = pfac2 * pfac;

    eqk[0] = eg[3] / eg[2] / eg[2] / pfac;
    eqk[1] = eg[4] / eg[1] / eg[2] / pfac;
    eqk[2] = eg[1] * eg[4] / eg[0] / eg[2];
    eqk[3] = eg[3] * eg[4] / eg[2] / eg[6];
    eqk[4] = eg[4] * eg[6] / eg[2] / eg[7];
    eqk[5] = eg[1] * eg[14] / eg[2] / eg[9];
    eqk[6] = eg[1] * eg[16] / eg[2] / eg[10];
    eqk[7] = eg[0] * eg[14] / eg[2] / eg[11];
    eqk[8] = eg[1] * eg[16] / eg[2] / eg[11];
    eqk[9] = eg[1] * eg[17] / eg[2] / eg[12];
    eqk[10] = eg[4] * eg[12] / eg[2] / eg[13];
    eqk[11] = eg[15] / eg[2] / eg[14] / pfac;
    eqk[12] = eg[4] * eg[14] / eg[2] / eg[16];
    eqk[13] = eg[1] * eg[15] / eg[2] / eg[16];
    eqk[14] = eg[4] * eg[16] / eg[2] / eg[17];
    eqk[15] = eg[4] * eg[17] / eg[2] / eg[18];
    eqk[16] = eg[1] * eg[24] / eg[2] / eg[19];
    eqk[17] = eg[10] * eg[14] / eg[2] / eg[19];
    eqk[18] = eg[1] * eg[25] / eg[2] / eg[20];
    eqk[19] = eg[12] * eg[16] / eg[2] / eg[21];
    eqk[20] = eg[12] * eg[17] / eg[2] / eg[22];
    eqk[21] = eg[4] * eg[22] / eg[2] / eg[23];
    eqk[22] = eg[1] * eg[14] * eg[14] / eg[2] / eg[24] * pfac;
    eqk[23] = eg[4] * eg[24] / eg[2] / eg[25];
    eqk[24] = eg[10] * eg[15] / eg[2] / eg[25];
    eqk[25] = eg[2] * eg[15] / eg[3] / eg[14];
    eqk[26] = eg[6] * eg[16] / eg[3] / eg[17];
    eqk[27] = eg[6] / eg[1] / eg[3] / pfac;
    eqk[28] = eqk[27];
    eqk[29] = eqk[27];
    eqk[30] = eqk[27];
    eqk[31] = eg[2] * eg[4] / eg[1] / eg[3];
    eqk[32] = eg[0] / eg[1] / eg[1] / pfac;
    eqk[33] = eqk[32];
    eqk[34] = eqk[32];
    eqk[35] = eqk[32];
    eqk[36] = eg[5] / eg[1] / eg[4] / pfac;
    eqk[37] = eg[2] * eg[5] / eg[1] / eg[6];
    eqk[38] = eg[0] * eg[3] / eg[1] / eg[6];
    eqk[39] = eg[4] * eg[4] / eg[1] / eg[6];
    eqk[40] = eg[0] * eg[6] / eg[1] / eg[7];
    eqk[41] = eg[4] * eg[5] / eg[1] / eg[7];
    eqk[42] = eg[0] * eg[8] / eg[1] / eg[9];
    eqk[43] = eg[12] / eg[1] / eg[10] / pfac;
    eqk[44] = eg[0] * eg[9] / eg[1] / eg[11];
    eqk[45] = eg[13] / eg[1] / eg[12] / pfac;
    eqk[46] = eg[0] * eg[12] / eg[1] / eg[13];
    eqk[47] = eg[17] / eg[1] / eg[16] / pfac;
    eqk[48] = eg[0] * eg[14] / eg[1] / eg[16];
    eqk[49] = eg[18] / eg[1] / eg[17] / pfac;
    eqk[50] = eg[0] * eg[16] / eg[1] / eg[17];
    eqk[51] = eg[0] * eg[17] / eg[1] / eg[18];
    eqk[52] = eg[4] * eg[12] / eg[1] / eg[18];
    eqk[53] = eg[5] * eg[11] / eg[1] / eg[18];
    eqk[54] = eg[20] / eg[1] / eg[19] / pfac;
    eqk[55] = eg[21] / eg[1] / eg[20] / pfac;
    eqk[56] = eg[0] * eg[19] / eg[1] / eg[20];
    eqk[57] = eg[22] / eg[1] / eg[21] / pfac;
    eqk[58] = eg[0] * eg[20] / eg[1] / eg[21];
    eqk[59] = eg[23] / eg[1] / eg[22] / pfac;
    eqk[60] = eg[0] * eg[21] / eg[1] / eg[22];
    eqk[61] = eg[0] * eg[22] / eg[1] / eg[23];
    eqk[62] = eg[11] * eg[14] / eg[1] / eg[24];
    eqk[63] = eg[0] * eg[24] / eg[1] / eg[25];
    eqk[64] = eg[12] * eg[14] / eg[1] / eg[25];
    eqk[65] = eg[17] / eg[0] / eg[14] / pfac;
    eqk[66] = eg[1] * eg[5] / eg[0] / eg[4];
    eqk[67] = eg[7] / eg[4] / eg[4] / pfac;
    eqk[68] = eg[2] * eg[5] / eg[4] / eg[4];
    eqk[69] = eg[3] * eg[5] / eg[4] / eg[6];
    eqk[147] = eqk[69];
    eqk[70] = eg[5] * eg[6] / eg[4] / eg[7];
    eqk[71] = eqk[70];
    eqk[72] = eg[1] * eg[14] / eg[4] / eg[8];
    eqk[73] = eg[1] * eg[16] / eg[4] / eg[9];
    eqk[74] = eg[1] * eg[17] / eg[4] / eg[10];
    eqk[75] = eg[5] * eg[9] / eg[4] / eg[10];
    eqk[76] = eg[1] * eg[17] / eg[4] / eg[11];
    eqk[77] = eg[5] * eg[10] / eg[4] / eg[12];
    eqk[78] = eg[5] * eg[11] / eg[4] / eg[12];
    eqk[79] = eg[5] * eg[12] / eg[4] / eg[13];
    eqk[80] = eg[1] * eg[15] / eg[4] / eg[14];
    eqk[81] = eg[5] * eg[14] / eg[4] / eg[16];
    eqk[82] = eg[5] * eg[16] / eg[4] / eg[17];
    eqk[83] = eg[5] * eg[17] / eg[4] / eg[18];
    eqk[84] = eg[1] * eg[25] / eg[4] / eg[19];
    eqk[85] = eg[12] * eg[14] / eg[4] / eg[19];
    eqk[86] = eg[5] * eg[19] / eg[4] / eg[20];
    eqk[87] = eg[5] * eg[20] / eg[4] / eg[21];
    eqk[88] = eg[5] * eg[22] / eg[4] / eg[23];
    eqk[89] = eg[5] * eg[24] / eg[4] / eg[25];
    eqk[90] = eg[3] * eg[7] / eg[6] / eg[6];
    eqk[91] = eqk[90];
    eqk[92] = eg[4] * eg[17] / eg[6] / eg[10];
    eqk[93] = eg[3] * eg[13] / eg[6] / eg[12];
    eqk[94] = eg[4] * eg[18] / eg[6] / eg[12];
    eqk[95] = eg[4] * eg[15] / eg[6] / eg[14];
    eqk[96] = eg[7] * eg[16] / eg[6] / eg[17];
    eqk[97] = eg[2] * eg[14] / eg[3] / eg[8];
    eqk[98] = eg[1] * eg[19] / eg[8] / eg[12];
    eqk[99] = eg[2] * eg[16] / eg[3] / eg[9];
    eqk[100] = eg[1] * eg[10] / eg[0] / eg[9];
    eqk[101] = eg[1] * eg[17] / eg[5] / eg[9];
    eqk[102] = eg[1] * eg[19] / eg[9] / eg[10];
    eqk[104] = eg[1] * eg[21] / eg[9] / eg[13];
    eqk[105] = eg[24] / eg[9] / eg[14] / pfac;
    eqk[106] = eg[14] * eg[16] / eg[9] / eg[15];
    eqk[107] = eg[1] * eg[25] / eg[9] / eg[17];
    eqk[108] = eg[14] * eg[19] / eg[9] / eg[24];
    eqk[110] = eg[1] * eg[12] / eg[0] / eg[10];
    eqk[111] = eg[0] * eg[19] / eg[10] / eg[10];
    eqk[112] = eg[1] * eg[21] / eg[10] / eg[12];
    eqk[113] = eg[12] * eg[12] / eg[10] / eg[13];
    eqk[114] = eg[25] / eg[10] / eg[14] / pfac;
    eqk[116] = eg[10] / eg[11];
    eqk[120] = eqk[116];
    eqk[123] = eqk[116];
    eqk[124] = eqk[116];
    eqk[117] = eg[1] * eg[4] * eg[14] / eg[3] / eg[11] * pfac;
    eqk[118] = eg[5] * eg[14] / eg[3] / eg[11];
    eqk[119] = eg[1] * eg[12] / eg[0] / eg[11];
    eqk[121] = eg[1] * eg[21] / eg[11] / eg[12];
    eqk[122] = eg[12] * eg[12] / eg[11] / eg[13];
    eqk[125] = eg[14] * eg[17] / eg[11] / eg[15];
    eqk[127] = eg[2] * eg[18] / eg[3] / eg[12];
    eqk[128] = eg[4] * eg[17] / eg[3] / eg[12];
    eqk[129] = eg[6] * eg[13] / eg[7] / eg[12];
    eqk[130] = eg[23] / eg[12] / eg[12] / pfac;
    eqk[131] = eg[1] * eg[22] / eg[12] / eg[12];
    eqk[132] = eg[13] * eg[14] / eg[12] / eg[16];
    eqk[133] = eg[13] * eg[16] / eg[12] / eg[17];
    eqk[134] = eg[13] * eg[20] / eg[12] / eg[21];
    eqk[135] = eg[13] * eg[22] / eg[12] / eg[23];
    eqk[136] = eg[1] * eg[14] / eg[16] * pfac;
    eqk[137] = eqk[136];
    eqk[138] = eg[6] * eg[14] / eg[3] / eg[16];
    eqk[139] = eg[6] * eg[17] / eg[3] / eg[18];
    eqk[140] = eg[16] * eg[17] / eg[3] / eg[20];
    eqk[141] = eg[0] * eg[19] / eg[21] * pfac;
    eqk[142] = eg[6] * eg[21] / eg[3] / eg[22];
    eqk[143] = eg[4] * eg[14] * eg[14] / eg[3] / eg[24] * pfac;
    eqk[144] = eg[14] * eg[14] * eg[19] / eg[24] / eg[24] * pfac;
    eqk[146] = eg[1] * eg[26] / eg[2] / eg[21];
    eqk[149] = eg[12] / eg[0] / eg[9] / pfac;
    eqk[151] = eg[2] * eg[17] / eg[3] / eg[10];
    eqk[155] = eg[6] * eg[19] / eg[3] / eg[20];
    eqk[156] = eg[26] / eg[1] / eg[25] / pfac;
    eqk[160] = eg[12] * eg[16] / eg[1] / eg[26];
    eqk[161] = eg[0] * eg[25] / eg[1] / eg[26];
    eqk[162] = eg[5] * eg[25] / eg[4] / eg[26];
    eqk[163] = eg[12] * eg[21] / eg[1] / eg[27];
    eqk[164] = eg[1] * eg[12] * eg[25] / eg[2] / eg[27] * pfac;
    eqk[166] = eg[27] / eg[12] / eg[20] / pfac;

    rb[1] = rf[1] / max(eqk[0],small);
    rb[2] = rf[2] / max(eqk[1],small);
    rb[3] = rf[3] / max(eqk[2],small);
    rb[4] = rf[4] / max(eqk[3],small);
    rb[5] = rf[5] / max(eqk[4],small);
    rb[6] = rf[6] / max(eqk[5],small);
    rb[7] = rf[7] / max(eqk[6],small);
    rb[8] = rf[8] / max(eqk[7],small);
    rb[9] = rf[9] / max(eqk[8],small);
    rb[10] = rf[10] / max(eqk[9],small);
    rb[11] = rf[11] / max(eqk[10],small);
    rb[12] = rf[12] / max(eqk[11],small);
    rb[13] = rf[13] / max(eqk[12],small);
    rb[14] = rf[14] / max(eqk[13],small);
    rb[15] = rf[15] / max(eqk[14],small);
    rb[16] = rf[16] / max(eqk[15],small);
    rb[17] = rf[17] / max(eqk[16],small);
    rb[18] = rf[18] / max(eqk[17],small);
    rb[19] = rf[19] / max(eqk[18],small);
    rb[20] = rf[20] / max(eqk[19],small);
    rb[21] = rf[21] / max(eqk[20],small);
    rb[22] = rf[22] / max(eqk[21],small);
    rb[23] = rf[23] / max(eqk[22],small);
    rb[24] = rf[24] / max(eqk[23],small);
    rb[25] = rf[25] / max(eqk[24],small);
    rb[26] = rf[26] / max(eqk[25],small);
    rb[27] = rf[27] / max(eqk[26],small);
    rb[28] = rf[28] / max(eqk[27],small);
    rb[29] = rf[29] / max(eqk[28],small);
    rb[30] = rf[30] / max(eqk[29],small);
    rb[31] = rf[31] / max(eqk[30],small);
    rb[32] = rf[32] / max(eqk[31],small);
    rb[33] = rf[33] / max(eqk[32],small);
    rb[34] = rf[34] / max(eqk[33],small);
    rb[35] = rf[35] / max(eqk[34],small);
    rb[36] = rf[36] / max(eqk[35],small);
    rb[37] = rf[37] / max(eqk[36],small);
    rb[38] = rf[38] / max(eqk[37],small);
    rb[39] = rf[39] / max(eqk[38],small);
    rb[40] = rf[40] / max(eqk[39],small);
    rb[41] = rf[41] / max(eqk[40],small);
    rb[42] = rf[42] / max(eqk[41],small);
    rb[43] = rf[43] / max(eqk[42],small);
    rb[44] = rf[44] / max(eqk[43],small);
    rb[45] = rf[45] / max(eqk[44],small);
    rb[46] = rf[46] / max(eqk[45],small);
    rb[47] = rf[47] / max(eqk[46],small);
    rb[48] = rf[48] / max(eqk[47],small);
    rb[49] = rf[49] / max(eqk[48],small);
    rb[50] = rf[50] / max(eqk[49],small);
    rb[51] = rf[51] / max(eqk[50],small);
    rb[52] = rf[52] / max(eqk[51],small);
    rb[53] = rf[53] / max(eqk[52],small);
    rb[54] = rf[54] / max(eqk[53],small);
    rb[55] = rf[55] / max(eqk[54],small);
    rb[56] = rf[56] / max(eqk[55],small);
    rb[57] = rf[57] / max(eqk[56],small);
    rb[58] = rf[58] / max(eqk[57],small);
    rb[59] = rf[59] / max(eqk[58],small);
    rb[60] = rf[60] / max(eqk[59],small);
    rb[61] = rf[61] / max(eqk[60],small);
    rb[62] = rf[62] / max(eqk[61],small);
    rb[63] = rf[63] / max(eqk[62],small);
    rb[64] = rf[64] / max(eqk[63],small);
    rb[65] = rf[65] / max(eqk[64],small);
    rb[66] = rf[66] / max(eqk[65],small);
    rb[67] = rf[67] / max(eqk[66],small);
    rb[68] = rf[68] / max(eqk[67],small);
    rb[69] = rf[69] / max(eqk[68],small);
    rb[70] = rf[70] / max(eqk[69],small);
    rb[71] = rf[71] / max(eqk[70],small);
    rb[72] = rf[72] / max(eqk[71],small);
    rb[73] = rf[73] / max(eqk[72],small);
    rb[74] = rf[74] / max(eqk[73],small);
    rb[75] = rf[75] / max(eqk[74],small);
    rb[76] = rf[76] / max(eqk[75],small);
    rb[77] = rf[77] / max(eqk[76],small);
    rb[78] = rf[78] / max(eqk[77],small);
    rb[79] = rf[79] / max(eqk[78],small);
    rb[80] = rf[80] / max(eqk[79],small);
    rb[81] = rf[81] / max(eqk[80],small);
    rb[82] = rf[82] / max(eqk[81],small);
    rb[83] = rf[83] / max(eqk[82],small);
    rb[84] = rf[84] / max(eqk[83],small);
    rb[85] = rf[85] / max(eqk[84],small);
    rb[86] = rf[86] / max(eqk[85],small);
    rb[87] = rf[87] / max(eqk[86],small);
    rb[88] = rf[88] / max(eqk[87],small);
    rb[89] = rf[89] / max(eqk[88],small);
    rb[90] = rf[90] / max(eqk[89],small);
    rb[91] = rf[91] / max(eqk[90],small);
    rb[92] = rf[92] / max(eqk[91],small);
    rb[93] = rf[93] / max(eqk[92],small);
    rb[94] = rf[94] / max(eqk[93],small);
    rb[95] = rf[95] / max(eqk[94],small);
    rb[96] = rf[96] / max(eqk[95],small);
    rb[97] = rf[97] / max(eqk[96],small);
    rb[98] = rf[98] / max(eqk[97],small);
    rb[99] = rf[99] / max(eqk[98],small);
    rb[100] = rf[100] / max(eqk[99],small);
    rb[101] = rf[101] / max(eqk[100],small);
    rb[102] = rf[102] / max(eqk[101],small);
    rb[103] = rf[103] / max(eqk[102],small);
    rf[103] = 0.;
    rb[105] = rf[105] / max(eqk[104],small);
    rb[106] = rf[106] / max(eqk[105],small);
    rb[107] = rf[107] / max(eqk[106],small);
    rb[108] = rf[108] / max(eqk[107],small);
    rb[109] = rf[109] / max(eqk[108],small);
    rf[109] = 0.;
    rb[110] = 0.;
    rb[111] = rf[111] / max(eqk[110],small);
    rb[112] = rf[112] / max(eqk[111],small);
    rf[112] = 0.;
    rb[113] = rf[113] / max(eqk[112],small);
    rb[114] = rf[114] / max(eqk[113],small);
    rb[115] = rf[115] / max(eqk[114],small);
    rb[117] = rf[117] / max(eqk[116],small);
    rb[118] = rf[118] / max(eqk[117],small);
    rb[119] = rf[119] / max(eqk[118],small);
    rb[120] = rf[120] / max(eqk[119],small);
    rb[121] = rf[121] / max(eqk[120],small);
    rb[122] = rf[122] / max(eqk[121],small);
    rb[123] = rf[123] / max(eqk[122],small);
    rb[124] = rf[124] / max(eqk[123],small);
    rb[125] = rf[125] / max(eqk[124],small);
    rb[126] = rf[126] / max(eqk[125],small);
    rb[128] = rf[128] / max(eqk[127],small);
    rb[129] = rf[129] / max(eqk[128],small);
    rb[130] = rf[130] / max(eqk[129],small);
    rb[131] = rf[131] / max(eqk[130],small);
    rb[132] = rf[132] / max(eqk[131],small);
    rb[133] = rf[133] / max(eqk[132],small);
    rb[134] = rf[134] / max(eqk[133],small);
    rb[135] = rf[135] / max(eqk[134],small);
    rb[136] = rf[136] / max(eqk[135],small);
    rb[137] = rf[137] / max(eqk[136],small);
    rb[138] = rf[138] / max(eqk[137],small);
    rb[139] = rf[139] / max(eqk[138],small);
    rb[140] = rf[140] / max(eqk[139],small);
    rb[142] = rf[142] / max(eqk[141],small);
    rb[143] = rf[143] / max(eqk[142],small);
    rb[144] = rf[144] / max(eqk[143],small);
    rb[145] = rf[145] / max(eqk[144],small);
    rf[145] = 0.;
    rb[146] = 0.f;
    rb[147] = rf[147] / max(eqk[146],small);
    rb[148] = rf[148] / max(eqk[147],small);
    rb[149] = 0.f;
    rb[150] = rf[150] / max(eqk[149],small);
    rb[151] = 0.f;
    rb[152] = rf[152] / max(eqk[151],small);
    rb[153] = 0.f;
    rb[154] = 0.f;
    rb[156] = rf[156] / max(eqk[155],small);
    rb[157] = rf[157] / max(eqk[156],small);
    rb[159] = 0.f;
    rb[161] = rf[161] / max(eqk[160],small);
    rf[161] = 0.;
    rb[162] = rf[162] / max(eqk[161],small);
    rb[163] = rf[163] / max(eqk[162],small);
    rb[164] = rf[164] / max(eqk[163],small);
    rb[165] = rf[165] / max(eqk[164],small);
    rb[167] = rf[167] / max(eqk[166],small);

    rklow[1] = exp(34.0312786 - 1509.65 / *t);
    rklow[2] = exp(59.9064331 - *alogt * 2.76 - 805.146668 / *t);
    rklow[3] = exp(76.9484824 - *alogt * 4.76 - 1227.84867 / *t);
    rklow[4] = exp(56.1662604 - *alogt * 2.57 - 213.867084 / *t);
    rklow[5] = exp(69.8660102 - *alogt * 4.8 - 2797.88467 / *t);
    rklow[6] = exp(93.4384048 - *alogt * 7.27 - 3633.22434 / *t);
    rklow[7] = exp(69.414025 - *alogt * 3.86 - 1670.67934 / *t);
    rklow[8] = exp(96.1977483 - *alogt * 7.62 - 3507.42017 / *t);
    rklow[9] = exp(95.0941235 - *alogt * 7.08 - 3364.00342 / *t);
    rklow[10] = exp(63.7931383 - *alogt * 3.42 - 42446.3259 / *t);
    rklow[11] = exp(42.2794408 - *alogt * .9 + 855.468335 / *t);
    rklow[12] = exp(65.4619238 - *alogt * 3.74 - 974.227469 / *t);
    rklow[13] = exp(76.9748493 - *alogt * 5.11 - 3570.32226 / *t);
    rklow[14] = exp(95.6297642 - *alogt * 7.03 - 1389.88444 / *t);
    rklow[15] = exp(117.889265 - *alogt * 9.3 - 49214.5901 / *t);
    rklow[16] = exp(59.1374013 - *alogt * 2.8 - 296.897834 / *t);
    rklow[17] = exp(96.7205025 - *alogt * 7.63 - 1939.39704 / *t);
    rklow[18] = exp(135.001549 - *alogt * 11.94 - 4916.3262 / *t);

    return 0;
} /* ratt_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int rdsmh_(double *t, double *tlog, double *smh)
{
    static double ti, tn[5];


/*  START PROLOGUE */

/*  SUBROUTINE DKSMH  (T, ICKWRK, RCKWRK, SMH)* */
/*     Returns the array of entropies minus enthalpies for the species. */
/*     It is normally not called directly by the user. */

/*  INPUT */
/*     T      - Temperature. */
/*                   cgs units - K */
/*                   Data type - real scalar */
/*     TLOG      - Log Temperature. */
/*  OUTPUT */
/*     SMH    - Entropy minus enthalpy for the species, */
/*              SMH(K) = S(K)/R - H(K)/RT. */
/*                   cgs units - none */
/*                   Data type - real array */
/*                   Dimension SMH(*) at least KK, the total number of */
/*                   species. */

/*  END PROLOGUE */

/* *****precision > double */
/* *****END precision > double */
/* *****precision > single */
/*        IMPLICIT REAL (A-H, O-Z), int (I-N) */
/* *****END precision > single */


    /* Parameter adjustments */
    --smh;

    /* Function Body */
    ti = 1. / *t;

    tn[0] = *tlog - 1.f;
    tn[1] = *t;
    tn[2] = tn[1] * *t;
    tn[3] = tn[2] * *t;
    tn[4] = tn[3] * *t;

    if (*t > 1e3) {

    smh[1] = ti * 950.158922 - 3.20502331 + tn[0] * 3.3372792 - tn[1] *
        2.47012365e-5 + tn[2] * 8.32427963e-8 - tn[3] *
        1.49638662e-11 + tn[4] * 1.00127688e-15;
    smh[2] = -.446682914 - ti * 25473.6599 + tn[0] * 2.50000001 - tn[1] *
        1.15421486e-11 + tn[2] * 2.69269913e-15 - tn[3] *
        3.94596029e-19 + tn[4] * 2.49098679e-23;
    smh[3] = 4.78433864 - ti * 29217.5791 + tn[0] * 2.56942078 - tn[1] *
        4.29870569e-5 + tn[2] * 6.99140982e-9 - tn[3] *
        8.34814992e-13 + tn[4] * 6.14168455e-17;
    smh[4] = ti * 1088.45772 + 5.45323129 + tn[0] * 3.28253784 + tn[1] *
        7.4154377e-4 - tn[2] * 1.26327778e-7 + tn[3] * 1.74558796e-11
        - tn[4] * 1.08358897e-15;
    smh[5] = 4.4766961 - ti * 3858.657 + tn[0] * 3.09288767 + tn[1] *
        2.74214858e-4 + tn[2] * 2.10842047e-8 - tn[3] * 7.3288463e-12
        + tn[4] * 5.8706188e-16;
    smh[6] = ti * 30004.2971 + 4.9667701 + tn[0] * 3.03399249 + tn[1] *
        .00108845902 - tn[2] * 2.73454197e-8 - tn[3] * 8.08683225e-12
        + tn[4] * 8.4100496e-16;
    smh[7] = 3.78510215 - ti * 111.856713 + tn[0] * 4.0172109 + tn[1] *
        .00111991007 - tn[2] * 1.05609692e-7 + tn[3] * 9.52053083e-12
        - tn[4] * 5.39542675e-16;
    smh[8] = ti * 17861.7877 + 2.91615662 + tn[0] * 4.16500285 + tn[1] *
        .00245415847 - tn[2] * 3.16898708e-7 + tn[3] * 3.09321655e-11
        - tn[4] * 1.43954153e-15;
    smh[9] = 4.80150373 - ti * 85451.2953 + tn[0] * 2.49266888 + tn[1] *
        2.39944642e-5 - tn[2] * 1.20722503e-8 + tn[3] *
        3.11909191e-12 - tn[4] * 2.43638946e-16;
    smh[10] = 5.48497999 - ti * 71012.4364 + tn[0] * 2.87846473 + tn[1] *
        4.85456841e-4 + tn[2] * 2.40742758e-8 - tn[3] *
        1.08906541e-11 + tn[4] * 8.80396915e-16;
    smh[11] = 6.17119324 - ti * 46263.604 + tn[0] * 2.87410113 + tn[1] *
        .00182819646 - tn[2] * 2.34824328e-7 + tn[3] * 2.16816291e-11
        - tn[4] * 9.38637835e-16;
    smh[12] = 8.62650169 - ti * 50925.9997 + tn[0] * 2.29203842 + tn[1] *
        .00232794319 - tn[2] * 3.35319912e-7 + tn[3] * 3.48255e-11 -
        tn[4] * 1.69858183e-15;
    smh[13] = 8.48007179 - ti * 16775.5843 + tn[0] * 2.28571772 + tn[1] *
        .00361995018 - tn[2] * 4.97857247e-7 + tn[3] * 4.9640387e-11
        - tn[4] * 2.33577197e-15;
    smh[14] = ti * 9468.34459 + 18.437318 + tn[0] * .074851495 + tn[1] *
        .00669547335 - tn[2] * 9.55476348e-7 + tn[3] * 1.01910446e-10
        - tn[4] * 5.0907615e-15;
    smh[15] = ti * 14151.8724 + 7.81868772 + tn[0] * 2.71518561 + tn[1] *
        .00103126372 - tn[2] * 1.66470962e-7 + tn[3] * 1.9171084e-11
        - tn[4] * 1.01823858e-15;
    smh[16] = ti * 48759.166 + 2.27163806 + tn[0] * 3.85746029 + tn[1] *
        .00220718513 - tn[2] * 3.69135673e-7 + tn[3] * 4.36241823e-11
        - tn[4] * 2.36042082e-15;
    smh[17] = 9.79834492 - ti * 4011.91815 + tn[0] * 2.77217438 + tn[1] *
        .00247847763 - tn[2] * 4.14076022e-7 + tn[3] * 4.90968148e-11
        - tn[4] * 2.66754356e-15;
    smh[18] = ti * 13995.8323 + 13.656323 + tn[0] * 1.76069008 + tn[1] *
        .00460000041 - tn[2] * 7.37098022e-7 + tn[3] * 8.38676767e-11
        - tn[4] * 4.4192782e-15;
    smh[19] = 2.929575 - ti * 127.83252 + tn[0] * 3.770799 + tn[1] *
        .0039357485 - tn[2] * 4.42730667e-7 + tn[3] * 3.28702583e-11
        - tn[4] * 1.056308e-15;
    smh[20] = -1.23028121 - ti * 25935.9992 + tn[0] * 4.14756964 + tn[1] *
         .00298083332 - tn[2] * 3.9549142e-7 + tn[3] * 3.89510143e-11
        - tn[4] * 1.80617607e-15;
    smh[21] = 7.78732378 - ti * 34612.8739 + tn[0] * 3.016724 + tn[1] *
        .0051651146 - tn[2] * 7.80137248e-7 + tn[3] * 8.480274e-11 -
        tn[4] * 4.3130352e-15;
    smh[22] = 10.3053693 - ti * 4939.88614 + tn[0] * 2.03611116 + tn[1] *
        .00732270755 - tn[2] * 1.11846319e-6 + tn[3] * 1.22685769e-10
        - tn[4] * 6.28530305e-15;
    smh[23] = 13.4624343 - ti * 12857.52 + tn[0] * 1.95465642 + tn[1] *
        .0086986361 - tn[2] * 1.33034445e-6 + tn[3] * 1.46014741e-10
        - tn[4] * 7.4820788e-15;
    smh[24] = ti * 11426.3932 + 15.1156107 + tn[0] * 1.0718815 + tn[1] *
        .0108426339 - tn[2] * 1.67093445e-6 + tn[3] * 1.84510001e-10
        - tn[4] * 9.5001445e-15;
    smh[25] = -3.9302595 - ti * 19327.215 + tn[0] * 5.6282058 + tn[1] *
        .00204267005 - tn[2] * 2.65575783e-7 + tn[3] * 2.38550433e-11
        - tn[4] * 9.703916e-16;
    smh[26] = ti * 7551.05311 + .632247205 + tn[0] * 4.51129732 + tn[1] *
        .00450179872 - tn[2] * 6.94899392e-7 + tn[3] * 7.69454902e-11
        - tn[4] * 3.974191e-15;
    smh[27] = -5.0320879 - ti * 490.32178 + tn[0] * 5.9756699 + tn[1] *
        .0040652957 - tn[2] * 4.5727075e-7 + tn[3] * 3.39192008e-11 -
        tn[4] * 1.08800855e-15;
    smh[28] = ti * 923.5703 - 13.31335 + tn[0] * 6.732257 + tn[1] *
        .00745417 - tn[2] * 8.24983167e-7 + tn[3] * 6.01001833e-11 -
        tn[4] * 1.883102e-15;

    } else {

    smh[1] = ti * 917.935173 + .683010238 + tn[0] * 2.34433112 + tn[1] *
        .00399026037 - tn[2] * 3.2463585e-6 + tn[3] * 1.67976745e-9 -
        tn[4] * 3.68805881e-13;
    smh[2] = -.446682853 - ti * 25473.6599 + tn[0] * 2.5 + tn[1] *
        3.52666409e-13 - tn[2] * 3.32653273e-16 + tn[3] *
        1.91734693e-19 - tn[4] * 4.63866166e-23;
    smh[3] = 2.05193346 - ti * 29122.2592 + tn[0] * 3.1682671 - tn[1] *
        .00163965942 + tn[2] * 1.10717733e-6 - tn[3] * 5.10672187e-10
        + tn[4] * 1.05632986e-13;
    smh[4] = ti * 1063.94356 + 3.65767573 + tn[0] * 3.78245636 - tn[1] *
        .00149836708 + tn[2] * 1.641217e-6 - tn[3] * 8.06774591e-10 +
        tn[4] * 1.62186419e-13;
    smh[5] = -.103925458 - ti * 3615.08056 + tn[0] * 3.99201543 - tn[1] *
        .00120065876 + tn[2] * 7.69656402e-7 - tn[3] * 3.23427778e-10
        + tn[4] * 6.8205735e-14;
    smh[6] = ti * 30293.7267 - .849032208 + tn[0] * 4.19864056 - tn[1] *
        .00101821705 + tn[2] * 1.08673369e-6 - tn[3] * 4.57330885e-10
        + tn[4] * 8.85989085e-14;
    smh[7] = 3.71666245 - ti * 294.80804 + tn[0] * 4.30179801 - tn[1] *
        .00237456025 + tn[2] * 3.52638152e-6 - tn[3] * 2.02303245e-9
        + tn[4] * 4.64612562e-13;
    smh[8] = ti * 17702.5821 + 3.43505074 + tn[0] * 4.27611269 - tn[1] *
        2.71411208e-4 + tn[2] * 2.78892835e-6 - tn[3] * 1.79809011e-9
        + tn[4] * 4.31227182e-13;
    smh[9] = 4.53130848 - ti * 85443.8832 + tn[0] * 2.55423955 - tn[1] *
        1.60768862e-4 + tn[2] * 1.22298708e-7 - tn[3] *
        6.10195741e-11 + tn[4] * 1.33260723e-14;
    smh[10] = 2.08401108 - ti * 70797.2934 + tn[0] * 3.48981665 + tn[1] *
        1.61917771e-4 - tn[2] * 2.81498442e-7 + tn[3] *
        2.63514439e-10 - tn[4] * 7.03045335e-14;
    smh[11] = 1.56253185 - ti * 46004.0401 + tn[0] * 3.76267867 + tn[1] *
        4.84436072e-4 + tn[2] * 4.65816402e-7 - tn[3] *
        3.20909294e-10 + tn[4] * 8.43708595e-14;
    smh[12] = -.769118967 - ti * 50496.8163 + tn[0] * 4.19860411 - tn[1] *
         .0011833071 + tn[2] * 1.37216037e-6 - tn[3] * 5.57346651e-10
        + tn[4] * 9.71573685e-14;
    smh[13] = 1.60456433 - ti * 16444.9988 + tn[0] * 3.6735904 + tn[1] *
        .00100547588 + tn[2] * 9.55036427e-7 - tn[3] * 5.72597854e-10
        + tn[4] * 1.27192867e-13;
    smh[14] = ti * 10246.6476 - 4.64130376 + tn[0] * 5.14987613 - tn[1] *
        .0068354894 + tn[2] * 8.19667665e-6 - tn[3] * 4.03952522e-9 +
        tn[4] * 8.3346978e-13;
    smh[15] = ti * 14344.086 + 3.50840928 + tn[0] * 3.57953347 - tn[1] *
        3.0517684e-4 + tn[2] * 1.69469055e-7 + tn[3] * 7.55838237e-11
        - tn[4] * 4.52212249e-14;
    smh[16] = ti * 48371.9697 + 9.90105222 + tn[0] * 2.35677352 + tn[1] *
        .00449229839 - tn[2] * 1.18726045e-6 + tn[3] * 2.04932518e-10
        - tn[4] * 7.1849774e-15;
    smh[17] = 3.39437243 - ti * 3839.56496 + tn[0] * 4.22118584 - tn[1] *
        .00162196266 + tn[2] * 2.29665743e-6 - tn[3] * 1.10953411e-9
        + tn[4] * 2.16884433e-13;
    smh[18] = ti * 14308.9567 + .6028129 + tn[0] * 4.79372315 - tn[1] *
        .00495416685 + tn[2] * 6.22033347e-6 - tn[3] * 3.16071051e-9
        + tn[4] * 6.5886326e-13;
    smh[19] = 13.152177 - ti * 978.6011 + tn[0] * 2.106204 + tn[1] *
        .0036082975 + tn[2] * 8.89745333e-7 - tn[3] * 6.14803e-10 +
        tn[4] * 1.037805e-13;
    smh[20] = 13.9397051 - ti * 26428.9807 + tn[0] * .808681094 + tn[1] *
        .0116807815 - tn[2] * 5.91953025e-6 + tn[3] * 2.33460364e-9 -
        tn[4] * 4.25036487e-13;
    smh[21] = 8.51054025 - ti * 34859.8468 + tn[0] * 3.21246645 + tn[1] *
        7.5739581e-4 + tn[2] * 4.32015687e-6 - tn[3] * 2.98048206e-9
        + tn[4] * 7.35754365e-13;
    smh[22] = 4.09733096 - ti * 5089.77593 + tn[0] * 3.95920148 - tn[1] *
        .00378526124 + tn[2] * 9.51650487e-6 - tn[3] * 5.76323961e-9
        + tn[4] * 1.34942187e-12;
    smh[23] = 4.70720924 - ti * 12841.6265 + tn[0] * 4.30646568 - tn[1] *
        .00209329446 + tn[2] * 8.28571345e-6 - tn[3] * 4.99272172e-9
        + tn[4] * 1.15254502e-12;
    smh[24] = ti * 11522.2055 + 2.66682316 + tn[0] * 4.29142492 - tn[1] *
        .00275077135 + tn[2] * 9.99063813e-6 - tn[3] * 5.90388571e-9
        + tn[4] * 1.34342886e-12;
    smh[25] = 12.490417 - ti * 20059.449 + tn[0] * 2.2517214 + tn[1] *
        .0088275105 - tn[2] * 3.95485017e-6 + tn[3] * 1.43964658e-9 -
        tn[4] * 2.53324055e-13;
    smh[26] = ti * 7042.91804 + 12.215648 + tn[0] * 2.1358363 + tn[1] *
        .00905943605 - tn[2] * 2.89912457e-6 + tn[3] * 7.7866464e-10
        - tn[4] * 1.00728807e-13;
    smh[27] = 9.5714535 - ti * 1521.4766 + tn[0] * 3.4090624 + tn[1] *
        .005369287 + tn[2] * 3.1524875e-7 + tn[3] * 5.96548592e-10 +
        tn[4] * 1.43369255e-13;
    smh[28] = 16.14534 - ti * 1074.826 + tn[0] * 1.493307 + tn[1] *
        .01046259 + tn[2] * 7.47799e-7 - tn[3] * 1.39076e-9 + tn[4] *
        3.579073e-13;
    }

    return 0;
} /* rdsmh_ */

/*                                                                      C */
/* ----------------------------------------------------------------------C */
/*                                                                      C */
/* Subroutine */ int ratx_(double *t, double *alogt, double *c__,
    double *rf, double *rb, double *rklow)
{
    /* Initialized data */

    static double small = 1e-50;

    /* System generated locals */
    double d__1;

    /* Builtin functions */
    //dol double log10(double *), exp(double), pow(double *,
     //dol    double *);

    /* Local variables */
    static int k;
    static double fc, pr, xn, ctb[167], flog, pcor, ctot, fclog, fcent,
        prlog, cprlog;


/* *****precision > double */
/* *****END precision > double */

    /* Parameter adjustments */
    --rklow;
    --rb;
    --rf;
    --c__;

    /* Function Body */

/*     third-body reactions */

    ctot = 0.f;
    for (k = 1; k <= 19; ++k) {
    ctot += c__[k];
    }

    ctb[0] = ctot + c__[1] * 1.4 + c__[6] * 14.4 + c__[10] + c__[11] * .75 +
        c__[12] * 2.6 + c__[16] * 2. + c__[15] * 2. + c__[18] * 3.;
    ctb[1] = ctot + c__[1] + c__[6] * 5. + c__[10] + c__[11] * .5 + c__[12] +
        c__[16] * 2. + c__[15] * 2. + c__[18] * 3.;
    ctb[11] = ctot + c__[1] + c__[4] * 5. + c__[6] * 5. + c__[10] + c__[11] *
        .5 + c__[12] * 2.5 + c__[16] * 2. + c__[15] * 2. + c__[18] * 3.;
    ctb[27] = ctot - c__[4] - c__[6] - c__[11] * .25 + c__[12] * .5 + c__[16]
        * .5 - c__[19] + c__[15] * 2. + c__[18] * 3.;
    ctb[32] = ctot - c__[1] - c__[6] + c__[10] - c__[12] + c__[16] * 2. + c__[
        15] * 2. + c__[18] * 3.;
    ctb[36] = ctot - c__[1] * .27 + c__[6] * 2.65 + c__[10] + c__[16] * 2. +
        c__[15] * 2. + c__[18] * 3.;
    ctb[43] = ctb[1];
    ctb[45] = ctot + c__[1] + c__[6] * 5. + c__[10] * 2. + c__[11] * .5 + c__[
        12] + c__[16] * 2. + c__[15] * 2. + c__[18] * 3.;
    ctb[47] = ctb[1];
    ctb[49] = ctb[1];
    ctb[54] = ctb[1];
    ctb[55] = ctb[1];
    ctb[57] = ctb[1];
    ctb[59] = ctb[1];
    ctb[65] = ctb[1];
    ctb[67] = ctb[1];
    ctb[105] = ctb[1];
    ctb[114] = ctb[1];
    ctb[130] = ctb[1];
    ctb[137] = ctot + c__[1] - c__[6] + c__[10] + c__[11] * .5 + c__[12] +
        c__[16] * 2.;
    ctb[141] = ctb[1];
    ctb[149] = ctb[1];
    ctb[156] = ctb[1];
    ctb[166] = ctb[1];

/*     If fall-off (pressure correction): */


    pr = rklow[1] * ctb[11] / rf[12];
    pcor = pr / (pr + 1.f);
    rf[12] *= pcor;
    rb[12] *= pcor;

    pr = rklow[2] * ctb[43] / rf[44];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 91.) * .438 + exp(-(*t) / 5836.) * .562 + exp(-8552. /
         *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[44] *= pcor;
    rb[44] *= pcor;

    pr = rklow[3] * ctb[45] / rf[46];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 74.) * .217 + exp(-(*t) / 2941.) * .783 + exp(-6964. /
         *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[46] *= pcor;
    rb[46] *= pcor;

    pr = rklow[4] * ctb[47] / rf[48];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 271.) * .2176 + exp(-(*t) / 2755.) * .7824 + exp(
        -6570. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[48] *= pcor;
    rb[48] *= pcor;

    pr = rklow[5] * ctb[49] / rf[50];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 94.) * .242 + exp(-(*t) / 1555.) * .758 + exp(-4200. /
         *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[50] *= pcor;
    rb[50] *= pcor;

    pr = rklow[6] * ctb[54] / rf[55];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 98.5) * .2493 + exp(-(*t) / 1302.) * .7507 + exp(
        -4167. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[55] *= pcor;
    rb[55] *= pcor;

    pr = rklow[7] * ctb[55] / rf[56];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 207.5) * .218 + exp(-(*t) / 2663.) * .782 + exp(
        -6095. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[56] *= pcor;
    rb[56] *= pcor;

    pr = rklow[8] * ctb[57] / rf[58];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 210.) * .0247 + exp(-(*t) / 984.) * .9753 + exp(
        -4374. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[58] *= pcor;
    rb[58] *= pcor;

    pr = rklow[9] * ctb[59] / rf[60];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 125.) * .1578 + exp(-(*t) / 2219.) * .8422 + exp(
        -6882. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[60] *= pcor;
    rb[60] *= pcor;

    pr = rklow[10] * ctb[65] / rf[66];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 197.) * .068 + exp(-(*t) / 1540.) * .932 + exp(
        -10300. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[66] *= pcor;
    rb[66] *= pcor;

    pr = rklow[11] * ctb[67] / rf[68];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 94.) * .2654 + exp(-(*t) / 1756.) * .7346 + exp(
        -5182. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[68] *= pcor;
    rb[68] *= pcor;

    pr = rklow[12] * ctb[105] / rf[106];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 237.) * .4243 + exp(-(*t) / 1652.) * .5757 + exp(
        -5069. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[106] *= pcor;
    rb[106] *= pcor;

    pr = rklow[13] * ctb[114] / rf[115];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 275.) * .4093 + exp(-(*t) / 1226.) * .5907 + exp(
        -5185. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[115] *= pcor;
    rb[115] *= pcor;

    pr = rklow[14] * ctb[130] / rf[131];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 73.2) * .381 + exp(-(*t) / 1180.) * .619 + exp(-9999.
        / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[131] *= pcor;
    rb[131] *= pcor;

    pr = rklow[15] * ctb[141] / rf[142];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 180.) * .2655 + exp(-(*t) / 1035.) * .7345 + exp(
        -5417. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[142] *= pcor;
    rb[142] *= pcor;

    pr = rklow[16] * ctb[149] / rf[150];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 122.) * .422 + exp(-(*t) / 2535.) * .578 + exp(-9365.
        / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[150] *= pcor;
    rb[150] *= pcor;

    pr = rklow[17] * ctb[156] / rf[157];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 201.) * .535 + exp(-(*t) / 1773.) * .465 + exp(-5333.
        / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[157] *= pcor;
    rb[157] *= pcor;

    pr = rklow[18] * ctb[166] / rf[167];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 1340.6) * .825 + exp(-(*t) / 6e4) * .175 + exp(
        -10139.8 / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b5, flog);
    pcor = fc * pcor;
    rf[167] *= pcor;
    rb[167] *= pcor;

    rf[1] = rf[1] * ctb[0] * c__[3] * c__[3];
    rf[2] = rf[2] * ctb[1] * c__[3] * c__[2];
    rf[3] = rf[3] * c__[3] * c__[1];
    rf[4] = rf[4] * c__[3] * c__[7];
    rf[5] = rf[5] * c__[3] * c__[8];
    rf[6] *= c__[3];
    rf[7] *= c__[3];
    rf[8] *= c__[3];
    rf[9] *= c__[3];
    rf[10] = rf[10] * c__[3] * c__[9];
    rf[11] = rf[11] * c__[3] * c__[10];
    rf[12] = rf[12] * c__[3] * c__[11];
    rf[13] *= c__[3];
    rf[14] *= c__[3];
    rf[15] = rf[15] * c__[3] * c__[13];
    rf[16] *= c__[3];
    rf[17] = rf[17] * c__[3] * c__[14];
    rf[18] = rf[18] * c__[3] * c__[14];
    rf[19] *= c__[3];
    rf[20] = rf[20] * c__[3] * c__[15];
    rf[21] *= c__[3];
    rf[22] = rf[22] * c__[3] * c__[16];
    rf[23] *= c__[3];
    rf[24] = rf[24] * c__[3] * c__[17];
    rf[25] = rf[25] * c__[3] * c__[17];
    rf[26] = rf[26] * c__[4] * c__[11];
    rf[27] = rf[27] * c__[4] * c__[13];
    rf[28] = rf[28] * ctb[27] * c__[2] * c__[4];
    rf[29] = rf[29] * c__[2] * c__[4] * c__[4];
    rf[30] = rf[30] * c__[2] * c__[4] * c__[6];
    rf[31] = rf[31] * c__[2] * c__[4] * c__[19];
    rf[32] = rf[32] * c__[2] * c__[4];
    rf[33] = rf[33] * ctb[32] * c__[2] * c__[2];
    rf[34] = rf[34] * c__[2] * c__[2] * c__[1];
    rf[35] = rf[35] * c__[2] * c__[2] * c__[6];
    rf[36] = rf[36] * c__[2] * c__[2] * c__[12];
    rf[37] = rf[37] * ctb[36] * c__[2] * c__[5];
    rf[38] = rf[38] * c__[2] * c__[7];
    rf[39] = rf[39] * c__[2] * c__[7];
    rf[40] = rf[40] * c__[2] * c__[7];
    rf[41] = rf[41] * c__[2] * c__[8];
    rf[42] = rf[42] * c__[2] * c__[8];
    rf[43] *= c__[2];
    rf[44] *= c__[2];
    rf[45] *= c__[2];
    rf[46] = rf[46] * c__[2] * c__[9];
    rf[47] = rf[47] * c__[2] * c__[10];
    rf[48] *= c__[2];
    rf[49] *= c__[2];
    rf[50] = rf[50] * c__[2] * c__[13];
    rf[51] = rf[51] * c__[2] * c__[13];
    rf[52] *= c__[2];
    rf[53] *= c__[2];
    rf[54] *= c__[2];
    rf[55] = rf[55] * c__[2] * c__[14];
    rf[56] *= c__[2];
    rf[57] *= c__[2];
    rf[58] = rf[58] * c__[2] * c__[15];
    rf[59] = rf[59] * c__[2] * c__[15];
    rf[60] *= c__[2];
    rf[61] *= c__[2];
    rf[62] = rf[62] * c__[2] * c__[16];
    rf[63] *= c__[2];
    rf[64] = rf[64] * c__[2] * c__[17];
    rf[65] = rf[65] * c__[2] * c__[17];
    rf[66] = rf[66] * c__[1] * c__[11];
    rf[67] = rf[67] * c__[5] * c__[1];
    rf[68] = rf[68] * c__[5] * c__[5];
    rf[69] = rf[69] * c__[5] * c__[5];
    rf[70] = rf[70] * c__[5] * c__[7];
    rf[71] = rf[71] * c__[5] * c__[8];
    rf[72] = rf[72] * c__[5] * c__[8];
    rf[73] *= c__[5];
    rf[74] *= c__[5];
    rf[75] *= c__[5];
    rf[76] *= c__[5];
    rf[77] *= c__[5];
    rf[78] = rf[78] * c__[5] * c__[9];
    rf[79] = rf[79] * c__[5] * c__[9];
    rf[80] = rf[80] * c__[5] * c__[10];
    rf[81] = rf[81] * c__[5] * c__[11];
    rf[82] *= c__[5];
    rf[83] = rf[83] * c__[5] * c__[13];
    rf[84] *= c__[5];
    rf[85] = rf[85] * c__[5] * c__[14];
    rf[86] = rf[86] * c__[5] * c__[14];
    rf[87] *= c__[5];
    rf[88] = rf[88] * c__[5] * c__[15];
    rf[89] = rf[89] * c__[5] * c__[16];
    rf[90] = rf[90] * c__[5] * c__[17];
    rf[91] = rf[91] * c__[7] * c__[7];
    rf[92] = rf[92] * c__[7] * c__[7];
    rf[93] *= c__[7];
    rf[94] = rf[94] * c__[7] * c__[9];
    rf[95] = rf[95] * c__[7] * c__[9];
    rf[96] = rf[96] * c__[7] * c__[11];
    rf[97] = rf[97] * c__[7] * c__[13];
    rf[98] *= c__[4];
    rf[99] *= c__[9];
    rf[100] *= c__[4];
    rf[101] *= c__[1];
    rf[102] *= c__[6];
    rf[105] *= c__[10];
    rf[106] *= c__[11];
    rf[107] *= c__[12];
    rf[108] *= c__[13];
    rf[110] *= c__[4];
    rf[111] *= c__[1];
    rf[113] *= c__[9];
    rf[114] *= c__[10];
    rf[115] *= c__[11];
    rf[117] *= c__[19];
    rf[118] *= c__[4];
    rf[119] *= c__[4];
    rf[120] *= c__[1];
    rf[121] *= c__[6];
    rf[122] *= c__[9];
    rf[123] *= c__[10];
    rf[124] *= c__[11];
    rf[125] *= c__[12];
    rf[126] *= c__[12];
    rf[128] = rf[128] * c__[9] * c__[4];
    rf[129] = rf[129] * c__[9] * c__[4];
    rf[130] = rf[130] * c__[9] * c__[8];
    rf[131] = rf[131] * c__[9] * c__[9];
    rf[132] = rf[132] * c__[9] * c__[9];
    rf[133] *= c__[9];
    rf[134] = rf[134] * c__[9] * c__[13];
    rf[135] = rf[135] * c__[9] * c__[15];
    rf[136] = rf[136] * c__[9] * c__[16];
    rf[137] *= c__[6];
    rf[138] *= ctb[137];
    rf[139] *= c__[4];
    rf[140] *= c__[4];
    rf[141] *= c__[4];
    rf[142] *= c__[15];
    rf[143] *= c__[4];
    rf[144] *= c__[4];
    rf[146] = rf[146] * c__[3] * c__[9];
    rf[147] = rf[147] * c__[3] * c__[15];
    rf[148] = rf[148] * c__[5] * c__[7];
    rf[149] = rf[149] * c__[5] * c__[9];
    rf[150] *= c__[1];
    rf[151] *= c__[4];
    rf[152] *= c__[4];
    rf[154] *= c__[6];
    rf[155] *= c__[4];
    rf[156] *= c__[4];
    rf[157] = rf[157] * c__[2] * c__[17];
    rf[159] *= c__[4];
    rf[162] *= c__[2];
    rf[163] *= c__[5];
    rf[164] = rf[164] * c__[18] * c__[2];
    rf[165] = rf[165] * c__[18] * c__[3];
    rf[166] = rf[166] * c__[18] * c__[3];
    rf[167] *= c__[9];
    rb[1] = rb[1] * ctb[0] * c__[4];
    rb[2] = rb[2] * ctb[1] * c__[5];
    rb[3] = rb[3] * c__[2] * c__[5];
    rb[4] = rb[4] * c__[5] * c__[4];
    rb[5] = rb[5] * c__[5] * c__[7];
    rb[6] = rb[6] * c__[2] * c__[11];
    rb[7] *= c__[2];
    rb[8] = rb[8] * c__[1] * c__[11];
    rb[9] *= c__[2];
    rb[10] = rb[10] * c__[2] * c__[13];
    rb[11] = rb[11] * c__[5] * c__[9];
    rb[12] *= c__[12];
    rb[13] = rb[13] * c__[5] * c__[11];
    rb[14] = rb[14] * c__[2] * c__[12];
    rb[15] *= c__[5];
    rb[16] = rb[16] * c__[5] * c__[13];
    rb[17] *= c__[2];
    rb[18] *= c__[11];
    rb[19] = rb[19] * c__[2] * c__[17];
    rb[20] *= c__[9];
    rb[21] = rb[21] * c__[9] * c__[13];
    rb[22] *= c__[5];
    rb[23] = rb[23] * c__[2] * c__[11] * c__[11];
    rb[24] *= c__[5];
    rb[25] *= c__[12];
    rb[26] = rb[26] * c__[3] * c__[12];
    rb[27] *= c__[7];
    rb[28] = rb[28] * ctb[27] * c__[7];
    rb[29] = rb[29] * c__[7] * c__[4];
    rb[30] = rb[30] * c__[7] * c__[6];
    rb[31] = rb[31] * c__[7] * c__[19];
    rb[32] = rb[32] * c__[3] * c__[5];
    rb[33] = rb[33] * ctb[32] * c__[1];
    rb[34] = rb[34] * c__[1] * c__[1];
    rb[35] = rb[35] * c__[1] * c__[6];
    rb[36] = rb[36] * c__[1] * c__[12];
    rb[37] = rb[37] * ctb[36] * c__[6];
    rb[38] = rb[38] * c__[3] * c__[6];
    rb[39] = rb[39] * c__[4] * c__[1];
    rb[40] = rb[40] * c__[5] * c__[5];
    rb[41] = rb[41] * c__[7] * c__[1];
    rb[42] = rb[42] * c__[5] * c__[6];
    rb[43] *= c__[1];
    rb[44] *= c__[9];
    rb[45] *= c__[1];
    rb[46] *= c__[10];
    rb[47] = rb[47] * c__[9] * c__[1];
    rb[48] *= c__[13];
    rb[49] = rb[49] * c__[1] * c__[11];
    rb[51] *= c__[1];
    rb[52] = rb[52] * c__[1] * c__[13];
    rb[53] = rb[53] * c__[5] * c__[9];
    rb[54] *= c__[6];
    rb[56] *= c__[15];
    rb[57] = rb[57] * c__[1] * c__[14];
    rb[59] *= c__[1];
    rb[60] *= c__[16];
    rb[61] = rb[61] * c__[1] * c__[15];
    rb[62] *= c__[1];
    rb[63] *= c__[11];
    rb[64] *= c__[1];
    rb[65] = rb[65] * c__[9] * c__[11];
    rb[66] *= c__[13];
    rb[67] = rb[67] * c__[2] * c__[6];
    rb[68] *= c__[8];
    rb[69] = rb[69] * c__[3] * c__[6];
    rb[70] = rb[70] * c__[4] * c__[6];
    rb[71] = rb[71] * c__[7] * c__[6];
    rb[72] = rb[72] * c__[7] * c__[6];
    rb[73] = rb[73] * c__[2] * c__[11];
    rb[74] *= c__[2];
    rb[75] = rb[75] * c__[2] * c__[13];
    rb[76] *= c__[6];
    rb[77] = rb[77] * c__[2] * c__[13];
    rb[78] *= c__[6];
    rb[79] *= c__[6];
    rb[80] = rb[80] * c__[9] * c__[6];
    rb[81] = rb[81] * c__[2] * c__[12];
    rb[82] = rb[82] * c__[6] * c__[11];
    rb[83] *= c__[6];
    rb[84] = rb[84] * c__[6] * c__[13];
    rb[85] = rb[85] * c__[2] * c__[17];
    rb[86] = rb[86] * c__[9] * c__[11];
    rb[87] = rb[87] * c__[6] * c__[14];
    rb[88] *= c__[6];
    rb[89] *= c__[6];
    rb[90] *= c__[6];
    rb[91] = rb[91] * c__[4] * c__[8];
    rb[92] = rb[92] * c__[4] * c__[8];
    rb[93] = rb[93] * c__[5] * c__[13];
    rb[94] = rb[94] * c__[4] * c__[10];
    rb[95] *= c__[5];
    rb[96] = rb[96] * c__[5] * c__[12];
    rb[97] *= c__[8];
    rb[98] = rb[98] * c__[3] * c__[11];
    rb[99] = rb[99] * c__[2] * c__[14];
    rb[100] *= c__[3];
    rb[101] *= c__[2];
    rb[102] = rb[102] * c__[2] * c__[13];
    rb[103] = rb[103] * c__[2] * c__[14];
    rb[105] = rb[105] * c__[2] * c__[15];
    rb[107] *= c__[11];
    rb[108] = rb[108] * c__[2] * c__[17];
    rb[109] = rb[109] * c__[11] * c__[14];
    rb[111] = rb[111] * c__[2] * c__[9];
    rb[112] = rb[112] * c__[1] * c__[14];
    rb[113] = rb[113] * c__[2] * c__[15];
    rb[114] = rb[114] * c__[9] * c__[9];
    rb[115] *= c__[17];
    rb[117] *= c__[19];
    rb[118] = rb[118] * c__[2] * c__[5] * c__[11];
    rb[119] = rb[119] * c__[11] * c__[6];
    rb[120] = rb[120] * c__[9] * c__[2];
    rb[121] *= c__[6];
    rb[122] = rb[122] * c__[2] * c__[15];
    rb[123] = rb[123] * c__[9] * c__[9];
    rb[124] *= c__[11];
    rb[125] *= c__[12];
    rb[126] = rb[126] * c__[11] * c__[13];
    rb[128] *= c__[3];
    rb[129] = rb[129] * c__[5] * c__[13];
    rb[130] = rb[130] * c__[7] * c__[10];
    rb[131] *= c__[16];
    rb[132] *= c__[2];
    rb[133] = rb[133] * c__[10] * c__[11];
    rb[134] *= c__[10];
    rb[135] *= c__[10];
    rb[136] *= c__[10];
    rb[137] = rb[137] * c__[2] * c__[11] * c__[6];
    rb[138] = rb[138] * ctb[137] * c__[2] * c__[11];
    rb[139] = rb[139] * c__[7] * c__[11];
    rb[140] = rb[140] * c__[7] * c__[13];
    rb[142] = rb[142] * c__[1] * c__[14];
    rb[143] = rb[143] * c__[7] * c__[15];
    rb[144] = rb[144] * c__[5] * c__[11] * c__[11];
    rb[145] = rb[145] * c__[11] * c__[11] * c__[14];
    rb[147] *= c__[2];
    rb[148] = rb[148] * c__[4] * c__[6];
    rb[150] *= c__[9];
    rb[152] = rb[152] * c__[3] * c__[13];
    rb[156] = rb[156] * c__[7] * c__[14];
    rb[161] *= c__[9];
    rb[162] = rb[162] * c__[17] * c__[1];
    rb[163] = rb[163] * c__[6] * c__[17];
    rb[164] = rb[164] * c__[15] * c__[9];
    rb[165] = rb[165] * c__[17] * c__[9] * c__[2];
    rb[167] *= c__[18];

    return 0;
} /* ratx_ */