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

//#include "f2c.h"

/* Table of constant values */

#include <cmath>            //dol
#include <algorithm>        //dol
using namespace std;        //dol
typedef int    integer;     //dol
typedef double doublereal;  //dol

static doublereal c_b2 = 10.;


/*     A 13-species reduced mechanism for lean atmospheric CH4/air */

/*     Contact: */
/*     Tianfeng Lu */
/*     University of Connecticut */
/*     Email: tlu@engr.uconn.edu */

/*     Reference: */
/*     R. Sankaran, E.R. Hawkes, J.H. Chen, T.F. Lu, C.K. Law, */
/*     "Structure of a spatially developing turbulent lean methane–air Bunsen flame," */
/*     Proceedings of the Combustion Institute 31 (2007) 1291–1298. */

/* Subroutine */ 
//int ckwyp_(doublereal *p, doublereal *t, doublereal *y, 
//	   integer *ickwrk, doublereal *rckwrk, doublereal *wdot)
//{

void getProblemSpecificRR(double rho, double ttemp, double pres,
                          double *yi, double *rrsp) {
    /* Initialized data */

    doublereal *p    = &pres;
    doublereal *t    = &ttemp;
    doublereal *y    = yi;
    doublereal *wdot = rrsp;

    static doublereal rf[73] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
    static doublereal rb[73] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    //double log(doublereal), exp(doublereal),  pow(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__[13];
    static integer i__, k, n;
    static doublereal f1, f2, f3, f4, a12, a13, a21, a24, a31, fc, a42, eg[17]
	    , ti, pr, xn, xq[4], ff1, ff2, ti2, tn1, tn2, tn3, tn4, tn5, aa12,
	     aa21, ctb[73], eqk[73], smh[17], rop[73], tmp, sum, den1, den2, 
	    abv1, abv2, abv3, den3, abv4, den4, pfac, flog, pcor, temp, ctot, 
	    pfac2, pfac3, fclog, fcent, alogt, small, prlog, rklow[5], cprlog;


    /* Parameter adjustments */
    --wdot;
    //--rckwrk;
    //--ickwrk;
    --y;

    /* Function Body */

    small = 1e-50;

/* 	 convert Y to C */

    c__[0] = max(y[1],small) * .496046521;
    c__[1] = max(y[2],small) * .992093043;
    c__[2] = max(y[3],small) * .0625023433;
    c__[3] = max(y[4],small) * .0312511716;
    c__[4] = max(y[5],small) * .0587980383;
    c__[5] = max(y[6],small) * .0555082499;
    c__[6] = max(y[7],small) * .0302968146;
    c__[7] = max(y[8],small) * .0665112065;
    c__[8] = max(y[9],small) * .0623323639;
    c__[9] = max(y[10],small) * .0357008335;
    c__[10] = max(y[11],small) * .0227221341;
    c__[11] = max(y[12],small) * .0333039255;
    c__[12] = max(y[13],small) * .0356972032;

    sum = 0.f;
    for (k = 1; k <= 13; ++k) {
	sum += c__[k - 1];
    }
    sum = *p / (sum * *t * 83145100.);

    for (k = 1; k <= 13; ++k) {
	c__[k - 1] *= sum;
    }

/* 	 forward reaction rates */

    alogt = log(*t);
    ti = 1. / *t;
    ti2 = ti * ti;

    rf[0] = ti * 1.2e17;
    rf[1] = ti * 5e17;
    rf[2] = exp(alogt * 2.67 + 10.8197783 - ti * 3165.23284);
    rf[3] = 2e13;
    rf[4] = 8e13;
    rf[5] = 1.5e13;
    rf[6] = 8.43e13;
    rf[7] = exp(alogt * 1.5 + 20.7430685 - ti * 4327.66334);
    rf[8] = exp(34.0312786 - ti * 1509.65);
    rf[9] = 3e13;
    rf[10] = 3e13;
    rf[11] = exp(31.2945828 - ti * 1781.387);
    rf[12] = 1e13;
    rf[13] = exp(28.5473118 - ti * 24053.7567);
    rf[14] = exp(32.2361913 - ti * 20128.6667);
    rf[15] = exp(42.4761511 - alogt * .86);
    tmp = exp(alogt * -1.72);
    rf[16] = tmp * 3e20;
    rf[18] = tmp * 3.75e20;
    rf[17] = exp(43.6851114 - alogt * .76);
    rf[19] = exp(32.0498617 - ti * 7252.86183);
    rf[20] = ti * 1e18;
    rf[21] = exp(39.0385861 - alogt * .6);
    rf[22] = exp(45.5408762 - alogt * 1.25);
    rf[23] = ti2 * 5.5e20;
    rf[24] = ti2 * 2.2e22;
    rf[25] = exp(29.0097872 - ti * 337.658384);
    rf[26] = exp(30.9632256 - ti * 537.435401);
    rf[27] = exp(32.5288609 - ti * 319.542584);
    rf[28] = exp(37.7576522 - alogt * .8);
    rf[29] = exp(37.0803784 - alogt * .63 - ti * 192.731984);
    rf[30] = exp(alogt * 1.62 + 20.3077504 - ti * 5454.86868);
    rf[31] = exp(alogt * .48 + 27.7171988 + ti * 130.836334);
    rf[32] = 7.34e13;
    rf[33] = exp(alogt * .454 + 27.014835 - ti * 1811.58);
    rf[34] = exp(alogt * 1.05 + 23.8587601 - ti * 1648.03459);
    rf[35] = 2e13;
    rf[36] = 1.2e13;
    rf[37] = 6e12;
    rf[38] = exp(alogt * 1.5 + 17.5767107 - ti * 40056.0467);
    rf[39] = exp(alogt * 1.51 + 19.190789 - ti * 1726.03317);
    rf[40] = exp(alogt * 2.4 + 10.482906 + ti * 1061.78717);
    rf[41] = exp(ti * 251.608334 + 30.9983169);
    rf[42] = 2e13;
    rf[43] = 3e13;
    rf[44] = exp(alogt * 1.6 + 17.8408622 - ti * 2727.43434);
    rf[45] = 2.501e13;
    rf[46] = exp(alogt * 1.6 + 18.4206807 - ti * 1570.036);
    rf[47] = exp(alogt * 1.228 + 17.6783433 - ti * 35.2251667);
    rf[48] = 5e13;
    rf[49] = exp(alogt * 1.18 + 21.9558261 + ti * 224.93785);
    rf[50] = 5e12;
    rf[51] = 2e13;
    rf[52] = 1e12;
    rf[53] = exp(32.6416564 - ti * 11875.9134);
    rf[54] = exp(30.2112379 - ti * 754.825001);
    rf[55] = exp(alogt * 2. + 13.1223634 - ti * 3638.25651);
    rf[56] = exp(alogt * 2. + 14.7156719 - ti * 4161.60184);
    rf[57] = exp(30.3390713 - ti * 301.930001);
    rf[58] = 2.8e13;
    rf[59] = 1.2e13;
    rf[60] = 7e13;
    rf[61] = 3e13;
    rf[62] = exp(ti * 286.833501 + 30.4036098);
    rf[63] = 9e12;
    rf[64] = 7e12;
    rf[65] = 1.4e13;
    rf[66] = exp(24.3067848 - ti * 4498.75701);
    rf[67] = 2.648e13;
    rf[68] = exp(alogt * 2.81 + 8.10772006 - ti * 2948.84967);
    tmp = exp(alogt * -1. - ti * 8554.68335);
    rf[69] = tmp * 2.244e18;
    rf[70] = tmp * 1.87e17;
    rf[71] = exp(29.6591694 - ti * 201.286667);
    rf[72] = exp(30.5213929 - ti * 452.895001);

/* 	 thermal data */

    tn1 = alogt - 1.f;
    tn2 = *t;
    tn3 = tn2 * *t;
    tn4 = tn3 * *t;
    tn5 = tn4 * *t;

    if (*t > 1e3) {

	smh[0] = ti * 950.158922 - 3.20502331 + tn1 * 3.3372792 - tn2 * 
		2.47012365e-5 + tn3 * 8.32427963e-8 - tn4 * 1.49638662e-11 + 
		tn5 * 1.00127688e-15;
	smh[1] = -.446682914 - ti * 25473.6599 + tn1 * 2.50000001 - tn2 * 
		1.15421486e-11 + tn3 * 2.69269913e-15 - tn4 * 3.94596029e-19 
		+ tn5 * 2.49098679e-23;
	smh[2] = 4.78433864 - ti * 29217.5791 + tn1 * 2.56942078 - tn2 * 
		4.29870569e-5 + tn3 * 6.99140982e-9 - tn4 * 8.34814992e-13 + 
		tn5 * 6.14168455e-17;
	smh[3] = ti * 1088.45772 + 5.45323129 + tn1 * 3.28253784 + tn2 * 
		7.4154377e-4 - tn3 * 1.26327778e-7 + tn4 * 1.74558796e-11 - 
		tn5 * 1.08358897e-15;
	smh[4] = 4.4766961 - ti * 3858.657 + tn1 * 3.09288767 + tn2 * 
		2.74214858e-4 + tn3 * 2.10842047e-8 - tn4 * 7.3288463e-12 + 
		tn5 * 5.8706188e-16;
	smh[5] = ti * 30004.2971 + 4.9667701 + tn1 * 3.03399249 + tn2 * 
		.00108845902 - tn3 * 2.73454197e-8 - tn4 * 8.08683225e-12 + 
		tn5 * 8.4100496e-16;
	smh[6] = 3.78510215 - ti * 111.856713 + tn1 * 4.0172109 + tn2 * 
		.00111991007 - tn3 * 1.05609692e-7 + tn4 * 9.52053083e-12 - 
		tn5 * 5.39542675e-16;
	smh[7] = 6.17119324 - ti * 46263.604 + tn1 * 2.87410113 + tn2 * 
		.00182819646 - tn3 * 2.34824328e-7 + tn4 * 2.16816291e-11 - 
		tn5 * 9.38637835e-16;
	smh[8] = 8.62650169 - ti * 50925.9997 + tn1 * 2.29203842 + tn2 * 
		.00232794319 - tn3 * 3.35319912e-7 + tn4 * 3.48255e-11 - tn5 *
		 1.69858183e-15;
	smh[9] = 8.48007179 - ti * 16775.5843 + tn1 * 2.28571772 + tn2 * 
		.00361995018 - tn3 * 4.97857247e-7 + tn4 * 4.9640387e-11 - 
		tn5 * 2.33577197e-15;
	smh[10] = ti * 9468.34459 + 18.437318 + tn1 * .074851495 + tn2 * 
		.00669547335 - tn3 * 9.55476348e-7 + tn4 * 1.01910446e-10 - 
		tn5 * 5.0907615e-15;
	smh[11] = ti * 14151.8724 + 7.81868772 + tn1 * 2.71518561 + tn2 * 
		.00103126372 - tn3 * 1.66470962e-7 + tn4 * 1.9171084e-11 - 
		tn5 * 1.01823858e-15;
	smh[12] = ti * 48759.166 + 2.27163806 + tn1 * 3.85746029 + tn2 * 
		.00220718513 - tn3 * 3.69135673e-7 + tn4 * 4.36241823e-11 - 
		tn5 * 2.36042082e-15;
	smh[13] = 9.79834492 - ti * 4011.91815 + tn1 * 2.77217438 + tn2 * 
		.00247847763 - tn3 * 4.14076022e-7 + tn4 * 4.90968148e-11 - 
		tn5 * 2.66754356e-15;
	smh[14] = ti * 13995.8323 + 13.656323 + tn1 * 1.76069008 + tn2 * 
		.00460000041 - tn3 * 7.37098022e-7 + tn4 * 8.38676767e-11 - 
		tn5 * 4.4192782e-15;
	smh[15] = ti * 3242.50627 + 5.81043215 + tn1 * 3.69266569 + tn2 * 
		.00432288399 - tn3 * 6.25168533e-7 + tn4 * 6.56028863e-11 - 
		tn5 * 3.24277101e-15;

    } else {

	smh[0] = ti * 917.935173 + .683010238 + tn1 * 2.34433112 + tn2 * 
		.00399026037 - tn3 * 3.2463585e-6 + tn4 * 1.67976745e-9 - tn5 
		* 3.68805881e-13;
	smh[1] = -.446682853 - ti * 25473.6599 + tn1 * 2.5 + tn2 * 
		3.52666409e-13 - tn3 * 3.32653273e-16 + tn4 * 1.91734693e-19 
		- tn5 * 4.63866166e-23;
	smh[2] = 2.05193346 - ti * 29122.2592 + tn1 * 3.1682671 - tn2 * 
		.00163965942 + tn3 * 1.10717733e-6 - tn4 * 5.10672187e-10 + 
		tn5 * 1.05632986e-13;
	smh[3] = ti * 1063.94356 + 3.65767573 + tn1 * 3.78245636 - tn2 * 
		.00149836708 + tn3 * 1.641217e-6 - tn4 * 8.06774591e-10 + tn5 
		* 1.62186419e-13;
	smh[4] = -.103925458 - ti * 3615.08056 + tn1 * 3.99201543 - tn2 * 
		.00120065876 + tn3 * 7.69656402e-7 - tn4 * 3.23427778e-10 + 
		tn5 * 6.8205735e-14;
	smh[5] = ti * 30293.7267 - .849032208 + tn1 * 4.19864056 - tn2 * 
		.00101821705 + tn3 * 1.08673369e-6 - tn4 * 4.57330885e-10 + 
		tn5 * 8.85989085e-14;
	smh[6] = 3.71666245 - ti * 294.80804 + tn1 * 4.30179801 - tn2 * 
		.00237456025 + tn3 * 3.52638152e-6 - tn4 * 2.02303245e-9 + 
		tn5 * 4.64612562e-13;
	smh[7] = 1.56253185 - ti * 46004.0401 + tn1 * 3.76267867 + tn2 * 
		4.84436072e-4 + tn3 * 4.65816402e-7 - tn4 * 3.20909294e-10 + 
		tn5 * 8.43708595e-14;
	smh[8] = -.769118967 - ti * 50496.8163 + tn1 * 4.19860411 - tn2 * 
		.0011833071 + tn3 * 1.37216037e-6 - tn4 * 5.57346651e-10 + 
		tn5 * 9.71573685e-14;
	smh[9] = 1.60456433 - ti * 16444.9988 + tn1 * 3.6735904 + tn2 * 
		.00100547588 + tn3 * 9.55036427e-7 - tn4 * 5.72597854e-10 + 
		tn5 * 1.27192867e-13;
	smh[10] = ti * 10246.6476 - 4.64130376 + tn1 * 5.14987613 - tn2 * 
		.0068354894 + tn3 * 8.19667665e-6 - tn4 * 4.03952522e-9 + tn5 
		* 8.3346978e-13;
	smh[11] = ti * 14344.086 + 3.50840928 + tn1 * 3.57953347 - tn2 * 
		3.0517684e-4 + tn3 * 1.69469055e-7 + tn4 * 7.55838237e-11 - 
		tn5 * 4.52212249e-14;
	smh[12] = ti * 48371.9697 + 9.90105222 + tn1 * 2.35677352 + tn2 * 
		.00449229839 - tn3 * 1.18726045e-6 + tn4 * 2.04932518e-10 - 
		tn5 * 7.1849774e-15;
	smh[13] = 3.39437243 - ti * 3839.56496 + tn1 * 4.22118584 - tn2 * 
		.00162196266 + tn3 * 2.29665743e-6 - tn4 * 1.10953411e-9 + 
		tn5 * 2.16884433e-13;
	smh[14] = ti * 14308.9567 + .6028129 + tn1 * 4.79372315 - tn2 * 
		.00495416685 + tn3 * 6.22033347e-6 - tn4 * 3.16071051e-9 + 
		tn5 * 6.5886326e-13;
	smh[15] = ti * 3193.91367 + 5.47302243 + tn1 * 3.86388918 + tn2 * 
		.00279836152 + tn3 * 9.88786318e-7 - tn4 * 8.711001e-10 + tn5 
		* 2.18483639e-13;
    }

/* 	 equilibrium constants */

    for (n = 1; n <= 16; ++n) {
/* Computing MIN */
	d__1 = smh[n - 1];
	eg[n - 1] = exp((min(d__1,230.26)));
    }

    pfac = 1013250. / (*t * 83145100.);
    pfac2 = pfac * pfac;
    pfac3 = pfac2 * pfac;

    eqk[0] = eg[3] / eg[2] / eg[2] / pfac;
    eqk[1] = eg[4] / eg[1] / eg[2] / pfac;
    eqk[2] = eg[1] * eg[4] / eg[0] / eg[2];
    eqk[3] = eg[3] * eg[4] / eg[2] / eg[6];
    eqk[4] = eg[1] * eg[13] / eg[2] / eg[7];
    eqk[5] = eg[0] * eg[11] / eg[2] / eg[8];
    eqk[6] = eg[1] * eg[14] / eg[2] / eg[9];
    eqk[7] = eg[4] * eg[9] / eg[2] / eg[10];
    eqk[8] = eg[12] / eg[2] / eg[11] / pfac;
    eqk[9] = eg[4] * eg[11] / eg[2] / eg[13];
    eqk[10] = eg[1] * eg[12] / eg[2] / eg[13];
    eqk[11] = eg[4] * eg[13] / eg[2] / eg[14];
    eqk[12] = eg[4] * eg[14] / eg[2] / eg[15];
    eqk[13] = eg[2] * eg[12] / eg[3] / eg[11];
    eqk[14] = eg[6] * eg[13] / eg[3] / eg[14];
    eqk[15] = eg[6] / eg[1] / eg[3] / pfac;
    eqk[16] = eqk[15];
    eqk[17] = eqk[15];
    eqk[18] = eqk[15];
    eqk[19] = eg[2] * eg[4] / eg[1] / eg[3];
    eqk[20] = eg[0] / eg[1] / eg[1] / pfac;
    eqk[21] = eqk[20];
    eqk[22] = eqk[20];
    eqk[23] = eqk[20];
    eqk[24] = eg[5] / eg[1] / eg[4] / pfac;
    eqk[25] = eg[2] * eg[5] / eg[1] / eg[6];
    eqk[26] = eg[0] * eg[3] / eg[1] / eg[6];
    eqk[27] = eg[4] * eg[4] / eg[1] / eg[6];
    eqk[28] = eg[9] / eg[1] / eg[7] / pfac;
    eqk[29] = eg[10] / eg[1] / eg[9] / pfac;
    eqk[30] = eg[0] * eg[9] / eg[1] / eg[10];
    eqk[31] = eg[14] / eg[1] / eg[13] / pfac;
    eqk[32] = eg[0] * eg[11] / eg[1] / eg[13];
    eqk[33] = eg[15] / eg[1] / eg[14] / pfac;
    eqk[34] = eg[0] * eg[13] / eg[1] / eg[14];
    eqk[35] = eg[0] * eg[14] / eg[1] / eg[15];
    eqk[36] = eg[4] * eg[9] / eg[1] / eg[15];
    eqk[37] = eg[5] * eg[8] / eg[1] / eg[15];
    eqk[38] = eg[14] / eg[0] / eg[11] / pfac;
    eqk[39] = eg[1] * eg[5] / eg[0] / eg[4];
    eqk[40] = eg[2] * eg[5] / eg[4] / eg[4];
    eqk[41] = eg[3] * eg[5] / eg[4] / eg[6];
    eqk[42] = eg[1] * eg[14] / eg[4] / eg[7];
    eqk[43] = eg[1] * eg[14] / eg[4] / eg[8];
    eqk[44] = eg[5] * eg[7] / eg[4] / eg[9];
    eqk[45] = eg[5] * eg[8] / eg[4] / eg[9];
    eqk[46] = eg[5] * eg[9] / eg[4] / eg[10];
    eqk[47] = eg[1] * eg[12] / eg[4] / eg[11];
    eqk[48] = eg[5] * eg[11] / eg[4] / eg[13];
    eqk[49] = eg[5] * eg[13] / eg[4] / eg[14];
    eqk[50] = eg[5] * eg[14] / eg[4] / eg[15];
    eqk[51] = eg[4] * eg[14] / eg[6] / eg[7];
    eqk[52] = eg[3] * eg[10] / eg[6] / eg[9];
    eqk[53] = eg[4] * eg[12] / eg[6] / eg[11];
    eqk[54] = eg[4] * eg[13] / eg[3] / eg[7];
    eqk[55] = eg[1] * eg[9] / eg[0] / eg[7];
    eqk[56] = eg[9] * eg[9] / eg[7] / eg[10];
    eqk[57] = eg[7] / eg[8];
    eqk[61] = eqk[57];
    eqk[63] = eqk[57];
    eqk[64] = eqk[57];
    eqk[58] = eg[1] * eg[4] * eg[11] / eg[3] / eg[8] * pfac;
    eqk[59] = eg[5] * eg[11] / eg[3] / eg[8];
    eqk[60] = eg[1] * eg[9] / eg[0] / eg[8];
    eqk[62] = eg[9] * eg[9] / eg[8] / eg[10];
    eqk[65] = eg[11] * eg[14] / eg[8] / eg[12];
    eqk[66] = eg[4] * eg[14] / eg[3] / eg[9];
    eqk[67] = eg[10] * eg[11] / eg[9] / eg[13];
    eqk[68] = eg[10] * eg[13] / eg[9] / eg[14];
    eqk[69] = eg[1] * eg[11] / eg[13] * pfac;
    eqk[70] = eqk[69];
    eqk[71] = eg[6] * eg[11] / eg[3] / eg[13];
    eqk[72] = eg[6] * eg[14] / eg[3] / eg[15];

/* 	Compute reverse reaction rates */

    for (i__ = 1; i__ <= 73; ++i__) {
/* Computing MAX */
	d__1 = eqk[i__ - 1];
	rb[i__ - 1] = rf[i__ - 1] / max(d__1,1e-100);
    }

/* 	 rates at low pressure limit */

    rklow[0] = exp(63.3329483 - alogt * 3.14 - 618.956501 / *t);
    rklow[1] = exp(76.8923562 - alogt * 4.76 - 1227.84867 / *t);
    rklow[2] = exp(55.5621468 - alogt * 2.57 - 717.083751 / *t);
    rklow[3] = exp(73.9217399 - alogt * 4.82 - 3286.00484 / *t);
    rklow[4] = exp(63.7931383 - alogt * 3.42 - 42446.3259 / *t);

/*      third-body concentrations */

    ctot = 0.f;
    for (k = 1; k <= 13; ++k) {
	ctot += c__[k - 1];
    }

    ctb[0] = ctot + c__[0] * 1.4 + c__[5] * 14.4 + c__[8] + c__[9] * .75 + 
	    c__[10] * 2.6;
    ctb[1] = ctot + c__[0] + c__[5] * 5. + c__[8] + c__[9] * .5 + c__[10];
    ctb[8] = ctot + c__[0] + c__[3] * 5. + c__[5] * 5. + c__[8] + c__[9] * .5 
	    + c__[10] * 2.5;
    ctb[15] = ctot - c__[3] - c__[5] - c__[9] * .25 + c__[10] * .5 - c__[12];
    ctb[20] = ctot - c__[0] - c__[5] + c__[8] - c__[10];
    ctb[24] = ctot - c__[0] * .27 + c__[5] * 2.65 + c__[8];
    ctb[28] = ctot + c__[0] + c__[5] * 5. + c__[8] + c__[9] * .5 + c__[10];
    ctb[29] = ctb[28];
    ctb[31] = ctb[28];
    ctb[33] = ctb[28];
    ctb[38] = ctb[28];
    ctb[70] = ctot + c__[0] - c__[5] + c__[8] + c__[9] * .5 + c__[10];

/*      fall-off reactions */


    pr = rklow[0] * ctb[28] / rf[28];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 78.) * .32 + exp(-(*t) / 1995.) * .68 + exp(-5590. / *
	    t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b2, flog);
    pcor = fc * pcor;
    rf[28] *= pcor;
    rb[28] *= pcor;

    pr = rklow[1] * ctb[29] / rf[29];
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
    fc = pow(c_b2, flog);
    pcor = fc * pcor;
    rf[29] *= pcor;
    rb[29] *= pcor;

    pr = rklow[2] * ctb[31] / rf[31];
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
    fc = pow(c_b2, flog);
    pcor = fc * pcor;
    rf[31] *= pcor;
    rb[31] *= pcor;

    pr = rklow[3] * ctb[33] / rf[33];
    pcor = pr / (pr + 1.f);
    d__1 = max(pr,small);
    prlog = log10(d__1);
    fcent = exp(-(*t) / 103.) * .2813 + exp(-(*t) / 1291.) * .7187 + exp(
	    -4160. / *t);
    d__1 = max(fcent,small);
    fclog = log10(d__1);
    xn = .75f - fclog * 1.27f;
    cprlog = prlog - (fclog * .67f + .4f);
/* Computing 2nd power */
    d__1 = cprlog / (xn - cprlog * .14f);
    flog = fclog / (d__1 * d__1 + 1.f);
    fc = pow(c_b2, flog);
    pcor = fc * pcor;
    rf[33] *= pcor;
    rb[33] *= pcor;

    pr = rklow[4] * ctb[38] / rf[38];
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
    fc = pow(c_b2, flog);
    pcor = fc * pcor;
    rf[38] *= pcor;
    rb[38] *= pcor;

    rf[0] = rf[0] * ctb[0] * c__[2] * c__[2];
    rf[1] = rf[1] * ctb[1] * c__[2] * c__[1];
    rf[2] = rf[2] * c__[2] * c__[0];
    rf[3] = rf[3] * c__[2] * c__[6];
    rf[4] *= c__[2];
    rf[5] *= c__[2];
    rf[6] = rf[6] * c__[2] * c__[7];
    rf[7] = rf[7] * c__[2] * c__[8];
    rf[8] = rf[8] * ctb[8] * c__[2] * c__[9];
    rf[9] *= c__[2];
    rf[10] *= c__[2];
    rf[11] = rf[11] * c__[2] * c__[11];
    rf[12] *= c__[2];
    rf[13] = rf[13] * c__[3] * c__[9];
    rf[14] = rf[14] * c__[3] * c__[11];
    rf[15] = rf[15] * ctb[15] * c__[1] * c__[3];
    rf[16] = rf[16] * c__[1] * c__[3] * c__[3];
    rf[17] = rf[17] * c__[1] * c__[3] * c__[5];
    rf[18] = rf[18] * c__[1] * c__[3] * c__[12];
    rf[19] = rf[19] * c__[1] * c__[3];
    rf[20] = rf[20] * ctb[20] * c__[1] * c__[1];
    rf[21] = rf[21] * c__[1] * c__[1] * c__[0];
    rf[22] = rf[22] * c__[1] * c__[1] * c__[5];
    rf[23] = rf[23] * c__[1] * c__[1] * c__[10];
    rf[24] = rf[24] * ctb[24] * c__[1] * c__[4];
    rf[25] = rf[25] * c__[1] * c__[6];
    rf[26] = rf[26] * c__[1] * c__[6];
    rf[27] = rf[27] * c__[1] * c__[6];
    rf[28] *= c__[1];
    rf[29] = rf[29] * c__[1] * c__[7];
    rf[30] = rf[30] * c__[1] * c__[8];
    rf[31] *= c__[1];
    rf[32] *= c__[1];
    rf[33] = rf[33] * c__[1] * c__[11];
    rf[34] = rf[34] * c__[1] * c__[11];
    rf[35] *= c__[1];
    rf[36] *= c__[1];
    rf[37] *= c__[1];
    rf[38] = rf[38] * c__[0] * c__[9];
    rf[39] = rf[39] * c__[4] * c__[0];
    rf[40] = rf[40] * c__[4] * c__[4];
    rf[41] = rf[41] * c__[4] * c__[6];
    rf[42] *= c__[4];
    rf[43] *= c__[4];
    rf[44] = rf[44] * c__[4] * c__[7];
    rf[45] = rf[45] * c__[4] * c__[7];
    rf[46] = rf[46] * c__[4] * c__[8];
    rf[47] = rf[47] * c__[4] * c__[9];
    rf[48] *= c__[4];
    rf[49] = rf[49] * c__[4] * c__[11];
    rf[50] *= c__[4];
    rf[51] *= c__[6];
    rf[52] = rf[52] * c__[6] * c__[7];
    rf[53] = rf[53] * c__[6] * c__[9];
    rf[54] *= c__[3];
    rf[55] *= c__[0];
    rf[56] *= c__[8];
    rf[57] *= c__[12];
    rf[58] *= c__[3];
    rf[59] *= c__[3];
    rf[60] *= c__[0];
    rf[61] *= c__[5];
    rf[62] *= c__[8];
    rf[63] *= c__[9];
    rf[64] *= c__[10];
    rf[65] *= c__[10];
    rf[66] = rf[66] * c__[7] * c__[3];
    rf[67] *= c__[7];
    rf[68] = rf[68] * c__[7] * c__[11];
    rf[69] *= c__[5];
    rf[70] *= ctb[70];
    rf[71] *= c__[3];
    rf[72] *= c__[3];
    rb[0] = rb[0] * ctb[0] * c__[3];
    rb[1] = rb[1] * ctb[1] * c__[4];
    rb[2] = rb[2] * c__[1] * c__[4];
    rb[3] = rb[3] * c__[4] * c__[3];
    rb[4] *= c__[1];
    rb[5] = rb[5] * c__[0] * c__[9];
    rb[6] = rb[6] * c__[1] * c__[11];
    rb[7] = rb[7] * c__[4] * c__[7];
    rb[8] = rb[8] * ctb[8] * c__[10];
    rb[9] = rb[9] * c__[4] * c__[9];
    rb[10] = rb[10] * c__[1] * c__[10];
    rb[11] *= c__[4];
    rb[12] = rb[12] * c__[4] * c__[11];
    rb[13] = rb[13] * c__[2] * c__[10];
    rb[14] *= c__[6];
    rb[15] = rb[15] * ctb[15] * c__[6];
    rb[16] = rb[16] * c__[6] * c__[3];
    rb[17] = rb[17] * c__[6] * c__[5];
    rb[18] = rb[18] * c__[6] * c__[12];
    rb[19] = rb[19] * c__[2] * c__[4];
    rb[20] = rb[20] * ctb[20] * c__[0];
    rb[21] = rb[21] * c__[0] * c__[0];
    rb[22] = rb[22] * c__[0] * c__[5];
    rb[23] = rb[23] * c__[0] * c__[10];
    rb[24] = rb[24] * ctb[24] * c__[5];
    rb[25] = rb[25] * c__[2] * c__[5];
    rb[26] = rb[26] * c__[3] * c__[0];
    rb[27] = rb[27] * c__[4] * c__[4];
    rb[28] *= c__[7];
    rb[29] *= c__[8];
    rb[30] = rb[30] * c__[7] * c__[0];
    rb[31] *= c__[11];
    rb[32] = rb[32] * c__[0] * c__[9];
    rb[34] *= c__[0];
    rb[35] = rb[35] * c__[0] * c__[11];
    rb[36] = rb[36] * c__[4] * c__[7];
    rb[37] *= c__[5];
    rb[38] *= c__[11];
    rb[39] = rb[39] * c__[1] * c__[5];
    rb[40] = rb[40] * c__[2] * c__[5];
    rb[41] = rb[41] * c__[3] * c__[5];
    rb[42] = rb[42] * c__[1] * c__[11];
    rb[43] = rb[43] * c__[1] * c__[11];
    rb[44] *= c__[5];
    rb[45] *= c__[5];
    rb[46] = rb[46] * c__[7] * c__[5];
    rb[47] = rb[47] * c__[1] * c__[10];
    rb[48] = rb[48] * c__[5] * c__[9];
    rb[49] *= c__[5];
    rb[50] = rb[50] * c__[5] * c__[11];
    rb[51] = rb[51] * c__[4] * c__[11];
    rb[52] = rb[52] * c__[3] * c__[8];
    rb[53] = rb[53] * c__[4] * c__[10];
    rb[54] *= c__[4];
    rb[55] = rb[55] * c__[1] * c__[7];
    rb[56] = rb[56] * c__[7] * c__[7];
    rb[57] *= c__[12];
    rb[58] = rb[58] * c__[1] * c__[4] * c__[9];
    rb[59] = rb[59] * c__[9] * c__[5];
    rb[60] = rb[60] * c__[7] * c__[1];
    rb[61] *= c__[5];
    rb[62] = rb[62] * c__[7] * c__[7];
    rb[63] *= c__[9];
    rb[64] *= c__[10];
    rb[65] = rb[65] * c__[9] * c__[11];
    rb[66] = rb[66] * c__[4] * c__[11];
    rb[67] = rb[67] * c__[8] * c__[9];
    rb[68] *= c__[8];
    rb[69] = rb[69] * c__[1] * c__[9] * c__[5];
    rb[70] = rb[70] * ctb[70] * c__[1] * c__[9];
    rb[71] = rb[71] * c__[6] * c__[9];
    rb[72] = rb[72] * c__[6] * c__[11];

/* 	 solving QSS species concentration */
/* 	------------------------------------------- */
/*     CH2 */
    abv1 = rb[28] + rb[42] + rf[44] + rb[51] + rb[55] + rb[56];
    den1 = rf[4] + rf[28] + rf[42] + rb[44] + rf[51] + rf[54] + rf[55] + rf[
	    56] + rb[57] + rb[61] + rb[63] + rb[64];

/*     CH2(S) */
    abv2 = rb[5] + rb[43] + rf[45] + rb[58] + rb[59] + rb[60] + rb[62] + rb[
	    65];
    den2 = rf[5] + rb[37] + rf[43] + rb[45] + rf[57] + rf[58] + rf[59] + rf[
	    60] + rf[61] + rf[62] + rf[63] + rf[64] + rf[65];

/*     HCO */
    abv3 = rb[9] + rb[10] + rf[11] + rf[14] + rb[31] + rb[32] + rf[34] + rb[
	    48] + rf[49] + rb[67] + rf[68] + rb[69] + rb[70] + rb[71];
    den3 = rb[4] + rf[9] + rf[10] + rb[11] + rb[14] + rf[31] + rf[32] + rb[34]
	     + rf[48] + rb[49] + rb[54] + rf[67] + rb[68] + rf[69] + rf[70] + 
	    rf[71];

/*     CH2OH */
    abv4 = rb[12] + rf[33] + rb[35] + rb[36] + rb[50] + rb[72];
    den4 = rf[12] + rb[33] + rf[35] + rf[36] + rf[37] + rf[50] + rf[72];

    f1 = abv1 / den1;
    a12 = (rf[57] + rf[61] + rf[63] + rf[64]) / den1;
    a13 = (rb[4] + rb[54]) / den1;

    f2 = abv2 / den2;
    a21 = (rb[57] + rb[61] + rb[63] + rb[64]) / den2;
    a24 = rf[37] / den2;

    f3 = abv3 / den3;
    a31 = (rf[4] + rf[54]) / den3;

    f4 = abv4 / den4;
    a42 = rb[37] / den4;

    temp = 1 - a13 * a31;
    ff1 = (f1 + a13 * f3) / temp;
    aa12 = a12 / temp;

    temp = 1 - a24 * a42;
    ff2 = (f2 + a24 * f4) / temp;
    aa21 = a21 / temp;

    xq[0] = (ff1 + aa12 * ff2) / (1 - aa12 * aa21);
    xq[1] = ff2 + aa21 * xq[0];
    xq[2] = f3 + a31 * xq[0];
    xq[3] = f4 + a42 * xq[1];

/* 	--------------------------------------------- */

/* 	 update rates of reactions involving QSS species */

    rf[4] *= xq[0];
    rb[4] *= xq[2];
    rf[5] *= xq[1];
    rf[9] *= xq[2];
    rf[10] *= xq[2];
    rb[11] *= xq[2];
    rf[12] *= xq[3];
    rb[14] *= xq[2];
    rf[28] *= xq[0];
    rf[31] *= xq[2];
    rf[32] *= xq[2];
    rb[33] *= xq[3];
    rb[34] *= xq[2];
    rf[35] *= xq[3];
    rf[36] *= xq[3];
    rf[37] *= xq[3];
    rb[37] *= xq[1];
    rf[42] *= xq[0];
    rf[43] *= xq[1];
    rb[44] *= xq[0];
    rb[45] *= xq[1];
    rf[48] *= xq[2];
    rb[49] *= xq[2];
    rf[50] *= xq[3];
    rf[51] *= xq[0];
    rf[54] *= xq[0];
    rb[54] *= xq[2];
    rf[55] *= xq[0];
    rf[56] *= xq[0];
    rf[57] *= xq[1];
    rb[57] *= xq[0];
    rf[58] *= xq[1];
    rf[59] *= xq[1];
    rf[60] *= xq[1];
    rf[61] *= xq[1];
    rb[61] *= xq[0];
    rf[62] *= xq[1];
    rf[63] *= xq[1];
    rb[63] *= xq[0];
    rf[64] *= xq[1];
    rb[64] *= xq[0];
    rf[65] *= xq[1];
    rf[67] *= xq[2];
    rb[68] *= xq[2];
    rf[69] *= xq[2];
    rf[70] *= xq[2];
    rf[71] *= xq[2];
    rf[72] *= xq[3];

/* 	Compute WDOT */

    for (i__ = 1; i__ <= 73; ++i__) {
	rop[i__ - 1] = rf[i__ - 1] - rb[i__ - 1];
    }

/* 	Combine duplicated reactions */

    rop[15] = rop[15] + rop[16] + rop[17] + rop[18];
    rop[20] = rop[20] + rop[21] + rop[22] + rop[23];
    rop[57] = rop[57] + rop[61] + rop[63] + rop[64];
    rop[69] += rop[70];

    wdot[1] = -rop[2] + rop[5] + rop[20] + rop[26] + rop[30] + rop[32] + rop[
	    34] + rop[35] - rop[38] - rop[39] - rop[55] - rop[60];
    wdot[2] = -rop[1] + rop[2] + rop[4] + rop[6] + rop[10] - rop[15] - rop[19]
	     - rop[20] - rop[20] - rop[24] - rop[25] - rop[26] - rop[27] - 
	    rop[28] - rop[29] - rop[30] - rop[31] - rop[32] - rop[33] - rop[
	    34] - rop[35] - rop[36] - rop[37] + rop[39] + rop[42] + rop[43] + 
	    rop[47] + rop[55] + rop[58] + rop[60] + rop[69];
    wdot[3] = -rop[0] - rop[0] - rop[1] - rop[2] - rop[3] - rop[4] - rop[5] - 
	    rop[6] - rop[7] - rop[8] - rop[9] - rop[10] - rop[11] - rop[12] + 
	    rop[13] + rop[19] + rop[25] + rop[40];
    wdot[4] = rop[0] + rop[3] - rop[13] - rop[14] - rop[15] - rop[19] + rop[
	    26] + rop[41] + rop[52] - rop[54] - rop[58] - rop[59] - rop[66] - 
	    rop[71] - rop[72];
    wdot[5] = rop[1] + rop[2] + rop[3] + rop[7] + rop[9] + rop[11] + rop[12] 
	    + rop[19] - rop[24] + rop[27] + rop[27] + rop[36] - rop[39] - rop[
	    40] - rop[40] - rop[41] - rop[42] - rop[43] - rop[44] - rop[45] - 
	    rop[46] - rop[47] - rop[48] - rop[49] - rop[50] + rop[51] + rop[
	    53] + rop[54] + rop[58] + rop[66];
    wdot[6] = rop[24] + rop[25] + rop[37] + rop[39] + rop[40] + rop[41] + rop[
	    44] + rop[45] + rop[46] + rop[48] + rop[49] + rop[50] + rop[59];
    wdot[7] = -rop[3] + rop[14] + rop[15] - rop[25] - rop[26] - rop[27] - rop[
	    41] - rop[51] - rop[52] - rop[53] + rop[71] + rop[72];
    wdot[8] = -rop[6] + rop[7] + rop[28] - rop[29] + rop[30] + rop[36] - rop[
	    44] - rop[45] + rop[46] - rop[52] + rop[55] + rop[56] + rop[56] + 
	    rop[60] + rop[62] + rop[62] - rop[66] - rop[67] - rop[68];
    wdot[9] = -rop[7] + rop[29] - rop[30] - rop[46] + rop[52] - rop[56] - rop[
	    62] + rop[67] + rop[68];
    wdot[10] = rop[5] - rop[8] + rop[9] - rop[13] + rop[32] - rop[38] - rop[
	    47] + rop[48] - rop[53] + rop[58] + rop[59] + rop[65] + rop[67] + 
	    rop[69] + rop[71];
    wdot[11] = rop[8] + rop[10] + rop[13] + rop[47] + rop[53] - rop[65];
    wdot[12] = rop[6] - rop[11] + rop[12] - rop[14] + rop[31] - rop[33] - rop[
	    34] + rop[35] + rop[38] + rop[42] + rop[43] - rop[49] + rop[50] + 
	    rop[51] + rop[65] + rop[66] - rop[68] + rop[72];
    wdot[13] = 0.f;

    return ;//0;
} /* ckwyp_ */

