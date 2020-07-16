/*
 * @file simple_dlr.cc
 * @brief Source file for two-step methane/hydrogen/air chemical mechanism
 */

using namespace std;
#include <vector>
#include <cmath>
#include <iostream>
#include<stdlib.h>

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) {
/** inputs:
 * rho [=] kg/m3
 * temp [=] K
 * yi are species mass fractions
 * outputs:
 * rrsp [=] kmol/m3*s
 */

//    CH4 +   2 O2 => CO2 + 2 H2O
//    H2  + 1/2 O2 => H2O

    int iH2   = 0;
    int iO2   = 1;
    int iH2O  = 2;
    int iCH4  = 3;
    int iCO2  = 4;
    int iN2   = 5;

    double dt = 6E-5;
    double k1 = 1000.0/dt;
    double k2 = 1000.0/dt;

    double cO2   = rho*yi[iO2] /31.9988;         // kmol/m3
    double cCH4  = rho*yi[iCH4]/16.0426;         // kmol/m3
    double cH2   = rho*yi[iH2] /2.0158;          // kmol/m3

    double r1 = k1 * cCH4 * cO2;
    double r2 = k2 * cH2  * cO2;

    rrsp[iH2]   = -r2;                        // kmol/m3*s
    rrsp[iCH4]  = -r1;
    rrsp[iO2]   = -2*r1 - 0.5*r2;
    rrsp[iH2O]  =  2*r1 + r2;
    rrsp[iCO2]  =  r1;
    rrsp[iN2]   =  0.0;

}

