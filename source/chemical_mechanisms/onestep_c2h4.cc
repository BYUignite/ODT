/**
 * @file onestep_c2h4.cc
 * @brief Source file for one step ethylene mechanism
 */

#include "onestep_c2h4.h"

////////////////////////////////////////////////////////////////////////////////
/*! getProblemSpecificRR function
 *
 */
void onestep_c2h4::getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const {

    int iO2   = 0;
    int iC2H4 = 1;
    int iCO2  = 2;
    int iH2O  = 3;
    int iN2   = 4;

    double rr_f = 2.0E12 * exp(-30000.0/(temp*1.9872155832));

    double cO2   = rho*yi[0]/31998.8;         // mol/cm3
    double cC2H4 = rho*yi[1]/28054.18;        // mol/cm3

    double rr = rr_f *
                pow(fabs(cO2),   1.65) *
                pow(fabs(cC2H4), 0.1+0.9*exp(-800*yi[iC2H4]));

    rrsp[iO2]   = -3000.0*rr;                        // kmol/m3*s
    rrsp[iC2H4] = -1000.0*rr;
    rrsp[iH2O]  = 2000.0*rr;
    rrsp[iCO2]  = 2000.0*rr;
    rrsp[iN2]   = 0.0;

}