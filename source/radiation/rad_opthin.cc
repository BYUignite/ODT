/**
 * @file rad_opthin.cc
 * Source file for class rad_opthin
 */

#include "rad_opthin.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiatative heat source (W/kg = W/m3 / rho).
 *  Use optically thin for Imode == 1, twoflux for Imod == 2.
 *  @param xMoleSp   \input species mole fractions
 *  @param temp      \input temperature profile
 *  @param fvSoot    \input optional soot volume fraction
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 */
void rad_opthin::getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                  const vector<double>          &temp,
                                  const double                  &pressure,
                                  const vector<double>          &fvSoot,
                                  vector<double>                &radSource){

    // warning for inconsistent boundary conditions (opthin only)
    if(domn->pram->TBClo != domn->pram->TBChi)
        cout << "\n#********** WARNING: TBClo != TBChi for opthin radiation; using TBChi" << endl;

    // make variables
    double TBChi4 = pow(domn->pram->TBChi, 4.0);

    vector<double> Kabs(domn->ngrd);
    radProps->get_planck_mean_coefs(xMoleSp, temp, pressure, fvSoot, Kabs);

    for(int i=0; i<domn->ngrd; i++)
        radSource[i] = -4*sigmaSB*Kabs[i]*(pow(temp[i],4.0) - TBChi4) / domn->rho->d[i];

}

