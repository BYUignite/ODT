/**
 * @file rad_twoflux.cc
 * Source file for class rad_twoflux
 */

#include "rad_twoflux.h"
#include "domain.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiatative heat source (W/kg = W/m3 / rho).
 *  Use optically thin for Imode == 1, twoflux for Imod == 2.
 *  @param xMoleSp   \input species mole fractions
 *  @param &temp     \input temperature profile
 *  @param fvSoot    \input optional soot volume fraction
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 */
void rad_twoflux::getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                   const vector<double>          &temp,
                                   const double                  &pressure,
                                   const vector<double>          &fvSoot,
                                   vector<double>                &radSource) {

    if(domn->pram->cCoord != 1){
        cout << "\n#********** Error: rad_twoflux defined only for cCoord = 1" << endl;
        exit(0);
    }

   int npts = domn->ngrd + 2;             // add in the two boundaries
   vector<double> qp(npts);
   vector<double> qm(npts);
   vector<double> dx(npts-1);
   vector<double> gasEmmTerm(npts);       // 2*kg*sigma*Tg^4

   //------------- Get the dx array

   dx[0] = 0.5*(domn->posf->d[1]-domn->posf->d[0]);
   for(int i=1; i<npts-2; i++)
       dx[i] = 0.5*(domn->posf->d[i+1]-domn->posf->d[i-1]);
   dx[npts-2] = 0.5*(domn->posf->d[domn->ngrdf-1]-domn->posf->d[domn->ngrdf-2]);

   //------------- Get the gas absorption coefficient

   vector<double> Kabs(domn->ngrd);
   radProps->get_planck_mean_coefs(xMoleSp, temp, pressure, fvSoot, Kabs);
   Kabs.insert(Kabs.begin(), Kabs[0]);
   Kabs.push_back(Kabs[Kabs.size()-1]);

   //------------- Get the gas emmision term

   for(int i=0; i<npts-2; i++)
       gasEmmTerm[i+1] = 2.0*Kabs[i+1] * sigmaSB*pow(temp[i],4.0);
   gasEmmTerm[0]      = gasEmmTerm[1];
   gasEmmTerm[npts-1] = gasEmmTerm[npts-2];

   //------------- Get radiative fluxes: qp, qm
   // Marching qp low to high, qm high to low

   double qpBClo = sigmaSB * pow(domn->pram->TBClo,4.0);
   double qmBChi = sigmaSB * pow(domn->pram->TBChi,4.0);

   qp[0] = qpBClo;
       for(int i=1; i<npts; i++)
       qp[i] = ( qp[i-1] + dx[i-1]*gasEmmTerm[i] ) / ( 1. + 2.*Kabs[i]*dx[i-1] );

   qm[npts-1] = qmBChi;
       for(int i=npts-2; i>=0; i--)
           qm[i] = ( qm[i+1] + dx[i]*gasEmmTerm[i] ) / ( 1. + 2.*Kabs[i]*dx[i] );

   //------------- Get gas radiative source: (div q) (=) W/kg (additive source to E-bal)

    for(int i=0, ip=1; i<npts-2; i++, ip++)
       radSource[i] = -2.*Kabs[ip] * (2.*sigmaSB*pow(temp[i],4.0) - qp[ip] - qm[ip]) /
                      domn->rho->d[i];
}

