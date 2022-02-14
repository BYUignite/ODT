/**
 * @file rad_opthin.cc
 * @brief Source file for class rad_opthin
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
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 */
void rad_opthin::getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                  const vector<double>          &temp,
                                  const double                  &pressure,
                                  vector<double>                &radSource){

    // warning for inconsistent boundary conditions (opthin only)
    if(domn->pram->TBClo != domn->pram->TBChi)
        cout << "\n#********** WARNING: TBClo != TBChi for opthin radiation; using TBChi" << endl;

    // make variables
    double TBChi4 = pow(domn->pram->TBChi, 4.0);

//    vector<double> Kabs(domn->ngrd);
//    radProps->get_planck_mean_coefs(xMoleSp, temp, pressure, Kabs);
//
//    for(int i=0; i<domn->ngrd; i++)
//        radSource[i] = -4*sigmaSB*Kabs[i]*(pow(temp[i],4.0) - TBChi4) / domn->rho->d[i];

    vector<double> Kabs(domn->ngrd);

    for(int i=0; i<domn->ngrd; i++) {

        double fvsoot = 0;
        if (domn->pram->Lsoot)
            fvsoot = domn->svar[1]->d[i] / domn->pram->Cmin;        // fvsoot = M1/rhoSoot = rhoGas*Ys/rhoSoot

        double xCH4 = (iRadIndx[0] < 0) ? 0 : xMoleSp.at(i).at(iRadIndx[0]);
        double xCO2 = (iRadIndx[1] < 0) ? 0 : xMoleSp.at(i).at(iRadIndx[1]);
        double xH2O = (iRadIndx[2] < 0) ? 0 : xMoleSp.at(i).at(iRadIndx[2]);
        double xCO  = (iRadIndx[3] < 0) ? 0 : xMoleSp.at(i).at(iRadIndx[3]);

        vector<double> K;       // temp storage for k values
        vector<double> A;       // temp storage for wts values

        radProps->get_k_a(K, A, temp[i], pressure, fvsoot, xH2O, xCO2, xCO, xCH4);

        for (int j=0; j<K.size(); j++)
            Kabs[i] = K[j]*A[j];

        radSource[i] = -4*sigmaSB*Kabs[i]*(pow(temp[i],4.0) - TBChi4) / domn->rho->d[i];
    }

}