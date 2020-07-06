/**
 * @file rad_fvdom.cc
 * @brief Source file for class rad_fvdom
 */

#include "rad_fvdom.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
/** Constructor function
**/

rad_fvdom::rad_fvdom(domain *p_domn) : radiation(p_domn) {

    if(domn->pram->cCoord == 3){
        cout << "\n#********** Error: rad_fvdom defined only for cCoord = 1 or 2" << endl;
        exit(0);
    }

    npsi = domn->pram->npsi;
    psi.resize(npsi, 0.0);
    dpsi = M_PI/npsi;
    for(int i=1; i<npsi; i++)
        psi[i] = psi[i-1] + dpsi;
    psi[npsi-1] = M_PI;

    ntheta = domn->pram->ntheta;
    if(domn->pram->cCoord==1 && ntheta != 1) {
        cout << "\n#********** WARNING: planar config needs ntheta=1 in FVDOM, resetting local value to 1";
        ntheta = 1;
    }
    theta.resize(ntheta, 0.0);
    dtheta = M_PI/2/ntheta;
    for(int i=1; i<ntheta; i++)
        theta[i] = theta[i-1] + dtheta;
    theta[ntheta-1] = M_PI/2;

    alpha = vector<vector<double> >(ntheta, vector<double>(npsi));
    beta  = vector<vector<double> >(ntheta, vector<double>(npsi));
    for(int i=0; i<ntheta; i++){
        for(int j=0; j<npsi; j++){
            if(domn->pram->cCoord == 1) {
                alpha[i][j] = cos(psi[j]);
                beta[i][j]  = sin(psi[j]);
            }
            else {    // cCoord == 2
                alpha[i][j] = sin(theta[i])*cos(psi[j]);
                beta[i][j]  = sin(theta[i])*sin(psi[j]);
            }
        }
    }

    IinfHi = sigmaSB * pow(domn->pram->TBChi,4.0)/M_PI;
    IinfLo = sigmaSB * pow(domn->pram->TBClo,4.0)/M_PI;


}
///////////////////////////////////////////////////////////////////////////////
/** Function computes radiatative heat source (W/kg = W/m3 / rho).
 *  Use finite volume discrete ordinates model
 *  @param xMoleSp   \input species mole fractions
 *  @param temp     \input temperature profile
 *  @param radSource \output gas radiation volumetric heat source (W/kg)
 */
void rad_fvdom::getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                 const vector<double>          &temp,
                                 const double                  &pressure,
                                 vector<double>                &radSource){

    vector<vector<double> >          Kabs(1, vector<double>(domn->ngrd));
    vector<vector<double> >          a(   1, vector<double>(domn->ngrd, 1.0));
    vector<vector<vector<double> > > I(domn->ngrd, vector<vector<double> >(ntheta, vector<double>(npsi))); // total intensity

    //-------------- Planck mean

    if(domn->pram->radCoefType == "PLANKMEAN") {
        radProps->get_planck_mean_coefs(xMoleSp, temp, pressure, Kabs[0]);
        if(domn->pram->cCoord == 1){
            get_I_planar(temp, Kabs[0], a[0], false, I);
            get_radSource_planar(I, radSource);
        }
        else{    // cCoord == 2
            get_I_cylindrical(temp, Kabs[0], a[0], false, I);
            get_radSource_cylindrical(I, radSource);
        }
    }
    //-------------- Spectral
    else {
        Kabs.resize(radProps->nGG, vector<double>(domn->ngrd));
        a.resize(   radProps->nGG, vector<double>(domn->ngrd));
        if (domn->pram->radCoefType == "WSGG")
            radProps->get_WSGG_coefs(xMoleSp, temp, pressure, Kabs, a);
        vector<vector<vector<double> > > Ik(domn->ngrd, vector<vector<double> >(ntheta, vector<double>(npsi))); // spectral intensity
        for(int k=0; k < radProps->nGG; k++){
            if(domn->pram->cCoord == 1)
                get_I_planar(temp, Kabs[k], a[k], k==0, Ik);
            else       // cCoord == 2
                get_I_cylindrical(temp, Kabs[k], a[k], k==0, Ik);
            for(int i=0; i<domn->ngrd; i++)
                for(int j=0; j<npsi; j++)
                    I[i][0][j] += Ik[i][0][j];
        }
        if(domn->pram->cCoord == 1)
            get_radSource_planar(I, radSource);
        else         // cCoord == 2
            get_radSource_cylindrical(I, radSource);
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes Radiative Intensity (W/m2*st).
 *  Use finite volume discrete ordinates model for planar configuration.
 *  @param temp     \input temperature profile
 *  @param Kabs     \input absorption coefficient
 *  @param a        \input weight factor
 *  @param I        \output intensity field
 */
void rad_fvdom::get_I_planar(const vector<double>             &temp,
                             const vector<double>             &Kabs,
                             const vector<double>             &a,
                             const bool                       LisClearGas,
                             vector<vector<vector<double> > > &I) {


    aIb.resize(domn->ngrd);
    for(int i=0; i<domn->ngrd; i++)
        aIb[i] = sigmaSB * pow(temp[i],4.0)/M_PI * a[i];

    domn->mesher->setGridDx(domn, V);

    int i,k;

    //-----------------------------------------

    if(LisClearGas) {
        for(k=0; k<domn->ngrd; k++){
            for(i=0; i<npsi/2; i++)
                I[k][0][i] = a[k]*IinfLo;
            for(i=npsi-1; i>npsi/2-1; i--)
                I[k][0][i] = a[k]*IinfHi;
        }
    }

    //-----------------------------------------

    else {


        for(i=0; i<npsi/2; i++){

            k = 0;
            I[k][0][i] = (V[k]*Kabs[k]*aIb[k] + alpha[0][i]*IinfLo) / (alpha[0][i] + V[k]*Kabs[k]);

            for(k=1; k<domn->ngrd; k++)
                I[k][0][i] = (V[k]*Kabs[k]*aIb[k] + alpha[0][i]*I[k-1][0][i]) / (alpha[0][i] + V[k]*Kabs[k]);
        }

        for(i=npsi-1; i>npsi/2-1; i--){

            k = domn->ngrd-1;
            I[k][0][i] = ( V[k]*Kabs[k]*aIb[k] - alpha[0][i]*IinfHi) / (-alpha[0][i] + V[k]*Kabs[k] );

            for(k=domn->ngrd-2; k>-1; k--)
                I[k][0][i] = ( V[k]*Kabs[k]*aIb[k] - alpha[0][i]*I[k+1][0][i] ) / (-alpha[0][i] + V[k]*Kabs[k] );
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Function radiatative heat source (W/kg = (W/m3) / rho).
 *  Use finite volume discrete ordinates model for planar configuration.
 *  Use Planck Mean absorption coefficients
 *  @param I         \input intensity field
 *  @param radSource \output heat source
 */

void rad_fvdom::get_radSource_planar(const vector<vector<vector<double> > > &I,
                                     vector<double>                         &radSource) {

    int & nx = domn->ngrd;
    vector<double> &x = domn->pos->d;
    q.resize(nx);

    int i,k;

    for(k=0; k<nx; k++){
        q[k] = 0.0;
        for(i=0; i<npsi; i++)
            q[k] += 2.0*M_PI*dpsi*I[k][0][i]*alpha[0][i]*beta[0][i];
    }

    //-----------------------------------------

    radSource.resize(nx);
    radSource[0]    = -(q[1] - q[0])  / (x[1] - x[0]);
    for(i=1; i<nx-1; i++)
        radSource[i] = -(q[i+1] - q[i-1]) / (x[i+1] - x[i-1]);
    radSource[nx-1] = -(q[nx-1] - q[nx-2]) / (x[nx-1] - x[nx-2]);    // W/m^3

    for(i=0; i<nx; i++)
        radSource[i] /= domn->rho->d[i];                             // W/kg
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes Radiative Intensity (W/m2*st).
 *  Use finite volume discrete ordinates model for cylindrical configuration.
 *  @param temp     \input temperature profile
 *  @param Kabs     \input absorption coefficient
 *  @param a        \input weight factor
 *  @param I        \output intensity field
 */
void rad_fvdom::get_I_cylindrical(const vector<double>             &temp,
                                  const vector<double>             &Kabs,
                                  const vector<double>             &a,
                                  const bool                       LisClearGas,
                                  vector<vector<vector<double> > > &I) {

    int            &nx = domn->ngrd;
    vector<double> &x  = domn->pos->d;
    vector<double> &xf = domn->posf->d;

    aIb.resize(domn->ngrd);
    for(int i=0; i<domn->ngrd; i++)
        aIb[i] = sigmaSB * pow(temp[i],4.0)/M_PI * a[i];

    domn->mesher->setGridDxc(domn, V, domn->pram->cCoord);
    for(int i=0; i<V.size(); i++)
        V[i] *= 0.5;
    domn->mesher->setGridDx(domn, Asn);

    int i,j,k;

    //-----------------------------------------

    if(LisClearGas) {
        for(int k=0; k<nx; k++)
            for(i=0; i<ntheta; i++) {
                for(j=0; j<npsi/2; j++)
                    I[k][i][j] = a[k]*IinfLo;
                for(j=npsi/2; j<npsi; j++)
                    I[k][j][i] = a[k]*IinfHi;
            }
    }

    ////////////////////////////////////////////////////////////////////////////////

    else{

        for(i=0; i<ntheta; i++) {

            //AAAAAAAAAAAAAAAAAA RHS, march <---- from BC, ψ goes clockwise from ----> dir

            j = npsi-1;

            k = nx-1;
            I[k][i][j] = (V[k]*Kabs[k]*aIb[k] - alpha[i][j]*abs(xf[k+1])*IinfHi) /
                       (-alpha[i][j]*abs(xf[k]) - alpha[i][j]*Asn[k] + V[k]*Kabs[k]);

            for(k=nx-2; k>nx/2; k--)
                I[k][i][j] = (V[k]*Kabs[k]*aIb[k] - alpha[i][j]*abs(xf[k+1])*I[k+1][i][j]) /
                           (-alpha[i][j]*abs(xf[k]) - alpha[i][j]*Asn[k] + V[k]*Kabs[k]);

            // Other j

            for(j=npsi-2; j>npsi/2-1; j--) {

                k = nx-1;
                I[k][i][j] = (V[k]*Kabs[k]*aIb[k] - alpha[i][j]*abs(xf[k+1])*IinfHi           + beta[i][j]*Asn[k]/dpsi*I[k][i][j+1]) /
                           (-alpha[i][j]*abs(xf[k]) - alpha[i][j]*Asn[k] + V[k]*Kabs[k]         + beta[i][j]*Asn[k]/dpsi);

                for(k=nx-2; k>nx/2; k--)
                    I[k][i][j] = (V[k]*Kabs[k]*aIb[k] - alpha[i][j]*abs(xf[k+1])*I[k+1][i][j]   + beta[i][j]*Asn[k]*I[k][i][j+1]/dpsi) /
                               (-alpha[i][j]*abs(xf[k]) - alpha[i][j]*Asn[k] + V[k]*Kabs[k]     + beta[i][j]*Asn[k]/dpsi);
            }

            //BBBBBBBBBBBBBBBBBBB LHS, march ----> from BC, ψ goes counterclockwise from ----> dir

            j = 0;

            k = 0;
            I[k][i][j] = (V[k]*Kabs[k]*aIb[k] + alpha[i][j]*abs(xf[k])*IinfLo) / \
                       ( alpha[i][j]*abs(xf[k+1]) + alpha[i][j]*Asn[k] + V[k]*Kabs[k]);

            for(k=1; k<nx/2; k++)
                I[k][i][j] = (V[k]*Kabs[k]*aIb[k] + alpha[i][j]*abs(xf[k])*I[k-1][i][j]) / \
                           ( alpha[i][j]*abs(xf[k+1]) + alpha[i][j]*Asn[k] + V[k]*Kabs[k]);

            // Other j

            for(j=1; j<npsi/2; j++) {

                k = 0;
                I[k][i][j] = (V[k]*Kabs[k]*aIb[k] + alpha[i][j]*abs(xf[k])*IinfLo           + beta[i][j]*Asn[k]/dpsi*I[k][i][j-1]) /
                           ( alpha[i][j]*abs(xf[k+1]) + alpha[i][j]*Asn[k] + V[k]*Kabs[k]     + beta[i][j]*Asn[k]/dpsi);

                for(k=1; k<nx/2; k++)
                    I[k][i][j] = (V[k]*Kabs[k]*aIb[k] + alpha[i][j]*abs(xf[k])*I[k-1][i][j]       + beta[i][j]*Asn[k]/dpsi*I[k][i][j-1]) /
                               ( alpha[i][j]*abs(xf[k+1]) + alpha[i][j]*Asn[k] + V[k]*Kabs[k]     + beta[i][j]*Asn[k]/dpsi);
            }

            //                                                                         ^
            //                                                                        /
            //CCCCCCCCCCCCCCCCCC RHS, march ----> from center, ψ goes clockwise from / dir

            for(j=npsi/2-1; j>=0; j--){

                k = nx/2;
                I[k][i][j] = I[k-1][i][j]; //+ (x[k]-x[k-1])*(I[k-1][i][j]-I[k-2][i][j])/(x[k-1]-x[k-2]); // extrapolation

                for(k=nx/2+1; k<nx; k++)
                    I[k][i][j] = (V[k]*Kabs[k]*aIb[k] + alpha[i][j]*abs(xf[k])*I[k-1][i][j]       + beta[i][j]*Asn[k]/dpsi*I[k][i][j+1]) /
                               ( alpha[i][j]*abs(xf[k+1]) - alpha[i][j]*Asn[k] + V[k]*Kabs[k]     + beta[i][j]*Asn[k]/dpsi);
            }

            //                                                                            ^
            //                                                                             |
            //DDDDDDDDDDDDDDDDD LHS, march <---- from center, ψ goes counterclockwise from \ dir

            for(j=npsi/2; j<npsi; j++){

                k = nx/2;
                I[k][i][j] = I[k+1][i][j]; //+ (x[k]-x[k+1])*(I[k+1][i][j]-I[k+2][i][j])/(x[k+1]-x[k+2]); // extrapolation

                for(k=nx/2-1; k>=0; k--)
                    I[k][i][j] = (V[k]*Kabs[k]*aIb[k] - alpha[i][j]*abs(xf[k+1])*I[k+1][i][j]   + beta[i][j]*Asn[k]/dpsi*I[k][i][j-1]) /
                               (-alpha[i][j]*abs(xf[k]) + alpha[i][j]*Asn[k] + V[k]*Kabs[k]     + beta[i][j]*Asn[k]/dpsi);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Function radiatative heat source (W/kg = (W/m3) / rho).
 *  Use finite volume discrete ordinates model for cylindrical configuration.
 *  Use Planck Mean absorption coefficients
 *  @param I         \input intensity field
 *  @param radSource \output heat source
 */

void rad_fvdom::get_radSource_cylindrical(const vector<vector<vector<double> > > &I,
                                          vector<double>                         &radSource) {

    int            &nx = domn->ngrd;
    vector<double> &x  = domn->pos->d;

    int i,j,k;

    q.resize(nx);
    for(k=0; k<nx; k++){
        q[k] = 0.0;
        for(i=0; i<ntheta; i++)
            for(j=0; j<npsi; j++)
                q[k] -= 4.0*dtheta*dpsi*I[k][i][j]*alpha[i][j]*beta[i][j];
    }

    //-----------------------------------------

    radSource.resize(nx);
    radSource[0]    = -(q[1] - q[0])  / (x[1] - x[0]);
    for(i=1; i<nx-1; i++)
        radSource[i] = -(q[i+1] - q[i-1]) / (x[i+1] - x[i-1]);
    radSource[nx-1] = -(q[nx-1] - q[nx-2]) / (x[nx-1] - x[nx-2]);    // W/m^3

    for(i=0; i<nx; i++)
        radSource[i] /= domn->rho->d[i];                             // W/kg

}
