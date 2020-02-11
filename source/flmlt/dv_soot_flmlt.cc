/**
 * @file dv_soot_flmlt.cc
 * Header file for class dv_soot_flmlt
 * @author Victoria B. Lansinger and David Lignell
 */

#include "dv_soot_flmlt.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

/*//////////////////////////////////////////////////////////////////////////////
 *! dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *
 * See page 191 of Lignell_thesis, equation 7.9.
 * Equation is dm/dt = C*dm/dZ + S.
 * This is hyperbolic and C is a "convective wind".
 * We'll upwind the discretization of dm/dZ depending on the sign of C.
 * Note, here S does not include the m_r factor and also the M_r/rho term
 *        that appears in the thesis. The M_r/rho is the reaction source term
 *        computed elsewhere. S without the m_r factor keeps S and C independent of m_r
 *        so that we can compute them once for kMe = 0. The m_r factor is then put
 *        in at the bottom.
 * As for ygas, we are just using a finite difference formulation here.
 */

void dv_soot_flmlt::getRhsMix(const vector<double> &gf,
                              const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(domn->ngrd, 0.0);

    //------------------

    static vector<double> D;
    static vector<double> mu;
    static vector<double> beta;
    static vector<double> C;
    static vector<double> S;
    double dp, dm;
    int i;

    //------------------

    if(kMe==0) {

        D.resize(domn->ngrd);
        mu.resize(domn->ngrd);
        beta.resize(domn->ngrd);
        C.resize(domn->ngrd);
        S.resize(domn->ngrd);

        //------------------

        for(i=0; i<domn->ngrd; i++){
            //domn->domc->setGasStateAtPt(i);   // keep commented for decoupled soot
            D[i] = domn->tran->thermalConductivity() / domn->rho->d[i] / domn->gas->cp_mass();
            mu[i] = domn->tran->viscosity();
            beta[i] = sqrt(domn->chi->d[i]*0.5/D[i]);
        }

        //------------------ Set C and S

        double dTdZ, drDbdZ, dmbTdZ, d2TdZ2;

        //-------- left point

        i = 0;
        dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
        dm = domn->pos->d.at(i)   - 0.0;
        dTdZ   = (domn->temp->d[i] - domn->strm->T0)/dm;
        drDbdZ = (domn->rho->d[i]*D[i]*beta[i] - 0.0)/dm;
        dmbTdZ = (mu[i]*beta[i]/domn->temp->d[i] - 0.0)/dm;
        d2TdZ2 = 2.0/(dp+dm)*( (domn->temp->d.at(i+1) - domn->temp->d.at(i))/dp -
                               (domn->temp->d.at(i)   - domn->strm->T0     )/dm );

        C[i] = 0.554*mu[i]/domn->rho->d[i]/domn->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/domn->rho->d[i]*drDbdZ;
        S[i] = 0.554/domn->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/domn->temp->d[i]*d2TdZ2 );

        //-------- interior points

        for(i=1; i<domn->ngrd-1; i++) {
            dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
            dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
            dTdZ   = (domn->temp->d[i+1] - domn->temp->d[i-1])/(dm+dp);
            drDbdZ = (domn->rho->d[i+1]*D[i+1]*beta[i+1] - domn->rho->d[i-1]*D[i-1]*beta[i-1])/(dm+dp);
            dmbTdZ = (mu[i+1]*beta[i+1]/domn->temp->d[i+1] - mu[i-1]*beta[i-1]/domn->temp->d[i-1])/(dm+dp);
            d2TdZ2 = 2.0/(dp+dm)*( (domn->temp->d[i+1] - domn->temp->d[i]  )/dp -
                                   (domn->temp->d[i]   - domn->temp->d[i-1])/dm );

            C[i] = 0.554*mu[i]/domn->rho->d[i]/domn->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/domn->rho->d[i]*drDbdZ;
            S[i] = 0.554/domn->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/domn->temp->d[i]*d2TdZ2 );
        }

        //-------- right point

        i = domn->ngrd-1;
        dp = 1.0                  - domn->pos->d.at(i);
        dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
        dTdZ   = (domn->strm->T1 - domn->temp->d[i])/dp;
        drDbdZ = (0.0 - domn->rho->d[i]*D[i]*beta[i])/dp;
        dmbTdZ = (0.0 - mu[i]*beta[i]/domn->temp->d[i])/dp;
        d2TdZ2 = 2.0/(dp+dm)*( (domn->strm->T1      - domn->temp->d.at(i)  )/dp -
                               (domn->temp->d.at(i) - domn->temp->d.at(i-1))/dm );

        C[i] = 0.554*mu[i]/domn->rho->d[i]/domn->temp->d[i]*beta[i]*beta[i]*dTdZ - beta[i]/domn->rho->d[i]*drDbdZ;
        S[i] = 0.554/domn->rho->d[i]*( beta[i]*dTdZ*dmbTdZ + beta[i]*beta[i]*mu[i]/domn->temp->d[i]*d2TdZ2 );

    }

    //------------------ Set mixing term
    // dm/dt = ( C*dm/dZ + m*S ) + RxnTerm
    // Upwind the C*dm/dZ part of the mixing term.

    double dmdZ;

    //-------- left point

    i = 0;
    dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
    dm = domn->pos->d.at(i)   - 0.0;
    dmdZ = (C[i] < 0.0) ? (d[i]-0.0)/dm : (d[i+1]-d[i])/dp;
    rhsMix[i] = C[i]*dmdZ + d[i]*S[i];              // d*S puts in the m_r factor on S.

    //-------- interior points

    for(i=1; i<domn->ngrd-1; i++) {
        dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
        dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
        dmdZ = (C[i] < 0.0) ? (d[i]-d[i-1])/dm : (d[i+1]-d[i])/dp;
        rhsMix[i] = C[i]*dmdZ + d[i]*S[i];
    }

    //-------- right point

    i = domn->ngrd-1;
    dp = 1.0                  - domn->pos->d.at(i);
    dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
    dmdZ = (C[i] < 0.0) ? (d[i]-d[i-1])/dm : (0.0-d[i])/dp;
    rhsMix[i] = C[i]*dmdZ + d[i]*S[i];

}












//void dv_soot_flmlt::getRhsMix(const vector<double> &gf,
//                        const vector<double> &dxc){
//
//    if(!L_transported) return;
//
//    rhsMix.resize(domn->ngrd, 0.0);
//
//    //------------------ Compute the mixing term
//
//    double dp, dm;
//    double d2mdZ2;
//    int i;
//
//    //-------- left point
//
//    i = 0;
//    dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
//    dm = domn->pos->d.at(i)   - 0.0;
//    d2mdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i) )/dp -
//                           (d.at(i)   - 0.0     )/dm );
//    rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2mdZ2;
//
//    //-------- interior points
//
//    for(i=1; i<domn->ngrd-1; i++) {
//        dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
//        dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
//        d2mdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i)  )/dp -
//                               (d.at(i)   - d.at(i-1))/dm );
//        rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2mdZ2;
//    }
//
//    //-------- right point
//
//    i = domn->ngrd-1;
//    dp = 1.0                  - domn->pos->d.at(i);
//    dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
//    d2mdZ2 = 2.0/(dp+dm)*( (0.0                 - d.at(i)  )/dp -
//                           (d.at(i)             - d.at(i-1))/dm );
//    rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2mdZ2;
//
//}

