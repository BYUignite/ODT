/**
 * @file dv_enth_hips.cc
 * Header file for class dv_enth_hips
 */

#include "dv_enth_hips.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/** dv_enth_hips  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_enth_hips::dv_enth_hips(domain      *line,
                           const string s,
                           const bool   Lt,
                           const bool   Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    rhsSrc.resize(domn->ngrd, 0.0);
    rhsMix.resize(domn->ngrd, 0.0);

    ScHips = 1.0;

}

////////////////////////////////////////////////////////////////////////////////
/** dv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *      This parameter is for the inerhited interface.
 *  This function just returns since the enthalpy source is zero in hips
 *     (well, as currently defined).
 */

void dv_enth_hips::getRhsSrc(const int ipt){
        return;
}

////////////////////////////////////////////////////////////////////////////////
/** dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *      These parameters are for the inherited interface (not used here).
 * Solving: dA/dt = (B-A)/τ,
 *          dB/dt = (A-B)/τ,
 * where A is enthalpy of one parcel, and B is enthalpy of its neighbor.
 * The solution is A = A_0 + (B_0-A_0)/2*(1-exp(-2*t/τ)),
 *                 B = B_0 + (A_0-B_0)/2*(1-exp(-2*t/τ)).
 * Note, we only need to solve for A at every point as long as we reference B.
 * Since we know the analytic solution, form the rate for use with Explicit Euler:
 * Rate_A = (A(t+Δt) - A(t)) / Δt = (B_0-A_0)/(2*Δt)*(1-exp(-2Δt/τ)).
 * Hence, Explicit Euler will give the exact solution.
 */

void dv_enth_hips::getRhsMix(const vector<double> &gf,
                             const vector<double> &dxc){

    if(!L_transported || domn->pram->LsimpleMix)
        return;

    int ime;
    int inb;

    for(int i=0; i<domn->ngrd; i++) {

        ime = domn->solv->pLoc[i];
        inb = (i%2 == 0) ? domn->solv->pLoc[i+1] : domn->solv->pLoc[i-1];

        rhsMix[ime] = (d.at(inb)-d.at(ime))/(2.0*domn->mimx->dt) *
                      (1.0 - exp(-2.0*domn->mimx->dt/domn->solv->tMix));
    }

}

