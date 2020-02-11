/**
 * @file dv_ygas_cold_hips.cc
 * Header file for class dv_ygas_cold_hips
 */


#include "dv_ygas_cold_hips.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int dv_ygas_cold_hips::nspc;

////////////////////////////////////////////////////////////////////////////////
/** dv_ygas  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_ygas_cold_hips::dv_ygas_cold_hips(domain     *line,
                           const      string s,
                           const bool Lt,
                           const bool Lo) {

    domn               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(domn->ngrd, 0.0);

    rhsSrc.resize(domn->ngrd, 0.0);
    rhsMix.resize(domn->ngrd, 0.0);

    nspc = 4;
    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.

    if(spName=="A") kMe = 0;
    if(spName=="B") kMe = 1;
    if(spName=="R") kMe = 2;
    if(spName=="P") kMe = 3;

    //Da1 = domn->io->inputFile["rxnParams"]["Da1"].as<double>();
    //Da2 = domn->io->inputFile["rxnParams"]["Da2"].as<double>();

    Da1 = domn->io->dvParams["Da1"].as<double>();
    Da2 = domn->io->dvParams["Da2"].as<double>();

    ScHips = 1.0;

}

////////////////////////////////////////////////////////////////////////////////
/** Source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  Gas temperature/rho needs to be set to use problem specific RR
 */

void dv_ygas_cold_hips::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    static vector<vector<double> > rrSpc(nspc, vector<double>(domn->ngrd));    // [nspc][ngrd]
    static vector<double>          yi(nspc);       // [nspc]
    static vector<double>          rr(nspc);       // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                 // to save cost, compute needed terms for all dv_ygas_cold_hips objects using this one.

        for(int i=iS; i<=iE; i++) {

            double A = domn->ysp[0]->d.at(i);
            double B = domn->ysp[1]->d.at(i);
            double R = domn->ysp[2]->d.at(i);
            double P = domn->ysp[3]->d.at(i);

            rrSpc.at(0).at(i) = -Da1*A*B    - 0.5*Da2*A*R;
            rrSpc.at(1).at(i) = -Da1*A*B;
            rrSpc.at(2).at(i) = 2.0*Da1*A*B - Da2*A*R;
            rrSpc.at(3).at(i) = 1.5*Da2*A*R;
        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSpc.at(kMe).at(i);
}

////////////////////////////////////////////////////////////////////////////////
/** dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 *      These parameters are for the inherited interface (not used here).
 * Solving: dA/dt = (B-A)/τ,
 *          dB/dt = (A-B)/τ,
 * where A is Y_k of one parcel, and B is Y_k of its neighbor.
 * The solution is A = A_0 + (B_0-A_0)/2*(1-exp(-2*t/τ)),
 *                 B = B_0 + (A_0-B_0)/2*(1-exp(-2*t/τ)).
 * Note, we only need to solve for A at every point as long as we reference B.
 * Since we know the analytic solution, form the rate for use with Explicit Euler:
 * Rate_A = (A(t+Δt) - A(t)) / Δt = (B_0-A_0)/(2*Δt)*(1-exp(-2Δt/τ)).
 * Hence, Explicit Euler will give the exact solution.
 */

void dv_ygas_cold_hips::getRhsMix(const vector<double> &gf,
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

