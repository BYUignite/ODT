/**
 * @file dv_ygas_flmlt.cc
 * Header file for class dv_ygas_flmlt
 */


#include "dv_ygas_flmlt.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int dv_ygas_flmlt::nspc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/*! dv_ygas  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_ygas_flmlt::dv_ygas_flmlt(domain  *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(domn->ngrd, 0.0);

    nspc = domn->gas->nSpecies();

    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.
    kMe                = domn->gas->speciesIndex(spName);

}

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 * Gas temperature needs to be set to use problem specific RR
 */

void dv_ygas_flmlt::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    static vector<vector<double> > rrSpc(nspc);    // [nspc][ngrd]
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

    if(kMe==0) {                 // to save cost, compute needed terms for all dv_ygas_flmlt objects using this one.

        for(int k=0; k<nspc; k++)
            rrSpc.at(k).resize(domn->ngrd);

        for(int i=iS; i<=iE; i++) {
#ifdef PROBLEMSPECIFICRR
            // make sure rho and T are set first (though it should be for the diffuser at least).
            for(int k=0; k<nspc; k++)
                yi.at(k) = domn->ysp[k]->d.at(i);
            getProblemSpecificRR(domn->rho->d.at(i), domn->temp->d.at(i), domn->pram->pres, &yi.at(0), &rr.at(0));
#else
            domn->domc->setGasStateAtPt(i);
            domn->gas->getNetProductionRates(&rr.at(0));
#endif
            for (int k=0; k<nspc; k++)
                rrSpc.at(k).at(i) = rr.at(k) * domn->gas->molecularWeight(k) / domn->rho->d.at(i);   // kmol/(m³ s)*(kg/kmol)*(kg/m3) = 1/s
        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSpc.at(kMe).at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 *
 * We are using a finite difference formulation here (different that the "ODT" formulation).
 */

void dv_ygas_flmlt::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(domn->ngrd, 0.0);

    //------------------ Compute the mixing term

    double dp, dm;
    double d2YdZ2;
    int i;

    //-------- left point

    i = 0;
    dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
    dm = domn->pos->d.at(i)   - 0.0;
    d2YdZ2 = 2.0/(dp+dm)*( (d.at(i+1)-d.at(i)             )/dp -
                           (d.at(i)  -domn->strm->y0[kMe] )/dm );
    rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2YdZ2;

    //-------- interior points

    for(i=1; i<domn->ngrd-1; i++) {
        dp = domn->pos->d.at(i+1) - domn->pos->d.at(i);
        dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
        d2YdZ2 = 2.0/(dp+dm)*( (d.at(i+1) - d.at(i)  )/dp -
                               (d.at(i)   - d.at(i-1))/dm );
        rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2YdZ2;
    }

    //-------- right point

    i = domn->ngrd-1;
    dp = 1.0                  - domn->pos->d.at(i);
    dm = domn->pos->d.at(i)   - domn->pos->d.at(i-1);
    d2YdZ2 = 2.0/(dp+dm)*( (domn->strm->y1[kMe] - d.at(i)  )/dp -
                           (d.at(i)             - d.at(i-1))/dm );
    rhsMix.at(i) = domn->chi->d.at(i)/2.0 * d2YdZ2;

}

