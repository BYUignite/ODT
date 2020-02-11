/**
 * @file dv_dvisc.cc
 * Header file for class dv_dvisc
 */


#include "dv_dvisc.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/*! dv_dvisc  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_dvisc::dv_dvisc(domain    *line,
                   const      string s,
                   const bool Lt,
                   const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, domn->pram->kvisc0 * domn->pram->rho0);

    if(Lt){
        *domn->io->ostrm << endl << "WARNING, you set dvisc to be transported. Resetting L_transported to false" << endl;
        L_transported = false;
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_dvisc setVar function
 *  @param ipt \input optional point to compute at
 */

void dv_dvisc::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd, domn->pram->kvisc0 * domn->pram->rho0);
    for(int i=0; i<domn->ngrd; i++) {
        domn->domc->setGasStateAtPt(i);
        d.at(i) = domn->tran->viscosity();
    }
}

