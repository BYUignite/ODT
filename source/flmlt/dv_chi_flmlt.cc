/**
 * @file dv_chi_flmlt.cc
 * Header file for class dv_chi_flmlt
 */


#include "dv_chi_flmlt.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_chi_flmlt  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_chi_flmlt::dv_chi_flmlt(domain  *line,
                 const               string s,
                 const bool          Lt,
                 const bool          Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/** Set scalar dissipation rate (chi) \cond
 *  @param ipt \input optional point to compute at
 *  Chi = chi0 * (1-(2*mixf - 1)^2)^2
 */

void dv_chi_flmlt::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be < 1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd);

    for(int i=0; i<domn->ngrd; i++)
        d.at(i) = domn->pram->chi0 * pow(1.0-pow(2.0*domn->pos->d.at(i)-1.0, 2.0), 2.0);

}

