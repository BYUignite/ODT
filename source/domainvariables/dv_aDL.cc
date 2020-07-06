/**
 * @file dv_aDL.cc
 * @brief Source file for class dv_aDL
 */

#include "dv_aDL.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_aDL  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_aDL::dv_aDL(domain     *line,
               const      string s,
               const bool Lt,
               const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! Set aDL
 *  @param ipt \input optional point to compute at
 */

void dv_aDL::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd,0.0);
}

