/**
 * @file dv_temp.cc
 * @brief Source file for class dv_temp
 */


#include "dv_temp.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_temp  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_temp::dv_temp(domain  *line,
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
/*! Set temperature from the gas state
 *  @param ipt \input optional point to compute at
 */

void dv_temp::setVar(const int ipt){

    d.resize(domn->ngrd);
    if(ipt == -1)
        for(int i=0; i<domn->ngrd; i++) {
            try {
                domn->domc->setGasStateAtPt(i);
            } catch (const CanteraError& e) {
                throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
            }
            d.at(i) = domn->gas->thermo()->temperature();
        }
    else {
        try {
            domn->domc->setGasStateAtPt(ipt);
        } catch (const CanteraError& e) {
            throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
        }
        d.at(ipt) = domn->gas->thermo()->temperature();
    }
}