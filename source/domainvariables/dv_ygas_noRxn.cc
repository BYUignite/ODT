/**
 * @file dv_ygas_noRxn.cc
 * Header file for class dv_ygas_noRxn
 */


#include "dv_ygas_noRxn.h"
#include "domain.h"


////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function.
 *  @param ipt \input optional point to compute source at.
 */

void dv_ygas_noRxn::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc = vector<double>(domn->ngrd, 0.0);

}

