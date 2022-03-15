/**
 * @file simple_dlr.h
 * @brief Header file for one step methane mechanism
 */

#ifndef ODT_SIMPLE_DLR_H
#define ODT_SIMPLE_DLR_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child simple_dlr of parent chemMech object.
 *
 */

class simple_dlr : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    /** inputs:
     * rho [=] kg/m3
     * temp [=] K
     * yi are species mass fractions
     * outputs:
     * rrsp [=] kmol/m3*s
     */

    //    CH4 +   2 O2 => CO2 + 2 H2O
    //    H2  + 1/2 O2 => H2O
    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    simple_dlr() : chemMech() {}
    ~simple_dlr() override = default;

};

#endif //ODT_SIMPLE_DLR_H