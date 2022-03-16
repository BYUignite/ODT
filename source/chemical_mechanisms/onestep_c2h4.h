/**
 * @file onestep_c2h4.h
 * @brief Header file for one step ethylene mechanism
 */

#ifndef ODT_ONESTEP_C2H4_H
#define ODT_ONESTEP_C2H4_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child onestep_c2h4 of parent chemMech object.
 *
 */

class onestep_c2h4 : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    /** inputs:
     * rho [=] kg/m3
     * temp [=] K
     * yi are species mass fractions
     * outputs:
     * rrsp [=] kmol/m3*s
     *
     */

    //    C2H4 + 3 O2 => 2 H2O + 2 CO2
    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    onestep_c2h4(domain *p_domn) : chemMech(p_domn) { domn = p_domn; canteraRR = false; }
    ~onestep_c2h4() override = default;

};

#endif //ODT_ONESTEP_C2H4_H