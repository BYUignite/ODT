/**
 * @file onestep_ch4.h
 * @brief Header file for one step methane mechanism
 */

#ifndef ODT_ONESTEP_CH4_H
#define ODT_ONESTEP_CH4_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child onestep_ch4 of parent chemMech object.
 *
 */

class onestep_ch4 : public chemMech {

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

    //    CH4 + 2 O2 => 2 H2O + CO2
    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    onestep_ch4() : chemMech() {}
    ~onestep_ch4() override = default;

};

#endif //ODT_ONESTEP_CH4_H