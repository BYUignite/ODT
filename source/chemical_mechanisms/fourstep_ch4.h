/**
 * @file fourstep_ch4.h
 * @brief Header file for one step ethylene mechanism
 */

#ifndef ODT_FOURSTEP_CH4_H
#define ODT_FOURSTEP_CH4_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child fourstep_ch4 of parent chemMech object.
 *
 */

class fourstep_ch4 : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:


    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    fourstep_ch4(domain *p_domn) : chemMech(p_domn) { domn = p_domn; canteraRR = false; }
    ~fourstep_ch4() override = default;

};

#endif //ODT_FOURSTEP_CH4_H