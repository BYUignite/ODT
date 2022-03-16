/**
 * @file ch4red.h
 * @brief Header file for one step ethylene mechanism
 */

#ifndef ODT_CH4RED_H
#define ODT_CH4RED_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child ch4red of parent chemMech object.
 *
 */

class ch4red : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:


    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    ch4red(domain *p_domn) : chemMech(p_domn) { domn = p_domn; canteraRR = false; }
    ~ch4red() override = default;

};

#endif //ODT_CH4RED_H