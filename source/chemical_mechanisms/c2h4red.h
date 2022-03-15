/**
 * @file c2h4red.h
 * @brief Header file for one step ethylene mechanism
 */

#ifndef ODT_C2H4RED_H
#define ODT_C2H4RED_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child c2h4red of parent chemMech object.
 *
 */

class c2h4red : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:


    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    c2h4red() : chemMech() {}
    ~c2h4red() override = default;

};

#endif //ODT_C2H4RED_H