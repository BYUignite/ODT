/**
 * @file chemNone.h
 * @brief Header file for chem object used for non-reacting ODT runs
 */

#ifndef ODT_CHEMNONE_H
#define ODT_CHEMNONE_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child chemNone of parent chemMech object.
 *
 */

class chemNone : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override { return; };

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    chemNone(domain *p_domn) : chemMech(p_domn) { domn = p_domn; canteraRR = false; }
    ~chemNone() override = default;

};

#endif //ODT_CHEMNONE_H