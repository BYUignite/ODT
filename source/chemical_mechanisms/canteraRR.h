/**
 * @file canteraRR.h
 * @brief Header file for reaction mechanisms handled by Cantera
 */

#ifndef ODT_CANTERARR_H
#define ODT_CANTERARR_H

#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child canteraRR of parent chemMech object.
 *
 */

class canteraRR : public chemMech {

    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const override {
        domn->gas->kinetics()->getNetProductionRates(&rrsp[0]);
    };

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    canteraRR(domain *p_domn) : chemMech(p_domn) { domn = p_domn; }
    ~canteraRR() override = default;

};

#endif //ODT_CANTERARR_H