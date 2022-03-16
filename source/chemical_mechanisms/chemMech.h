/**
 * @file chemMech.h
 * @brief Header file for class \ref chemMech
 */

#ifndef ODT_CHEMMECH_H
#define ODT_CHEMMECH_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing chemMech object
 *
 */

class chemMech {

    //////////////////// DATA MEMBERS //////////////////////

public:

    domain         *domn;          ///< pointer to domain object

    bool canteraRR = true;

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    virtual void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rrsp) const = 0;

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////

public:

    chemMech(domain *p_domn){ domn = p_domn; }
    virtual ~chemMech() = default;
};

#endif //ODT_CHEMMECH_H