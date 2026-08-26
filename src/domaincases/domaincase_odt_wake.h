/**
 * @file domaincase_odt_wake.h
 * @brief Header file for class domaincase_odt_wake
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_wake of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_wake : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_wake(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_wake(){}

};


////////////////////////////////////////////////////////////////////////////////

