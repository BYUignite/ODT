/**
 * @file domaincase_premixedFlameBurner.h
 * Header file for class domaincase_premixedFlameBurner
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_premixedFlame of parent domaincase object.
 *
 *  @author David O. Lignell
 *  @author Victoria B. Stephens
 */

class domaincase_premixedFlameBurner : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        double mixf_reactants;


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_premixedFlameBurner(){}
        virtual void init(domain *p_domn);
        ~domaincase_premixedFlameBurner(){}

};


////////////////////////////////////////////////////////////////////////////////

