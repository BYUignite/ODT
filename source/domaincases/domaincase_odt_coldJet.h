/**
 * @file domaincase_odt_coldJet.h
 * Header file for class domaincase_odt_coldJet
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_coldJet of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_coldJet : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_coldJet(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_coldJet(){}

};


////////////////////////////////////////////////////////////////////////////////

