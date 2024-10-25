/**
 * @file domaincase_odt_RT.h
 * @brief Header file for class domaincase_odt_RT
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_RT of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_RT : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_RT(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_RT(){}

};


////////////////////////////////////////////////////////////////////////////////

