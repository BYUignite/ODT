/**
 * @file domaincase_odt_channel.h
 * Header file for class domaincase_odt_channel
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_channel of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_channel : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_channel(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_channel(){}

};


////////////////////////////////////////////////////////////////////////////////

