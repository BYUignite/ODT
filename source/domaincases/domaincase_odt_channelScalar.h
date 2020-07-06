/**
 * @file domaincase_odt_channelScalar.h
 * @brief Header file for class domaincase_odt_channelScalar
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_channelScalar of parent domaincase object.
 *
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class domaincase_odt_channelScalar : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_channelScalar(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_channelScalar(){}

};

////////////////////////////////////////////////////////////////////////////////

