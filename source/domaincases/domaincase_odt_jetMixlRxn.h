/**
 * @file domaincase_odt_jetMixlRxn.h
 * Header file for class domaincase_odt_jetMixlRxn
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_jetMixlRxn of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_jetMixlRxn : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_jetMixlRxn(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_jetMixlRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

