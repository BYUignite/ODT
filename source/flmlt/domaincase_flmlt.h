/**
 * @file domaincase_flmlt.h
 * Header file for class domaincase_flmlt
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_flmlt of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_flmlt : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        bool LisFlameD;


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_flmlt(){}
        virtual void init(domain *p_domn);
        ~domaincase_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

