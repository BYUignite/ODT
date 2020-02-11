/**
 * @file domaincase_flmltX.h
 * Header file for class domaincase_flmltX
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_flmltX of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_flmltX : public domaincase {


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

        domaincase_flmltX(){}
        virtual void init(domain *p_domn);
        ~domaincase_flmltX(){}

};


////////////////////////////////////////////////////////////////////////////////

