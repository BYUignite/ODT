/**
 * @file domaincase_odt_jetFlame.h
 * Header file for class domaincase_odt_jetFlame
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_jetFlame of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_jetFlame : public domaincase {


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

        domaincase_odt_jetFlame(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_jetFlame(){}

};


////////////////////////////////////////////////////////////////////////////////

