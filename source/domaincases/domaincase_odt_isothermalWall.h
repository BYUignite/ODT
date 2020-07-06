/**
 * @file domaincase_odt_isothermalWall.h
 * @brief Header file for class domaincase_odt_isothermalWall
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_odt_isothermalWall of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_odt_isothermalWall : public domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setCaseSpecificVars();
        virtual void setGasStateAtPt(const int &ipt);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_odt_isothermalWall(){}
        virtual void init(domain *p_domn);
        ~domaincase_odt_isothermalWall(){}

};


////////////////////////////////////////////////////////////////////////////////

