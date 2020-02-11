/**
 * @file domaincase_hips_comb.h
 * Header file for class domaincase_hips_comb
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_hips_comb of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_hips_comb : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setGasStateAtPt(const int &ipt);
        virtual void setCaseSpecificVars();
        virtual void setCaseSpecificVars_cvode(const int &ipt);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_hips_comb(){}
        virtual void init(domain *p_domn);
        ~domaincase_hips_comb(){}

};


////////////////////////////////////////////////////////////////////////////////

