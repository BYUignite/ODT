/**
 * @file domaincase_hips_simpleRxn.h
 * Header file for class domaincase_hips_simpleRxn
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

class domaincase_hips_simpleRxn : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_hips_simpleRxn(){}
        virtual void init(domain *p_domn);
        ~domaincase_hips_simpleRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

