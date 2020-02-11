/**
 * @file domaincase_hips.h
 * Header file for class domaincase_hips
 */

#pragma once

#include "domaincase.h"
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child domaincase_hips of parent domaincase object.
 *
 *  @author David O. Lignell
 */

class domaincase_hips : public domaincase {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase_hips(){}
        virtual void init(domain *p_domn);
        ~domaincase_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

