/**
 * @file domaincase.h
 * Header file for class domaincase
 */

#pragma once

#include <vector>

class domain;

using std::vector;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing base domaincase object.
 *  Derived domaincase objects will be made.
 *
 *  @author David O. Lignell
 */

class domaincase {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain                        *domn;                  ///< pointer to domain object (parent)

        vector<double>                 inlet_cell_dv_props;   ///< list of all dv properties for inserted inlet cell for channel suction/blowing case

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setGasStateAtPt(const int &ipt){}
        virtual void setCaseSpecificVars(){}
        virtual void setCaseSpecificVars_cvode(const int &ipt){}
        virtual void enforceMassFractions();
        virtual void enforceSootMom();

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        domaincase(){}
        virtual void init(domain *p_domn) { domn = p_domn; }
        virtual ~domaincase(){}

};

////////////////////////////////////////////////////////////////////////////////

