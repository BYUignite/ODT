/**
 * @file dv_rho.h
 * Header file for class dv_rho
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_rho of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_rho : public dv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_rho(){}
        dv_rho(domain      *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~dv_rho(){}

};


////////////////////////////////////////////////////////////////////////////////

