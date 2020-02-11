/**
 * @file dv_dvisc.h
 * Header file for class dv_dvisc
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_dvisc of parent lv object.
 *  Dynamic viscosity
 *
 *  @author David O. Lignell
 */

class dv_dvisc : public dv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void setVar(const int ipt=-1);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_dvisc(){}
        dv_dvisc(domain      *line,
                 const string s,
                 const bool   Lt,
                 const bool   Lo=true);

        virtual ~dv_dvisc(){}

};


////////////////////////////////////////////////////////////////////////////////

