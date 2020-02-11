/**
 * @file dv_chi_flmlt.h
 * Header file for class dv_chi_flmlt
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_chi_flmlt of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_chi_flmlt : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_chi_flmlt(){}
        dv_chi_flmlt(domain      *line,
               const string       s,
               const bool         Lt,
               const bool         Lo=true);

        virtual ~dv_chi_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

