/**
 * @file dv_aDL.h
 * Header file for class dv_aDL
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_mixf of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_aDL : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_aDL(){}
        dv_aDL(domain      *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~dv_aDL(){}

};


////////////////////////////////////////////////////////////////////////////////

