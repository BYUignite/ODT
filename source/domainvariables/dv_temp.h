/**
 * @file dv_temp.h
 * Header file for class dv_temp
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_temp of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_temp : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_temp(){}
        dv_temp(domain     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_temp(){}

};


////////////////////////////////////////////////////////////////////////////////

