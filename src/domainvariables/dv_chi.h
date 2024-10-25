/**
 * @file dv_chi.h
 * @brief Header file for class dv_chi
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_chi of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_chi : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_chi(){}
        dv_chi(domain      *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~dv_chi(){}

};


////////////////////////////////////////////////////////////////////////////////

