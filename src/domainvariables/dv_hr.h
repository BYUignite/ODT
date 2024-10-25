/**
 * @file dv_hr.h
 * @brief Header file for class dv_hr
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_hr of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_hr : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_hr(){}
        dv_hr(domain       *line,
              const string s,
              const bool   Lt,
              const bool   Lo=true);

        virtual ~dv_hr(){}

};


////////////////////////////////////////////////////////////////////////////////

