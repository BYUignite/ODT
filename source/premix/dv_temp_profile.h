/**
 * @file dv_temp_profile.h
 * Header file for class dv_temp_profile
 */

#pragma once

#include "dv.h"
#include "interp_linear.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_temp_profile of parent lv object.
 *  Set a profile temperature profile to interpolate
 *  @author David O. Lignell
 */

class dv_temp_profile : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        vector<double> x;       ///< fixed temperature profile: x values
        vector<double> T;       ///< fixed temperature profile: T values: T(x)
        Linear_interp Linterp;  ///< interpolator

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_temp_profile(){}
        dv_temp_profile(domain     *line,
                const string s,
                const string p_bcType,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_temp_profile(){}
};


////////////////////////////////////////////////////////////////////////////////

