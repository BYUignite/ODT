/**
 * @file dv_enth_flmltX.h
 * Header file for class dv_enth_flmltX
 */

#pragma once

#include "dv_enth.h"
#include <string>
#include <vector>
#include "interp_linear.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_enth_flmltX of parent dv_enth object.
 *  Inherit everything, but edit the setVar for use with heat loss.
 *
 *  @author David O. Lignell
 */

class dv_enth_flmltX : public dv_enth {


    //////////////////// DATA MEMBERS     /////////////////

    public:
        Linear_interp Linterp_a;
        Linear_interp Linterp_s;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);
        virtual void set_hsens(vector<double> &hsens);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_enth_flmltX(){}
        dv_enth_flmltX(domain     *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~dv_enth_flmltX(){ }

};

////////////////////////////////////////////////////////////////////////////////

