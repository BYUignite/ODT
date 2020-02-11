/**
 * @file dv_enth_flmlt.h
 * Header file for class dv_enth_flmlt
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>
#include "interp_linear.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_enth_flmlt of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_enth_flmlt : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    public:
        Linear_interp Linterp_s;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);
        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);
        virtual void set_hsens(vector<double> &hsens);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_enth_flmlt(){}
        dv_enth_flmlt(domain      *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~dv_enth_flmlt(){ }

};

////////////////////////////////////////////////////////////////////////////////

