/**
 * @file dv_enth_hips.h
 * Header file for class dv_enth_hips
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_enth_hips of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_enth_hips : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_enth_hips(){}
        dv_enth_hips(domain      *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~dv_enth_hips(){ }

};


////////////////////////////////////////////////////////////////////////////////

