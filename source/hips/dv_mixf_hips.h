/**
 * @file dv_mixf_hips.h
 * Header file for class dv_mixf_hips
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_mixf_hips of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_mixf_hips : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_mixf_hips(){}
        dv_mixf_hips(domain      *line,
                     const string s,
                     const bool   Lt,
                     const bool   Lo=true);

        virtual ~dv_mixf_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

