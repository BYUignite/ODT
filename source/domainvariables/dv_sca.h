/**
 * @file dv_sca.h
 * @brief Header file for class dv_sca
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_sca of parent lv object.
 *
 *  @author David O. Lignell
 *  @author Marten Klein
 */

class dv_sca : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        double constantSource;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        //virtual void setVar(const int ipt=-1);
        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_sca(){}
        dv_sca(domain     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_sca(){}

};

////////////////////////////////////////////////////////////////////////////////

