/**
 * @file dv_mixf.h
 * Header file for class dv_mixf
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

class dv_mixf : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        double Dmf;     ///<  default mixture fraction diffusivity (read from input file)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);
        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_mixf(){}
        dv_mixf(domain      *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_mixf(){}

};


////////////////////////////////////////////////////////////////////////////////

