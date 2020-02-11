/**
 * @file dv_enth.h
 * Header file for class dv_enth
 */

#pragma once

#include "dv.h"
#include "radiation/radiation.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_enth of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_enth : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        int         nspc;
        int         iMe;
        radiation  *rad;
        bool        LdoSpeciesFlux;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_enth(){ rad = 0; }
        dv_enth(domain      *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_enth(){ if(rad) delete rad; }

};


////////////////////////////////////////////////////////////////////////////////

