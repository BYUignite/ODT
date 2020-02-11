/**
 * @file dv_soot_flmlt.h
 * Header file for class dv_soot_flmlt
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_flmlt of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_flmlt : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_flmlt(domain       *line,
                      const string  s,
                      const bool    Lt,
                      const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {}

        virtual ~dv_soot_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////



