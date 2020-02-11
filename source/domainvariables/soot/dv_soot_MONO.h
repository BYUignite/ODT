/**
 * @file dv_soot_MONO.h
 * Header file for class dv_soot_MONO
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_MONO of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_MONO : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double>        wts;        ///< weights of the particle size distribution
        vector<double>        absc;       ///< abscissas of the particle size distribution

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_MONO(domain      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {
            wts.resize(1);
            absc.resize(1);
        }

        virtual ~dv_soot_MONO(){}

};


////////////////////////////////////////////////////////////////////////////////



