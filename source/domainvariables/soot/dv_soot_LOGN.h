/**
 * @file dv_soot_LOGN.h
 * Header file for class dv_soot_LOGN
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_LOGN of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_LOGN : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:
        double M0;
        double M1;
        double M2;

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double Mk(const double &k);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_LOGN(domain      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {}

        virtual ~dv_soot_LOGN(){}

};


////////////////////////////////////////////////////////////////////////////////



