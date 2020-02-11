/**
 * @file dv_soot_QMOM.h
 * Header file for class dv_soot_QMOM
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_QMOM of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_QMOM : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double>        wts;        ///< weights from inversion algorithm
        vector<double>        absc;       ///< abscissas from inversion algoritm

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double  Mk(double exp);
        void    getWtsAbs(vector<double> M, vector<double> &wts, vector<double> &abs);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_QMOM(domain      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {
            wts.resize(nsvar/2);
            absc.resize(nsvar/2);
        }

        virtual ~dv_soot_QMOM(){}

};


////////////////////////////////////////////////////////////////////////////////



