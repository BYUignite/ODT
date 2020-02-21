/**
 * @file rad_twoflux.h
 * Header file for class radiation
 */

#pragma once

#include <vector>
#include <string>
#include "radiation.h"

class domain;

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing radiation models: optically thin or two flux.
 *  Assumes gray gases
 *
 */

class rad_twoflux : virtual public radiation {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

    ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public:

        virtual void getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                      const vector<double>          &temp,
                                      const double                  &pressure,
                                      vector<double>                &radSource);
    private:

    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        rad_twoflux(domain *p_domn) : radiation(p_domn) {};   // constructor

        virtual ~rad_twoflux(){}

};

