/**
 * @file rad_opthin.h
 * @brief Header file for class radiation
 */

#pragma once

#include "radiation.h"
#include <vector>
#include <string>

class domain;

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing radiation models: optically thin or two flux.
 *  Assumes gray gases
 *
 */

class rad_opthin : virtual public radiation {

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

        rad_opthin(domain *p_domn) : radiation(p_domn) {};   // constructor

        virtual ~rad_opthin(){}

};

