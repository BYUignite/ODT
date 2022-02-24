/**
 * @file radiation.h
 * @brief Header file for class radiation
 */

#pragma once

#include <vector>
#include <string>
#include "radPropModel.h"
#include "rad_planck_mean.h"
#include "rad_rcslw.h"
#include "rad_wsgg.h"

class domain;

using namespace std;

// TODO:
// add pressure to the getRadHeatSource call in enthalpy variable
// verify destructor for child classes for radProps
// check for missing species according to chosen radProps model

///////////////////////////////////////////////////////////////////////////////
/** Class implementing radiation models
 *
 *  Base class function for various RTE solution methods.
 *
 *  Child classes:
 *      rad_opthin.h    OPTHIN: optically thin approximation
 *      rad_twoflux.h   TWOFLUX: two-flux approximation, i.e. Schuster-Schwarzchild approximation
 *      rad_exact.h     EXACT: exact analytical solution with no scattering
 *      rad_fvdom.h     FVDOM: finite volume discrete ordinates method
 *      rad_raytrace.h  RAYTRACE: ray tracing solution
 *
 */

class radiation {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        domain          *domn;                  ///< pointer to domain
        radPropModel    *radProps;              ///< tools for getting k's and a's, etc. // defines radProps // radlib

        double          sigmaSB = 5.670E-8;     ///< Stefan Boltzman const (W/m2*K4)

        int             nRadSp;                 ///< number of radiating species
        vector<int>     iRadIndx;               ///< index of radiating species
        int             nGG;                    ///< number of grey gases for WSGG approaches

    ////////////////////// MEMBER FUNCTIONS /////////////////////

    public:

        virtual void getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                      const vector<double>          &temp,
                                      const double                  &pressure,
                                      vector<double>                &radSource) = 0;

    private:

    ////////////////////// CONSTRUCTOR FUNCTIONS /////////////////////

    public:

        radiation(domain *p_domn);   // constructor

        virtual ~radiation(){delete radProps;}          // destructor

};