/**
 * @file radiation.h
 * @brief Header file for class radiation
 */

#pragma once

#include <vector>
#include <string>
//#include "radiationProperties.h"
#include "rad.h"
#include "rad_planck_mean.h"
#include "rad_rcslw.h"
#include "rad_wsgg.h"

class domain;

using namespace std;

// TODO:
// add pressure to the getRadHeatSource call in enthalpy variable
// need a radiative property class object
// use constructor instead of init function here.
// verify destructor for child classes for radProps

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
 *  @author Victoria B. Stephens 
 *  @author David O. Lignell
 *
 */

class radiation {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        domain                      *domn;          ///< pointer to domain
//        radiationProperties         *radProps;      ///< tools for getting k's and a's, etc. // defines radProps
        rad                         *radProps;      ///< tools for getting k's and a's, etc. // defines radProps // radlib

        double                      sigmaSB;        ///< Stefan Boltzman const

        int                              nRadSp;     ///< number of radiating species
        vector<int>                      iRadIndx;   ///< index of radiating species
        int                              nGG;        ///< number of grey gases for WSGG approaches

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