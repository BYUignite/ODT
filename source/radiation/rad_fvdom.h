/**
 * @file rad_fvdom.h
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

class rad_fvdom : virtual public radiation {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

    private:

        int            npsi;                     ///< number of psi ray angles
        double         dpsi;                     ///< psi increment
        vector<double> psi;                      ///< psi ray angles
        int            ntheta;                   ///< number of theta ray angles
        double         dtheta;                   ///< theta increment
        vector<double> theta;                    ///< theta ray angles
        vector<double> aIb;                      ///< black intensity * weight a
        double         IinfHi;                   ///< surrounding intensity
        double         IinfLo;                   ///< surrounding intensity
        vector<double> q;                        ///< heat flux at each grid point
        vector<vector<double> > alpha;           ///< geometric coefficients
        vector<vector<double> > beta;            ///< geometric coefficients
        vector<double> V;                        ///< cell volumes
        vector<double> Asn;                      ///< south/north face areas


    ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public:

        virtual void getRadHeatSource(const vector<vector<double> > &xMoleSp,
                                      const vector<double>          &temp,
                                      const double                  &pressure,
                                      vector<double>                &radSource);


    private:

        virtual void get_I_planar(     const vector<double>             &temp,
                                       const vector<double>             &Kabs,
                                       const vector<double>             &a,
                                       const bool                       LisClearGas,
                                       vector<vector<vector<double> > > &I);

        virtual void get_I_cylindrical(const vector<double>             &temp,
                                       const vector<double>             &Kabs,
                                       const vector<double>             &a,
                                       const bool                       LisClearGas,
                                       vector<vector<vector<double> > > &I);

        virtual void get_radSource_planar(     const vector<vector<vector<double> > > &I,
                                               vector<double>                         &radSource);

        virtual void get_radSource_cylindrical(const vector<vector<vector<double> > > &I,
                                               vector<double>                         &radSource);


    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        rad_fvdom(domain *p_domn);   // constructor

        virtual ~rad_fvdom(){}

};

