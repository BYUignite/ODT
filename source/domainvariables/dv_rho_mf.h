/**
 * @file dv_rho_mf.h
 * @brief Header file for class dv_rho_mf
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_rho_mf of parent lv object.
 *  Density is based on simple mixing between two streams with densities rho0, rho1
 *  (Assumes ideal gas mixing with const MW, cp.)
 *
 *  @author David O. Lignell
 */

class dv_rho_mf : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        double rho0;         ///< read from input file (streams section)
        double rho1;         ///< read from input file (streams section)
        double temp0;        ///< read from input file (streams section)
        double temp1;        ///< read from input file (streams section)
        double tempFlame;    ///< read from input file (streams section)
        double Zst;          ///< read from input file (streams section)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void setVar(const int ipt=-1);

        virtual void merge2cells(const int    imrg,
                                 const double m2,
                                 const double m1,
                                 const bool   LconstVolume=false);

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_rho_mf(){}
        dv_rho_mf(domain      *line,
                  const string s,
                  const bool   Lt,
                  const bool   Lo=true);

        virtual ~dv_rho_mf(){}

};


////////////////////////////////////////////////////////////////////////////////

