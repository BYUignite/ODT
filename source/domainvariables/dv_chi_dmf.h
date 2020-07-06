/**
 * @file dv_chi_dmf.h
 * @brief Header file for class dv_chi_dmf -- uses constant Dmf for computing chi
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_chi_dmf of parent lv object.  This uses Dmf, constant mixture fraction diffusivity, for use without cantera.
 *
 *  @author David O. Lignell, John C. Hewson
 */

class dv_chi_dmf : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:
        double Dmf;     ///<  default mixture fraction diffusivity (read from input file)

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:
        virtual void setVar(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_chi_dmf(){}
        dv_chi_dmf(domain   *line,
               const string s,
               const bool   Lt,
               const bool   Lo=true);

        virtual ~dv_chi_dmf(){}

};


////////////////////////////////////////////////////////////////////////////////

