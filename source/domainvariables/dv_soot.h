/**
 * @file dv_soot.h
 * Header file for class dv_soot
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

//#include "constants.h"
#include "sootModel.h"
#include "state.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot of parent dv object.
 *  This is a virtual base class.
 */

class dv_soot : public dv {

    //////////////////// DATA MEMBERS //////////////////////

protected:

    static int              N;                      ///< iterator for assigning kMe values

    int                     kMe;                    ///< kMe
    static int              nsvar;                  ///< number of soot variables

    //-----------

    static vector<double>  *yi;                     ///< pointer to species mass fractions

    //-----------

    static int   i_c2h2;                 ///< soot gas species indices
    static int   i_o2;
    static int   i_o;
    static int   i_h;
    static int   i_h2;
    static int   i_oh;
    static int   i_h2o;
    static int   i_co;
    static int   i_c10h8;       ///< soot PAH species indices
    static int   i_c12h8;
    static int   i_c12h10;
    static int   i_c14h10;
    static int   i_c16h10;
    static int   i_c18h10;
//    static int              i_elem_c;               ///< index if element C
//    static int              i_elem_h;               ///< index if element H
    static vector<int>      i_pah;                  ///< vector of PAH species indicies
//    static vector<double>   MW_sp;                  ///< vector of molecular weights
//    static double           MW_c;                   ///< mw of carbon 12.011;
//    static double           MW_h;                   ///< mw of hydrogen 1.0079;

    double                  Cmin;                   ///< number of carbons in soot nucleation size
    double                  rhoSoot;                ///< soot density kg/m3
    double                  b_coag;                 ///< coagulation rate parameter (see Lignell thesis page 58.)

    soot::nucleationMech    nucleation_mech;        ///< soot nucleation chemistry flag
    soot::growthMech        growth_mech;            ///< soot growth chemistry flag
    soot::oxidationMech     oxidation_mech;         ///< soot oxidation chemistry flag
    soot::coagulationMech   coagulation_mech;       ///< soot coagulation mechanism flag
    soot::psdMech           psd_mech;               ///< soot PSD method flag

    soot::sootModel         *SM;                     ///< pointer to soot model object (sootlib)
    soot::state             *S;                      ///< pointer to soot state object (sootlib)

    //////////////////// MEMBER FUNCTIONS /////////////////

public:

    virtual void getRhsSrc(const int ipt=-1);
    virtual void getRhsMix(const vector<double> &gf,
                           const vector<double> &dxc);

private:

    void setFlux(const vector<double> &gf,
                 const vector<double> &dxc);


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

public:

    dv_soot(domain      *line,
            const string s,
            const bool   Lt,
            const bool   Lo=true);

    virtual ~dv_soot(){}

};

////////////////////////////////////////////////////////////////////////////////