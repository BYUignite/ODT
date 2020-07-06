/**
 * @file dv_rho_const.cc
 * @brief Source file for class dv_rho_const
 */


#include "dv_rho_const.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho_const  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_rho_const::dv_rho_const(domain    *line,
                           const      string s,
                           const bool Lt,
                           const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, domn->pram->rho0);

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho_const merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg    \input merge cells imrg and imrg+1
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void dv_rho_const::merge2cells(const int    imrg,
                               const double m1,
                               const double m2,
                               const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho_const setVar function
 *  @param ipt \input optional point to compute at
 */

void dv_rho_const::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd, domn->pram->rho0);
}

