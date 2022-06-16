/**
 * @file dv_rho.cc
 * @brief Source file for class dv_rho
 */


#include "dv_rho.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_rho::dv_rho(domain    *line,
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
/*! dv_rho merger2cells function
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

void dv_rho::merge2cells(const int    imrg,
                         const double m1,
                         const double m2,
                         const bool   LconstVolume) {

    try {
        domn->domc->setGasStateAtPt(imrg);
    } catch (const CanteraError& e) {
        throw odtCanteraError(STR_TRACE, "setGasStateAtPt", e);
    }
    d.at(imrg) = domn->gas->thermo()->density();
    d.erase(d.begin() + imrg+1);

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho setVar function
 *  @param ipt \input optional point to compute at
 */

void dv_rho::setVar(const int ipt){

    d.resize(domn->ngrd, domn->pram->rho0);
    if(ipt == -1)
        for(int i=0; i<domn->ngrd; i++) {
            try {
                domn->domc->setGasStateAtPt(i);
            } catch (const CanteraError& e) {
                cout << domn->gas->thermo()->density() << endl; //vbsdb
                throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
            }
            d.at(i) = domn->gas->thermo()->density();
        }
    else {
        try {
            domn->domc->setGasStateAtPt(ipt);
        } catch (const CanteraError& e) {
            cout << domn->gas->thermo()->density() << endl; //vbsdb
            throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
        }
        d.at(ipt) = domn->gas->thermo()->density();
    }
}