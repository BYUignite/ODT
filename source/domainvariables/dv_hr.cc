/**
 * @file dv_hr.cc
 * @brief Source file for class dv_hr
 */

#include "dv_hr.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

//void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/*! dv_hr  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_hr::dv_hr(domain *line,
        const       string s,
        const bool  Lt,
        const bool  Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! Set heat release rate (J/m3*s) from the gas state
 *  @param ipt \input optional point to compute at
 */

void dv_hr::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd);

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> rr(nsp);
    vector<double> yi(nsp);
    vector<double> hsp(nsp);

    for(int i=0; i<domn->ngrd; i++){
        try {
            domn->domc->setGasStateAtPt(i);
        } catch (const CanteraError& e) {
            throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
        }
        // catch extremes of temperature due to diff-diff
        double temperatureHere = domn->gas->thermo()->temperature();
        temperatureHere = ( temperatureHere < 250.0 ) ? 250.0 : temperatureHere;

        domn->gas->thermo()->getMassFractions(&yi[0]);
        //        getProblemSpecificRR(domn->gas->density(), domn->gas->temperature(), domn->pram->pres, &yi.at(0), &rr.at(0));
        domn->chem->getProblemSpecificRR(domn->gas->thermo()->density(), temperatureHere, domn->pram->pres, &yi.at(0), &rr.at(0));

        domn->gas->thermo()->getEnthalpy_RT(&hsp.at(0));               // non-dimensional enthalpy
        d.at(i) = 0.0;
        for(int k=0; k<nsp; k++)
            d.at(i) -= rr.at(k)*(hsp.at(k)*temperatureHere*GasConstant);
    }
}