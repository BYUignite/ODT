/**
 * @file dv_chi.cc
 * @brief Source file for class dv_chi
 */


#include "dv_chi.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_chi  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_chi::dv_chi(domain   *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/** Set scalar dissipation rate (chi) \cond
 *  @param ipt \input optional point to compute at
 *  Chi = 2D(df/dx)^2
 *  Compute the derivative as follows:
 *        *       *                 *
 *       i-1      i                i+1
 *            a            b
 *  Linearly interpolate the derivatives at a and b to point i
 *  f' = ( d2*fa' + d1*fb' ) / (d1+d2) where d1 is dist between i and a
 *  and d2 is distance between b and i \endcond
 */

void dv_chi::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd);
 
    domn->mixf->setVar();     // this may be redundant

    vector<double> gradZ(domn->ngrd);
    vector<double> Dthm(domn->ngrd);

    //-------------- Get thermal diffusivity

    for(int i=0; i<domn->ngrd; i++) {
        domn->domc->setGasStateAtPt(i);
        double tcond = domn->trn->thermalConductivity();
        double cp    = domn->gas->cp_mass();
        Dthm.at(i) = tcond/domn->rho->d.at(i)/cp;
    }

    //------------- Compute chi

    gradZ.at(0) = (domn->mixf->d.at(1)-domn->mixf->d.at(0))/(domn->pos->d.at(1)-domn->pos->d.at(0));
    double d1, d2;
    for(int i=1; i<domn->ngrd-1; i++) {
        d1 = 0.5*(domn->pos->d.at(i)  -domn->pos->d.at(i-1));
        d2 = 0.5*(domn->pos->d.at(i+1)-domn->pos->d.at(i));
        gradZ.at(i) = (d2*d2*(domn->mixf->d.at(i)  -domn->mixf->d.at(i-1)) +
                    d1*d1*(domn->mixf->d.at(i+1)-domn->mixf->d.at(i)))/(2*d1*d2*(d1+d2));
    }
    gradZ.at(domn->ngrd-1) = (domn->mixf->d.at(domn->ngrd-1)-domn->mixf->d.at(domn->ngrd-2)) /
                          (domn->pos->d.at(domn->ngrd-1)-domn->pos->d.at(domn->ngrd-2));

    for(int i=0; i<domn->ngrd; i++)
        d.at(i) = 2.0*Dthm.at(i)*gradZ.at(i)*gradZ.at(i);

}

