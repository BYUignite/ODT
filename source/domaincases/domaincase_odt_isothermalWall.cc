/**
 * @file domaincase_odt_isothermalWall.cc
 * Header file for class domaincase_odt_isothermalWall
 */

#include "domaincase_odt_isothermalWall.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_uvw.h"
#include "dv_enth.h"
#include "dv_temp.h"

////////////////////////////////////////////////////////////////////////////////
/** Initialization
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_isothermalWall::init(domain *p_domn){

    domn = p_domn;

    domn->v.push_back(new dv_pos(  domn, "pos",   false, true));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf( domn, "posf",  false, true));
    domn->v.push_back(new dv_rho(  domn, "rho",   false, true));
    domn->v.push_back(new dv_dvisc(domn, "dvisc", false, false));
    domn->v.push_back(new dv_temp( domn, "temp",  false, true ));
    domn->v.push_back(new dv_enth( domn, "enth",  true,  true ));
    domn->v.push_back(new dv_uvw(  domn, "uvel",  true,  true));
    domn->v.push_back(new dv_uvw(  domn, "vvel",  true,  true));
    domn->v.push_back(new dv_uvw(  domn, "wvel",  true,  true));

    int ii = 0;
    domn->pos   = domn->v.at(ii++);
    domn->posf  = domn->v.at(ii++);
    domn->rho   = domn->v.at(ii++);
    domn->dvisc = domn->v.at(ii++);
    domn->temp  = domn->v.at(ii++);
    domn->enth  = domn->v.at(ii++);
    domn->uvel  = domn->v.at(ii++);
    domn->vvel  = domn->v.at(ii++);
    domn->wvel  = domn->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->temp);
    domn->mesher->init(domn, phi);


    //--------------------  Initialize velocity field to a curve fit of a laminar boundary layer solution

    double nu = domn->pram->kvisc0;    // m2/s
    double vdiff    = domn->io->initParams["BL_vdiff"].as<double>();
    double delta    = domn->io->initParams["BL_width"].as<double>();

    double etafac = sqrt(vdiff/nu/delta);
    double eta;
    for(int i=0; i<domn->uvel->d.size(); i++) {
        eta = etafac * domn->pos->d.at(i);
        if(eta < 7.4)
            domn->uvel->d[i] = vdiff * (-1.55926E-4*pow(eta,5.0) +
                                         3.70924E-3*pow(eta,4.0) -
                                         2.84820E-2*pow(eta,3.0) +
                                         4.77151E-2*pow(eta,2.0) +
                                         3.05872E-1*eta +
                                         2.02326E-3 );
        else
            domn->uvel->d[i] = vdiff;
    }

    //------------------- set the gas composition and pressure, which won't change during the sim.
    // Also set temperature to initialize the enthalpy field.
    // Hard coding the composition to be air

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    yi[domn->gas->speciesIndex("O2")] = 0.233;
    yi[domn->gas->speciesIndex("N2")] = 0.767;
    domn->gas->setState_TPY(domn->pram->TBChi, domn->pram->pres, &yi.at(0));
    domn->enth->d = vector<double>(domn->ngrd, domn->gas->enthalpy_mass());
    domn->temp->d = vector<double>(domn->ngrd, domn->pram->TBChi);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_isothermalWall::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->temp->setVar();

    domn->enth->LagSrc = false;     // reset to false; in enth source the src is computed if false, then set to true on subsequent calls.
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_isothermalWall::setGasStateAtPt(const int &ipt) {

    domn->gas->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);

}
