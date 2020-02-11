/**
 * @file domaincase_odt_channel.cc
 * Header file for class domaincase_odt_channel
 */


#include "domaincase_odt_channel.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho_const.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"

////////////////////////////////////////////////////////////////////////////////
/** Initialization
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_channel::init(domain *p_domn){

    domn = p_domn;

    domn->v.push_back(new dv_pos(        domn, "pos",   false, true));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(       domn, "posf",  false, true));
    domn->v.push_back(new dv_rho_const(  domn, "rho",   false, false));
    domn->v.push_back(new dv_dvisc_const(domn, "dvisc", false, false));
    domn->v.push_back(new dv_uvw(        domn, "uvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "vvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "wvel",  true,  true));

    domn->pos   = domn->v.at(0);
    domn->posf  = domn->v.at(1);
    domn->rho   = domn->v.at(2);
    domn->dvisc = domn->v.at(3);
    domn->uvel  = domn->v.at(4);
    domn->vvel  = domn->v.at(5);
    domn->wvel  = domn->v.at(6);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    domn->mesher->init(domn, phi);

    //------------------- set inlet_cell_dv_props for inlet cell inserted for suction/blowing case

    inlet_cell_dv_props.resize(domn->v.size());

    inlet_cell_dv_props[0] = -1;                                     // pos:  set elsewhere
    inlet_cell_dv_props[1] = -domn->pram->domainLength/2.0;          // posf:
    inlet_cell_dv_props[2] = domn->pram->rho0;                       // rho:
    inlet_cell_dv_props[3] = domn->pram->kvisc0 * domn->pram->rho0;  // dvisc:
    inlet_cell_dv_props[4] = domn->pram->uBClo;                      // uvel:
    inlet_cell_dv_props[5] = domn->pram->vBClo;                      // vvel:
    inlet_cell_dv_props[6] = domn->pram->wBClo;                      // wvel:


    //------------------- default velocity values (0.0) are fine, along with rho, dvisc.

    //for(int i=0; i<domn->uvel->d.size(); i++)
    //  domn->uvel->d[i] = 10*domn->pos->d.at(i);
      //domn->uvel->d[i] = 10*0.016/4.0/0.002*(1.0-domn->pos->d.at(i)*domn->pos->d.at(i));

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_channel::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();

}
