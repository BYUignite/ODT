/**
 * @file domaincase_odt_channelScalar.cc
 * Header file for class domaincase_odt_channelScalar
 */


#include "domaincase_odt_channelScalar.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho_const.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"
#include "dv_sca.h"

////////////////////////////////////////////////////////////////////////////////
/** Initialization
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_channelScalar::init(domain *p_domn){

    domn = p_domn;

    domn->v.push_back(new dv_pos(        domn, "pos",   false, true));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(       domn, "posf",  false, true));
    domn->v.push_back(new dv_rho_const(  domn, "rho",   false, false));
    domn->v.push_back(new dv_dvisc_const(domn, "dvisc", false, false));
    domn->v.push_back(new dv_uvw(        domn, "uvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "vvel",  true,  true));
    domn->v.push_back(new dv_uvw(        domn, "wvel",  true,  true));
    domn->v.push_back(new dv_sca(        domn, "sca",   true,  true));

    int k=0;
    domn->pos   = domn->v.at(k++);
    domn->posf  = domn->v.at(k++);
    domn->rho   = domn->v.at(k++);
    domn->dvisc = domn->v.at(k++);
    domn->uvel  = domn->v.at(k++);
    domn->vvel  = domn->v.at(k++);
    domn->wvel  = domn->v.at(k++);
    domn->sca   = domn->v.at(k++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->sca);
    domn->mesher->init(domn, phi);

    //------------------- default vel. & sca. values (0.0) are fine, along with rho, dvisc.

    for(int i=0; i<domn->ngrd; i++) {
        //---------- nonzero velocity profile (mind the BCs)
        //domn->uvel->d[i] = 10*domn->pos->d.at(i); //doldbg
        //domn->uvel->d[i] = 10*0.016/4.0/0.002*(1.0-domn->pos->d.at(i)*domn->pos->d.at(i)); //doldbg
        //---------- non-zero scalar profile (mind the BCs)
        double xm = (domn->pos->d.at(i) - domn->posf->d.at(0))/domn->Ldomain();
        domn->sca->d.at(i) = domn->pram->sBClo * (1.0-xm) + domn->pram->sBChi * xm;
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_channelScalar::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();

}
