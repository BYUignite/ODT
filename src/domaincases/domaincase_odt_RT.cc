/**
 * @file domaincase_odt_RT.cc
 * @brief Source file for class domaincase_odt_RT
 */

#include "domaincase_odt_RT.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho_mf.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"
#include "dv_mixf.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_RT initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_RT::init(domain *p_domn) {

    domn = p_domn;

    domn->v.push_back(new dv_pos(         domn, "pos",     false, true  ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(        domn, "posf",    false, true  ));
    domn->v.push_back(new dv_uvw(         domn, "uvel",    true,  true  ));
    domn->v.push_back(new dv_uvw(         domn, "vvel",    true,  true  ));
    domn->v.push_back(new dv_uvw(         domn, "wvel",    true,  true  ));
    domn->v.push_back(new dv_mixf(        domn, "mixf",    true,  true  ));
    domn->v.push_back(new dv_rho_mf(      domn, "rho",     false, true  ));
    domn->v.push_back(new dv_dvisc_const( domn, "dvisc",   false, false ));

    int ii = 0;
    domn->pos    = domn->v.at(ii++);
    domn->posf   = domn->v.at(ii++);
    domn->uvel   = domn->v.at(ii++);
    domn->vvel   = domn->v.at(ii++);
    domn->wvel   = domn->v.at(ii++);
    domn->mixf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->mixf);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    double Zcent = domn->io->initParams["Zcent"].as<double>();
    double Ztran = domn->io->initParams["Ztran"].as<double>();

    Zcent += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd)); // Zcent=dist from center --> dist from left

    for(int i=0; i<domn->ngrd; i++)
        domn->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/Ztran*(domn->pos->d.at(i)-Zcent)));
        // domn->mixf->d.at(i) = (domn->pos->d.at(i) >= 0.5*(domn->pos->d[0]+domn->pos->d[domn->ngrd-1])) ? 1.0 : 0.0;

    //--------------------

    domn->rho->setVar();
    domn->dvisc->setVar();

    //-------------------

    if(domn->pram->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver should be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_RT::setCaseSpecificVars() {
    domn->rho->setVar();
    domn->dvisc->setVar();
}


