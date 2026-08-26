/**
 * @file domaincase_odt_wake.cc
 * @brief Source file for class domaincase_odt_wake
 */


#include "domaincase_odt_wake.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho_const.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"
#include "dv_ygas_noRxn.h"
#include "dv_mixf.h"
#include "dv_chi.h"
#include "dv_aDL.h"

#include "interp_linear.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_wake initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_wake::init(domain *p_domn) {

    domn = p_domn;

    domn->v.push_back(new dv_pos(         domn, "pos",     false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(        domn, "posf",    false, true ));
    domn->v.push_back(new dv_rho_const(   domn, "rho",     false, true ));
    domn->v.push_back(new dv_dvisc_const( domn, "dvisc",   false, true ));
    domn->v.push_back(new dv_uvw(         domn, "uvel",    true,  true ));
    domn->v.push_back(new dv_uvw(         domn, "vvel",    true,  true ));
    domn->v.push_back(new dv_uvw(         domn, "wvel",    true,  true ));

    int ii = 0;
    domn->pos    = domn->v.at(ii++);
    domn->posf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);
    domn->uvel   = domn->v.at(ii++);
    domn->vvel   = domn->v.at(ii++);
    domn->wvel   = domn->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    double r0   = domn->io->initParams["r0"].as<double>();   // half width half max
    double uinf = domn->io->initParams["uinf"].as<double>(); // free string velocity
    double du   = domn->io->initParams["du"].as<double>();   // velocity difference (uinf - wake center)

    for(int i=0; i<domn->ngrd; i++)
        domn->uvel->d[i] = uinf - du*exp( -0.5/r0/r0*domn->pos->d[i]*domn->pos->d[i] );

    //--------------------

    domn->rho->setVar();
    domn->dvisc->setVar();

    //-------------------

    if(domn->pram->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver needs to be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_wake::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();
}
