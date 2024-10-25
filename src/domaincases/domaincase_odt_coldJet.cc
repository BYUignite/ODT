/**
 * @file domaincase_odt_coldJet.cc
 * @brief Source file for class domaincase_odt_coldJet
 */


#include "domaincase_odt_coldJet.h"
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
/** domaincase_odt_coldJet initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_coldJet::init(domain *p_domn) {

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

    double djeti      = domn->io->initParams["djeti"].as<double>();
    double vel_min    = domn->io->initParams["vel_min"].as<double>();
    double vel_max    = domn->io->initParams["vel_max"].as<double>();
    double delta_vel  = domn->io->initParams["delta_vel"].as<double>();
    double vel_diff   = vel_max - vel_min;

    if(djeti >= 0.0){
        double vyc1 = -0.5*djeti;
        double vyc2 =  0.5*djeti;
        for(int i=0; i<domn->ngrd; i++){
            domn->uvel->d.at(i) = vel_diff * 0.5 * (1.0 + tanh(2.0 / delta_vel * (domn->pos->d.at(i) - vyc1))) *
                0.5 * (1.0 + tanh(2.0 / delta_vel * (vyc2 - domn->pos->d.at(i)))) + vel_min;

            //domn->uvel->d.at(i) = 10.0 + 1.0*domn->pos->d.at(i); //doldb
        }
    }

    //--------- hack in a vprof profile

    else {
        vector<double> xprof, uprof;
        for(int i=0; i<domn->io->initParams["vprof"].size(); i++){
            xprof.push_back(domn->io->initParams["vprof"][i][0].as<double>());
            uprof.push_back(domn->io->initParams["vprof"][i][1].as<double>());
        }
        Linear_interp Linterp(xprof, uprof);
        for(int i=0; i<domn->ngrd; i++)
            domn->uvel->d.at(i) = Linterp.interp(domn->pos->d.at(i));
    }

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
void domaincase_odt_coldJet::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();
}

