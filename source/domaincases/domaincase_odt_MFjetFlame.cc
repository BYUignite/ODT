/**
 * @file domaincase_odt_MFjetFlame.cc
 * @brief Source file for class domaincase_odt_MFjetFlame
 */

#include "domaincase_odt_MFjetFlame.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
/* customize the mixture fraction profile here */
#include "dv_rho_mf.h"
#include "dv_dvisc_const.h"
#include "dv_uvw.h"
#include "dv_mixf.h"
#include "dv_chi_dmf.h"

#include "interp_linear.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_MFjetFlame initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_MFjetFlame::init(domain *p_domn) {

    domn = p_domn;

    domn->v.push_back(new dv_pos(   domn, "pos",     false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(  domn, "posf",    false, true ));
    domn->v.push_back(new dv_rho_mf(   domn, "rho",     false, true ));
    domn->v.push_back(new dv_dvisc_const( domn, "dvisc",   false, false ));
    domn->v.push_back(new dv_uvw(   domn, "uvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "vvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "wvel",    true,  true ));
    domn->v.push_back(new dv_mixf(  domn, "mixf",    true, true ));
    domn->v.push_back(new dv_chi_dmf(   domn, "chi",     false, true ));

    int ii = 0;
    domn->pos    = domn->v.at(ii++);
    domn->posf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);
    domn->uvel   = domn->v.at(ii++);
    domn->vvel   = domn->v.at(ii++);
    domn->wvel   = domn->v.at(ii++);
    domn->mixf   = domn->v.at(ii++);
    domn->chi    = domn->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->mixf);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    double d_f    = domn->io->initParams["d_f"].as<double>();
    double d_p    = domn->io->initParams["d_p"].as<double>();
    double Z_p    = domn->io->initParams["Z_p"].as<double>();
    double dTrans = domn->io->initParams["dTrans"].as<double>();
    double U_f    = domn->io->initParams["U_f"].as<double>();
    double U_p    = domn->io->initParams["U_p"].as<double>();
    double U_a    = domn->io->initParams["U_a"].as<double>();

    //--------------------

    double r_f = d_f/2.0;
    double r_p = d_p/2.0;
    double p1;                 // dummy for setting profiles
    double p2;                 // dummy for setting profiles

    for(int i=0; i<domn->ngrd; i++){
        p1 = 0.5*(1.0+tanh(2.0/dTrans*(domn->pos->d.at(i)+r_f))) *
             0.5*(1.0+tanh(2.0/dTrans*(r_f-domn->pos->d.at(i))));
        p2 = 0.5*(1.0+tanh(2.0/dTrans*(domn->pos->d.at(i)+r_p))) *
             0.5*(1.0+tanh(2.0/dTrans*(r_p-domn->pos->d.at(i))));
        domn->mixf->d.at(i) = p1*(1.0-Z_p) + p2*Z_p;
        if(U_f > 0)
            domn->uvel->d.at(i) = p1*(U_f-U_p-U_a) + p2*(U_p-U_a) + U_a;
    }

    if(U_f < 0) {
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

    if(domn->pram->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver should be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_MFjetFlame::setGasStateAtPt(const int &ipt) {

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_MFjetFlame::setCaseSpecificVars() {

    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->chi->setVar();

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_MFjetFlame::setCaseSpecificVars_cvode(const int &ipt) {
    domn->rho->setVar(ipt);
}
