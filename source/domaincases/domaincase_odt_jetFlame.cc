/**
 * @file domaincase_odt_jetFlame.cc
 * @brief Source file for class domaincase_odt_jetFlame
 */

#include "domaincase_odt_jetFlame.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_uvw.h"
#include "dv_enth.h"
#include "dv_temp.h"
#include "dv_ygas.h"
#include "dv_mixf.h"
#include "dv_chi.h"
#include "dv_hr.h"
#include "dv_aDL.h"
#include "dv_soot.h"

#include "interp_linear.h"
#include "odtExceptions.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_jetFlame initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_jetFlame::init(domain *p_domn) {

    domn = p_domn;

    LisFlameD = domn->io->initParams["LisFlameD"] ? domn->io->initParams["LisFlameD"].as<bool>() : false;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("O"));
    gammas[3] = 0.0;
    if(LisFlameD) gammas[2] = 0.0;
    domn->strm->init(domn, gammas);

    domn->v.push_back(new dv_pos(   domn, "pos",     false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(  domn, "posf",    false, true ));
    domn->v.push_back(new dv_rho(   domn, "rho",     false, true ));
    domn->v.push_back(new dv_dvisc( domn, "dvisc",   false, true ));
    domn->v.push_back(new dv_uvw(   domn, "uvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "vvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "wvel",    true,  true ));
    domn->v.push_back(new dv_temp(  domn, "temp",    false, true ));
    domn->v.push_back(new dv_mixf(  domn, "mixf",    false, true ));
    domn->v.push_back(new dv_chi(   domn, "chi",     false, true ));
    domn->v.push_back(new dv_hr(    domn, "hr",      false, true ));
    if (domn->pram->LdoDL)
        domn->v.push_back(new dv_aDL(   domn, "aDL",     false, false ));
    for(int k=0; k<domn->gas->thermo()->nSpecies(); k++)
        domn->v.push_back(new dv_ygas(domn, "y_"+domn->gas->thermo()->speciesName(k), true, true ));
    domn->v.push_back(new dv_enth(  domn, "enth",    true,  true ));  // enth AFTER ygas for enth flux (see dv_enth)
    if (domn->pram->Lsoot) {
        for(int k=0; k<domn->pram->nsvar; k++) {
            domn->v.push_back(new dv_soot(domn, "M"+to_string(k), true, true ));
        }
    }

    int ii = 0;
    domn->pos    = domn->v.at(ii++);
    domn->posf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);
    domn->uvel   = domn->v.at(ii++);
    domn->vvel   = domn->v.at(ii++);
    domn->wvel   = domn->v.at(ii++);
    domn->temp   = domn->v.at(ii++);
    domn->mixf   = domn->v.at(ii++);
    domn->chi    = domn->v.at(ii++);
    domn->hr     = domn->v.at(ii++);
    if (domn->pram->LdoDL)
        domn->aDL    = domn->v.at(ii++);
    domn->ysp = domn->v.begin()+ii;       // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->thermo()->nSpecies();
    domn->enth   = domn->v.at(ii++);
    if (domn->pram->Lsoot == true) {
        domn->svar = domn->v.begin()+ii;    // access as domn->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += domn->pram->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->temp);
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

    //------------------- set initial soot profile
    if (domn->pram->Lsoot) {
        if (domn->pram->PSD_method == "QMOM" || domn->pram->PSD_method == "MOMIC" || domn->pram->PSD_method == "LOGN") {
            double M0 = 1.0E0;
            double sigL = 3.0;
            double mavg = 1.0E-21;
            for (int k=0; k<domn->pram->nsvar; k++) {
                for(int j=0; j<domn->ngrd; j++) {
                    domn->svar[k]->d[j] = M0 * pow(mavg, k) * exp(0.5 * pow(k,2) * pow(sigL,2));
                }
            }
        }
    }

    //--------------------

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    for(int i=0; i<domn->ngrd; i++) {
        if ( domn->pram->Lignition ) {
            /* non-reacting mixing for ignition */
            domn->strm->getMixingState(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        }
        else {
            /* products of combustion for equilibrium initial state */
            domn->strm->getProdOfCompleteComb(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        }
        for(int k=0; k<nsp; k++)
            domn->ysp[k]->d.at(i) = ysp.at(k);
    }

    enforceMassFractions();

    domn->rho->setVar();
    domn->dvisc->setVar();

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_jetFlame::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++) {
        yi.at(k) = domn->ysp[k]->d.at(ipt);
    }

    try {
        domn->gas->thermo()->setState_PY(domn->pram->pres, &yi.at(0));
    } catch (const CanteraError& c) {
        throw odtCanteraError(STR_TRACE, "setState_PY",c);
    }

    try {
        domn->gas->thermo()->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);
    } catch (const CanteraError& c) {
        throw odtCanteraError(STR_TRACE, "setState_HP",c);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_jetFlame::setCaseSpecificVars() {

    enforceMassFractions();
    domn->enth->setVar();
    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->temp->setVar();
    domn->chi->setVar();

    domn->enth->LagSrc = false;     // reset to false; in enth source the src is computed if false, then set to true on subsequent calls.
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_jetFlame::setCaseSpecificVars_cvode(const int &ipt) {

    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}