/**
 * @file domaincase_odt_jetMixlRxn.cc
 * @brief Source file for class domaincase_odt_jetMixlRxn
 */

#include "domaincase_odt_jetMixlRxn.h"
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

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_jetMixlRxn initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_jetMixlRxn::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("O"));
    gammas[3] = 0.0;
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
    for(int k=0; k<domn->gas->thermo()->nSpecies(); k++)
        domn->v.push_back(new dv_ygas(domn, "y_"+domn->gas->thermo()->speciesName(k), true, true ));
    domn->v.push_back(new dv_enth(  domn, "enth",    true,  true ));  // enth AFTER ygas for enth flux (see dv_enth)

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
    domn->ysp = domn->v.begin()+ii;       // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->thermo()->nSpecies();
    domn->enth   = domn->v.at(ii++);

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    phi.push_back(domn->temp);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    string jetOrMixl  = domn->io->initParams["jetOrMixl"].as<string>();
    double delta_mixf = domn->io->initParams["delta_mixf"].as<double>();
    double fyc1       = domn->io->initParams["fyc1"].as<double>();
    double delta_vel  = domn->io->initParams["delta_vel"].as<double>();
    double vyc1       = domn->io->initParams["vyc1"].as<double>();
    double vel_min    = domn->io->initParams["vel_min"].as<double>();
    double vel_max    = domn->io->initParams["vel_max"].as<double>();
    double vel_diff   = vel_max - vel_min;

    //--------------------

    if(jetOrMixl == "MIXL") {    // mixing layer

        fyc1 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd)); // fyc1=dist from center --> dist from left
        for(int i=0; i<domn->ngrd; i++)
            domn->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/delta_mixf*(domn->pos->d.at(i)-fyc1)));

        vyc1 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd));  // vyc1=dist from center --> dist from left
        for(int i=0; i<domn->ngrd; i++)
            domn->uvel->d.at(i) = vel_diff*0.5*(1.0+tanh(2.0/delta_vel*(domn->pos->d.at(i)-vyc1))) + vel_min;
    }
    //--------------------

    else if(jetOrMixl == "JET") { // jet

        double fyc2       = domn->io->initParams["fyc2"].as<double>();
        double vyc2       = domn->io->initParams["vyc2"].as<double>();
        fyc1 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd)); // fyc1=dist from center --> dist from left
        fyc2 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd)); // fyc2=dist from center --> dist from left
        for(int i=0; i<domn->ngrd; i++){
            domn->mixf->d.at(i) = 0.5*(1.0+tanh(2.0/delta_mixf*(domn->pos->d.at(i)-fyc1))) *
                0.5*(1.0+tanh(2.0/delta_mixf*(fyc2-domn->pos->d.at(i))));
        }

        vyc1 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd));  // vyc1=dist from center --> dist from left
        vyc2 += 0.5*(domn->posf->d.at(0)+domn->posf->d.at(domn->ngrd));  // vyc2=dist from center --> dist from left
        for(int i=0; i<domn->ngrd; i++){
            domn->uvel->d.at(i) = vel_diff * 0.5 * (1.0 + tanh(2.0 / delta_vel * (domn->pos->d.at(i) - vyc1))) *
                0.5 * (1.0 + tanh(2.0 / delta_vel * (vyc2 - domn->pos->d.at(i)))) + vel_min;
        }
    }
    else {
        *domn->io->ostrm << endl << "Error setting the domain: option " << jetOrMixl
            << " not recognized " << endl;
        exit(0);
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

void domaincase_odt_jetMixlRxn::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = domn->ysp[k]->d.at(ipt);
    domn->gas->thermo()->setState_PY(domn->pram->pres, &yi.at(0));
    domn->gas->thermo()->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_jetMixlRxn::setCaseSpecificVars() {

    enforceMassFractions();
    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->temp->setVar();

    domn->enth->LagSrc = false;     // reset to false; in enth source the src is computed if false, then set to true on subsequent calls.
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_jetMixlRxn::setCaseSpecificVars_cvode(const int &ipt) {
    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}