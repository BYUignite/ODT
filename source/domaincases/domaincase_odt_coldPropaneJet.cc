/**
 * @file domaincase_odt_coldPropaneJet.cc
 * @brief Source file for class domaincase_odt_coldPropaneJet
 */

#include "domaincase_odt_coldPropaneJet.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_uvw.h"
#include "dv_ygas_noRxn.h"
#include "dv_mixf.h"
#include "dv_chi.h"
#include "dv_aDL.h"

#include "interp_linear.h"

#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_odt_coldPropaneJet initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_odt_coldPropaneJet::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 1.0;
    domn->strm->init(domn,gammas);

    domn->v.push_back(new dv_pos(   domn, "pos",     false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(  domn, "posf",    false, true ));
    domn->v.push_back(new dv_rho(   domn, "rho",     false, true ));
    domn->v.push_back(new dv_dvisc( domn, "dvisc",   false, true ));
    domn->v.push_back(new dv_uvw(   domn, "uvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "vvel",    true,  true ));
    domn->v.push_back(new dv_uvw(   domn, "wvel",    true,  true ));
    domn->v.push_back(new dv_mixf(  domn, "mixf",    false, true ));
    domn->v.push_back(new dv_chi(   domn, "chi",     false, true ));
    domn->v.push_back(new dv_aDL(   domn, "aDL",     false, false ));
    for(int k=0; k<domn->gas->thermo()->nSpecies(); k++)
        domn->v.push_back(new dv_ygas_noRxn(domn, "y_"+domn->gas->thermo()->speciesName(k), true, true ));

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
    domn->aDL    = domn->v.at(ii++);
    domn->ysp = domn->v.begin()+ii;       // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->thermo()->nSpecies();

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->uvel);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    double djeti = domn->io->initParams["djeti"].as<double>();
    double djeto = domn->io->initParams["djeto"].as<double>();

    vector<double> xprof, uprof;
    for(int i=0; i<domn->io->initParams["vprof"].size(); i++){
        xprof.push_back(domn->io->initParams["vprof"][i][0].as<double>());
        uprof.push_back(domn->io->initParams["vprof"][i][1].as<double>());
    }

    Linear_interp Linterp(xprof, uprof);
    for(int i=0; i<domn->ngrd; i++)
        domn->uvel->d.at(i) = Linterp.interp(domn->pos->d.at(i));

    for(int i=0; i<domn->ngrd; i++)
        domn->mixf->d.at(i) = (domn->pos->d.at(i) >= -djeti/2.0 && domn->posf->d.at(i) <= djeti/2.0) ? 1.0 : 0.0;

    //--------------------

    for(int i=0; i<domn->ngrd; i++)
        for(int k=0; k<domn->gas->thermo()->nSpecies(); k++)
            domn->ysp[k]->d.at(i) = domn->mixf->d.at(i) * domn->strm->y1[k] + (1.0-domn->mixf->d.at(i))*domn->strm->y0[k];

    enforceMassFractions();

    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->aDL->setVar();

    //-------------------

    if(domn->pram->Lsolver!="EXPLICIT") {
            cout << "\nError Lsolver needs to be EXPLICIT for this case (no rxn sources)" << endl;
            exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_odt_coldPropaneJet::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = domn->ysp[k]->d.at(ipt);
    domn->gas->thermo()->setState_TPY( domn->strm->T0, domn->pram->pres, &yi[0] );

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_odt_coldPropaneJet::setCaseSpecificVars() {

    enforceMassFractions();
    domn->rho->setVar();
    domn->dvisc->setVar();
}