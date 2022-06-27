/**
 * @file domaincase_premixedFlameBurner.cc
 * Source file for class domaincase_premixedFlameBurner
 */

#include "domaincase_premixedFlameBurner.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_temp.h"
#include "dv_hr.h"
#include "dv_ygas.h"
#include "dv_enth.h"
#include "dv_soot.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_premixedFlameBurner initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_premixedFlameBurner::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->thermo()->atomicWeight (domn->gas->thermo()->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->thermo()->atomicWeight (domn->gas->thermo()->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->thermo()->atomicWeight(domn->gas->thermo()->elementIndex("O"));
    gammas[3] = 0.0;
    domn->strm = new streams();
    domn->strm->init(domn, gammas);
    domn->LstrmSet = true;

    domn->v.push_back(new dv_pos(     domn, "pos",   false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(    domn, "posf",  false, true ));
    domn->v.push_back(new dv_rho(     domn, "rho",   false, true ));
    domn->v.push_back(new dv_dvisc(   domn, "dvisc", false, false ));
    domn->v.push_back(new dv_temp(    domn, "temp",  false, true ));
    domn->v.push_back(new dv_hr(      domn, "hr",    false, true ));
    for(int k=0; k<domn->gas->thermo()->nSpecies(); k++)
        domn->v.push_back(new dv_ygas(domn, "y_"+domn->gas->thermo()->speciesName(k), true, true ));
    domn->v.push_back(new dv_enth(    domn, "enth", true,  true ));  // enth AFTER ygas for enth flux (see dv_enth)
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
    domn->temp   = domn->v.at(ii++);
    domn->hr     = domn->v.at(ii++);
    domn->ysp = domn->v.begin()+ii;       // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->thermo()->nSpecies();
    domn->enth   = domn->v.at(ii++);
    if (domn->pram->Lsoot == true) {
        domn->svar = domn->v.begin()+ii;    // access as domn->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += domn->pram->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->temp);
    domn->mesher = new meshManager();
    domn->mesher->init(domn, phi);
    domn->LmeshSet = true;

    //------------------- set initial profiles

    // calculate initial mixf at user specified value of phi

    double phi0      = domn->io->initParams["phi0"].as<double>();
    double fracBurnt = domn->io->initParams["fracBurnt"].as<double>();

    double FA_st = domn->strm->mixfStoic/(1-domn->strm->mixfStoic);

    mixf_reactants = phi0*FA_st/(1+phi0*FA_st);

    // initialize composition profiles, set burnt for fracBurnt

    int nsp = domn->gas->thermo()->nSpecies();
    vector<double> ysp(nsp);               // dummy storage

    for(int i=0; i<domn->ngrd; i++) {
        if(domn->pos->d.at(i) < (1-fracBurnt)*domn->pram->domainLength)
            domn->strm->getMixingState(mixf_reactants, ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        else {
            if(domn->pram->chemMech=="onestep_c2h4")
                domn->strm->getProdOfCompleteComb(mixf_reactants, ysp, domn->enth->d.at(i), domn->temp->d.at(i));
            else
                domn->strm->getEquilibrium_HP(mixf_reactants, ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        }
        for(int k=0; k<nsp; k++)
            domn->ysp[k]->d.at(i) = ysp.at(k);
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

    enforceMassFractions();

    domn->rho->setVar();
    domn->dvisc->setVar();

    //------------------- set needed boundary conditions

    //----------- burner side (lo) Dirichlet conditions

    domn->strm->getMixingState(mixf_reactants, ysp, domn->enth->bcLo, domn->temp->bcLo);
    // this also sets the gas object state

    domn->rho->bcLo = domn->gas->thermo()->density();
    domn->dvisc->bcLo = domn->gas->transport()->viscosity();
    for(int k=0; k<domn->gas->thermo()->nSpecies(); ++k)
        domn->ysp[k]->bcLo = ysp[k];

    if(domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            domn->svar[k]->bcLo = 0.0;

    //----------- outflow side (hi) Neumann conditions

    for(int k=0; k<domn->gas->thermo()->nSpecies(); ++k)
        domn->ysp[k]->bcHi = 0.0;

    domn->enth->bcHi = 0.0;

    if(domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            domn->svar[k]->bcHi = 0.0;

    //------------------- set inlet_cell_dv_props for inlet cell inserted on left as inlet condition

    inlet_cell_dv_props.resize(domn->v.size());

    ii = 0;
    inlet_cell_dv_props[ii++] = -1;                     // pos:  set explicitly in micromixer_premix::updateGrid
    inlet_cell_dv_props[ii++] = 0.0;                    // posf
    inlet_cell_dv_props[ii++] = domn->rho->bcLo;        // rho
    inlet_cell_dv_props[ii++] = domn->dvisc->bcLo;      // dvisc
    inlet_cell_dv_props[ii++] = domn->temp->bcLo;       // temp
    inlet_cell_dv_props[ii++] = 0.0;                    // hr
    for(int k=0; k<domn->gas->thermo()->nSpecies(); k++, ii++)
        inlet_cell_dv_props[ii] = domn->ysp[k]->bcLo;   // gas species
    inlet_cell_dv_props[ii++] = domn->enth->bcLo;       // enthalpy
    if (domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            inlet_cell_dv_props[ii++] = domn->svar[k]->bcLo;  // soot

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_premixedFlameBurner::setCaseSpecificVars() {

    enforceSootMom();
    enforceMassFractions();
    domn->enth->setVar();
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
void domaincase_premixedFlameBurner::setCaseSpecificVars_cvode(const int &ipt) {
    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}