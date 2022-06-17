/**
 * @file domaincase_premixedFlameBurner_fixed_T.cc
 * Source file for class domaincase_premixedFlameBurner_fixed_T
 */

#include "domaincase_premixedFlameBurner_fixed_T.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_temp_profile.h"
#include "dv_hr.h"
#include "dv_ygas.h"
#include "soot/dv_soot.h"
#include "soot/dv_soot_MONO.h"
#include "soot/dv_soot_LOGN.h"
#include "soot/dv_soot_QMOM.h"
//#include "soot/dv_soot_CQMOM.h"
#include "soot/dv_soot_MOMIC.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_premixedFlameBurner_fixed_T initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_premixedFlameBurner_fixed_T::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->atomicWeight(domn->gas->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->atomicWeight(domn->gas->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->atomicWeight(domn->gas->elementIndex("O"));
    gammas[3] = 0.0;
    domn->strm = new streams(domn, gammas); domn->LstrmSet = true;

    domn->v.push_back(new dv_pos(          domn, "pos",     "--", false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(         domn, "posf",    "--", false, true ));
    domn->v.push_back(new dv_rho(          domn, "rho",     "--", false, true ));
    domn->v.push_back(new dv_dvisc(        domn, "dvisc",   "--", false, false ));
    domn->v.push_back(new dv_temp_profile( domn, "temp",    "--", false, true ));
    domn->v.push_back(new dv_hr(           domn, "hr",      "--", false, true ));
    for(int k=0; k<domn->gas->nSpecies(); k++)
        domn->v.push_back(new dv_ygas(     domn, "y_"+domn->gas->speciesName(k), "DN", true, true ));

    // Add soot moments to variable list
    if (domn->pram->Lsoot) {

        string PSD_method = domn->io->sootParams["PSD_method"].as<string>();
        stringstream ss;

        if (PSD_method == "MONO") {
            domn->pram->nsvar = 2;
            for(int k=0; k<2; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_MONO(domn, "M"+ss.str(), "DN", true, true ));
            }
        }
        else if (PSD_method == "LOGN") {
            domn->pram->nsvar = 3;
            for(int k=0; k<3; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_LOGN(domn, "M"+ss.str(), "DN", true, true ));
            }
        }
        else if (PSD_method == "QMOM") {
            for(int k=0; k<domn->pram->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_QMOM(domn, "M"+ss.str(), "DN", true, true ));
            }
        }
        else if (PSD_method == "MOMIC") {
            for(int k=0; k<domn->pram->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_MOMIC(domn, "M"+ss.str(), "DN", true, true ));
            }
        }   // end MOMIC
        //else if (PSD_method == "CQMOM") {
        //    domn->pram->nsvar = 2*domn->pram->nsvar_v/2*domn->pram->nsvar_s/2 + domn->pram->nsvar_v/2;  // stored and accessed by column
        //    for (int k=0; k<domn->pram->nsvar_v; k++) {
        //        domn->v.push_back(new dv_soot_flmlt_CQMOM(domn, "M"+to_string(k)+","+'0', "DN", true, true ));          // first s column: M00, M10, M20, etc.
        //    }
        //    for (int k=1; k<domn->pram->nsvar_s; k++) {
        //        for (int j=0; j<domn->pram->nsvar_v/2; j++) {
        //            domn->v.push_back(new dv_soot_flmlt_CQMOM(domn, "M"+to_string(j)+","+to_string(k), "DN", true, true ));       // other s columns: M01, M11, M02, M12, etc.
        //        }
        //    }
        //}   // end CQMOM

    }

    int ii = 0;
    domn->pos    = domn->v.at(ii++);
    domn->posf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);
    domn->temp   = domn->v.at(ii++);
    domn->hr     = domn->v.at(ii++);
    domn->ysp = domn->v.begin()+ii;       // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->nSpecies();
    if (domn->pram->Lsoot == true) {
        domn->svar     = domn->v.begin()+ii;    // access as domn->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += domn->pram->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->temp);
    domn->mesher = new meshManager(domn, phi); domn->LmeshSet = true;

    //------------------- set initial profiles

    // calculate initial mixf at user specified value of phi

    double phi0      = domn->io->initParams["phi0"].as<double>();
    double fracBurnt = domn->io->initParams["fracBurnt"].as<double>();

    double FA_st = domn->strm->mixfStoic/(1-domn->strm->mixfStoic);

    mixf_reactants = phi0*FA_st/(1+phi0*FA_st);

    domn->temp->setVar();

    // initialize composition profiles, set burnt for fracBurnt

    int nsp = domn->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    double hdummy;
    double Tdummy;

    for(int i=0; i<domn->ngrd; i++) {
        if(domn->pos->d.at(i) < (1-fracBurnt)*domn->pram->domainLength)
            domn->strm->getMixingState(mixf_reactants, ysp, hdummy, Tdummy);
        else {
            if(domn->pram->chemMechFile=="onestep_c2h4.xml")
                domn->strm->getProdOfCompleteComb(mixf_reactants, ysp, hdummy, Tdummy);
            else
                domn->strm->getEquilibrium_TP(mixf_reactants, domn->temp->d.at(i), ysp, hdummy);
        }
        for(int k=0; k<nsp; k++)
            domn->ysp[k]->d.at(i) = ysp.at(k);
    }

    enforceMassFractions();

    domn->rho->setVar();
    domn->dvisc->setVar();

    //------------------- set needed boundary conditions

    //----------- burner side (lo) Dirichlet conditions

    domn->strm->getMixingState(mixf_reactants, ysp, hdummy, Tdummy); // this also sets the gas object state

    domn->rho->bcLo = domn->gas->density();
    domn->dvisc->bcLo = domn->tran->viscosity();
    for(int k=0; k<domn->gas->nSpecies(); ++k)
        domn->ysp[k]->bcLo = ysp[k];

    if(domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            domn->svar[k]->bcLo = 0.0;

    domn->temp->bcLo = domn->temp->d[0];

    //----------- outflow side (hi) Neumann conditions

    for(int k=0; k<domn->gas->nSpecies(); ++k)
        domn->ysp[k]->bcHi = 0.0;

    if(domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            domn->svar[k]->bcHi = 0.0;

    domn->temp->bcHi = domn->temp->d.back();

    //------------------- set inlet_cell_dv_props for inlet cell inserted on left as inlet condition

    inlet_cell_dv_props.resize(domn->v.size());

    ii = 0;
    inlet_cell_dv_props[ii++] = -1;                     // pos:  set explicitly in micromixer_premix::updateGrid
    inlet_cell_dv_props[ii++] = 0.0;                    // posf
    inlet_cell_dv_props[ii++] = domn->rho->bcLo;        // rho
    inlet_cell_dv_props[ii++] = domn->dvisc->bcLo;      // dvisc
    inlet_cell_dv_props[ii++] = domn->temp->bcLo;       // temp
    inlet_cell_dv_props[ii++] = 0.0;                    // hr
    for(int k=0; k<domn->gas->nSpecies(); k++, ii++)   
        inlet_cell_dv_props[ii] = domn->ysp[k]->bcLo;   // gas species
    if (domn->pram->Lsoot)
        for(int k=0; k<domn->pram->nsvar; k++)
            inlet_cell_dv_props[ii++] = domn->svar[k]->bcLo;  // soot

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_premixedFlameBurner_fixed_T::setCaseSpecificVars() {

    domn->temp->setVar();
    enforceSootMom();
    enforceMassFractions();
    domn->rho->setVar();
    domn->dvisc->setVar();

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_premixedFlameBurner_fixed_T::setCaseSpecificVars_cvode(const int &ipt) {
    domn->rho->setVar(ipt);
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 *  If ipt=-1 set based on the left boundary.
 *  If ipt=ngrd+1=ngrdf set based on the right boundary.
 */

void domaincase_premixedFlameBurner_fixed_T::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);

    double T;
    if(ipt==-1)
        T = domn->temp->bcLo;
    else if(ipt==domn->ngrd)
        T = domn->temp->bcHi;
    else{
        domn->temp->setVar(ipt);
        T = domn->temp->d.at(ipt);
    }
    for(int k=0; k<nsp; k++)
        yi[k] = (ipt==-1 ? domn->ysp[k]->bcLo : (ipt==domn->ngrdf ? domn->ysp[k]->bcHi : domn->ysp[k]->d.at(ipt)) );

    domn->gas->setState_TPY(T, domn->pram->pres, &yi.at(0));
}
