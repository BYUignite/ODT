/**
 * @file domaincase_hips_comb.cc
 * Header file for class domaincase_hips_comb
 */

#include "domaincase_hips_comb.h"
#include "domain.h"
#include "dv.h"
#include "dv_rho.h"
#include "dv_temp.h"
#include "dv_mixf_hips.h"
#include "dv_enth_hips.h"
#include "dv_ygas_hips.h"

////////////////////////////////////////////////////////////////////////////////
/** domaincase_hips_comb initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_hips_comb::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->atomicWeight(domn->gas->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->atomicWeight(domn->gas->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->atomicWeight(domn->gas->elementIndex("O"));
    gammas[3] = 0.0;
    domn->strm->init(domn, gammas);

    domn->v.push_back(new dv_mixf_hips(    domn, "mixf", true,  true ));    // use same order as for Sc below
    domn->v.push_back(new dv_rho(          domn, "rho",  false, true ));
    domn->v.push_back(new dv_temp(         domn, "temp", false, true ));
    domn->v.push_back(new dv_enth_hips(    domn, "enth", true,  true ));
    for(int k=0; k<domn->gas->nSpecies(); k++)
        domn->v.push_back(new dv_ygas_hips(domn, "y_"+domn->gas->speciesName(k), true, true ));

    int ii = 0;
    domn->mixf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->temp   = domn->v.at(ii++);
    domn->enth   = domn->v.at(ii++);
    domn->ysp    = domn->v.begin()+ii;          // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->nSpecies();

    //------------- set Sc of transported scalars: input file needs to be in the same order as above: mixf, enth, species

    if(domn->pram->LScHips){
        domn->v[0]->ScHips = domn->io->scalarSc[0].as<double>();
        domn->v[1]->ScHips = domn->io->scalarSc[1].as<double>();
        for(int k=0; k<domn->gas->nSpecies(); k++)
            domn->v[k+2]->ScHips = domn->io->scalarSc[k+2].as<double>();
    }

    //-------------------- initialize profiles

    double premixed_mixf      = domn->io->initParams["premixed_mixf"].as<double>();
    double frac_burnt         = domn->io->initParams["frac_burnt"].as<double>();

    for(int i=0; i<domn->ngrd; i++)
        domn->mixf->d.at(i) = premixed_mixf;
        //domn->mixf->d.at(i) = i<domn->ngrd/2 ? 0.055 : 0.055;

    int nsp = domn->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage

    for(int i=0; i<domn->ngrd; i++) {
        if(i < domn->ngrd*(1.0-frac_burnt))
            domn->strm->getMixingState(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        else
            domn->strm->getProdOfCompleteComb(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        for(int k=0; k<nsp; k++)
            domn->ysp[k]->d.at(i) = ysp.at(k);
    }




    //for(int i=0; i<domn->ngrd; i++)
    //    //domn->mixf->d.at(i) = double(i)/(domn->ngrd-1);
    //    if(i<256-23)
    //        domn->mixf->d.at(i) = 0.0;
    //    else if(i < 256-10)
    //        domn->mixf->d.at(i) = double(i-256+23)/12;
    //    else
    //        domn->mixf->d.at(i) = 1.0;
    //
    //int nsp = domn->gas->nSpecies();
    //vector<double> ysp(nsp);               // dummy storage
    //
    //for(int i=0; i<domn->ngrd; i++) {
    //    domn->strm->getProdOfCompleteComb(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
    //    for(int k=0; k<nsp; k++)
    //        domn->ysp[k]->d.at(i) = ysp.at(k);
    //}



    enforceMassFractions();

    domn->rho->setVar();

    //------------------- set minimial mesher

    vector<dv*> phi;
    domn->mesher->init(domn, phi);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_hips_comb::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = domn->ysp[k]->d.at(ipt);

    domn->gas->setState_PY(domn->pram->pres, &yi.at(0));
    domn->gas->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_hips_comb::setCaseSpecificVars() {

    enforceMassFractions();
    domn->rho->setVar();
    domn->temp->setVar();
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_hips_comb::setCaseSpecificVars_cvode(const int &ipt) {

    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}
