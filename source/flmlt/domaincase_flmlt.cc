/**
 * @file domaincase_flmlt.cc
 * Header file for class domaincase_flmlt
 */

#include "domaincase_flmlt.h"
#include "domain.h"
#include "dv.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_enth_flmlt.h"
#include "dv_temp.h"
#include "dv_ygas_flmlt.h"
#include "dv_chi_flmlt.h"
#include "dv_soot_flmlt_MONO.h"
#include "dv_soot_flmlt_MOMIC.h"
#include "dv_soot_flmlt_LOGN.h"
#include "dv_soot_flmlt_QMOM.h"
//#include "dv_soot_flmlt_CQMOM.h"

#include "interp_linear.h"

#include <cmath>
#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_flmlt initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_flmlt::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->atomicWeight(domn->gas->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->atomicWeight(domn->gas->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->atomicWeight(domn->gas->elementIndex("O"));
    gammas[3] = 0.0;
    domn->strm->init(domn, gammas);

    domn->v.push_back(new dv_pos(           domn, "pos",     false, true ));   // last are: L_transported, L_output
    domn->v.push_back(new dv_posf(          domn, "posf",    false, true ));
    domn->v.push_back(new dv_rho(           domn, "rho",     false, true ));
    domn->v.push_back(new dv_dvisc(         domn, "dvisc",   false, true ));
    domn->v.push_back(new dv_temp(          domn, "temp",    false, true ));
    domn->v.push_back(new dv_chi_flmlt(     domn, "chi",     false, true ));
    for(int k=0; k<domn->gas->nSpecies(); k++)
        domn->v.push_back(new dv_ygas_flmlt(domn, "y_"+domn->gas->speciesName(k), true, true ));
    domn->v.push_back(new dv_enth_flmlt(    domn, "enth",    false,  true ));

    // Add soot moments to variable list
    if (domn->pram->Lsoot) {

        string PSD_method = domn->io->sootParams["PSD_method"].as<string>();
        stringstream ss;

        if (PSD_method == "MONO") {
            domn->pram->nsvar = 2;
            for(int k=0; k<2; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_flmlt_MONO(domn, "M"+ss.str(), true, true ));
            }
        }
        else if (PSD_method == "LOGN") {
            domn->pram->nsvar = 3;
            for(int k=0; k<3; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_flmlt_LOGN(domn, "M"+ss.str(), true, true ));
            }
        }
        else if (PSD_method == "QMOM") {
            for(int k=0; k<domn->pram->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_flmlt_QMOM(domn, "M"+ss.str(), true, true ));
            }
        }
        else if (PSD_method == "MOMIC") {
            for(int k=0; k<domn->pram->nsvar; k++) {
                ss.str(""); ss.clear(); ss << k;
                domn->v.push_back(new dv_soot_flmlt_MOMIC(domn, "M"+ss.str(), true, true ));
            }
        }   // end MOMIC
        //else if (PSD_method == "CQMOM") {
        //    domn->pram->nsvar = 2*domn->pram->nsvar_v/2*domn->pram->nsvar_s/2 + domn->pram->nsvar_v/2;  // stored and accessed by column
        //    for (int k=0; k<domn->pram->nsvar_v; k++) {
        //        domn->v.push_back(new dv_soot_flmlt_CQMOM(domn, "M"+to_string(k)+","+'0', true, true ));          // first s column: M00, M10, M20, etc.
        //    }
        //    for (int k=1; k<domn->pram->nsvar_s; k++) {
        //        for (int j=0; j<domn->pram->nsvar_v/2; j++) {
        //            domn->v.push_back(new dv_soot_flmlt_CQMOM(domn, "M"+to_string(j)+","+to_string(k), true, true ));       // other s columns: M01, M11, M02, M12, etc.
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
    domn->chi    = domn->v.at(ii++);
    domn->ysp    = domn->v.begin()+ii;          // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->nSpecies();
    domn->enth   = domn->v.at(ii++);
    if (domn->pram->Lsoot == true) {
        domn->svar     = domn->v.begin()+ii;    // access as domn->svar[k]->d[i], etc. where k is the species starting from 0.
        ii += domn->pram->nsvar;
    }

    //------------------- set variables used for mesh adaption

    vector<dv*> phi;
    phi.push_back(domn->temp);
    domn->mesher->init(domn, phi);

    //------------------- set profiles

    if(domn->pram->domainLength != 1.0) {
        cout << endl << "ERROR: for flmlt, domainLength is mixf and must be 1" << endl;
        exit(0);
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

    double dx = domn->pram->domainLength / domn->ngrd;
    domn->posf->d.at(0) = 0.0;
    for(int i=1; i<domn->ngrdf; i++)
        domn->posf->d.at(i) = domn->posf->d.at(i-1) + dx;
    domn->posf->d.at(domn->ngrd) = domn->pram->domainLength;
    domn->pos->setVar();

    int nsp = domn->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage
    for(int i=0; i<domn->ngrd; i++) {
        if ( domn->pram->Lignition ) {
            /* non-reacting mixing for ignition */
            domn->strm->getMixingState(domn->pos->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        }
        else {
            /* products of combustion for equilibrium initial state */
            domn->strm->getProdOfCompleteComb(domn->pos->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
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

void domaincase_flmlt::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++) {
        yi.at(k) = domn->ysp[k]->d.at(ipt);
    }

    domn->gas->setState_PY(domn->pram->pres, &yi.at(0));
    domn->gas->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_flmlt::setCaseSpecificVars() {

    enforceSootMom();
    enforceMassFractions();
    domn->enth->setVar();
    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->temp->setVar();
    domn->chi->setVar();
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_flmlt::setCaseSpecificVars_cvode(const int &ipt) {

    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}
