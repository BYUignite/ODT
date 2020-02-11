/**
 * @file domaincase_hips.cc
 * Header file for class domaincase_hips
 */

#include "domaincase_hips.h"
#include "domain.h"
#include "dv.h"
#include "dv_mixf_hips.h"
#include <sstream>
#include <iomanip>
#include <string>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_hips initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_hips::init(domain *p_domn) {

    domn = p_domn;

    if(domn->pram->LScHips){
        for(int i=0; i<domn->io->scalarSc.size(); i++){
            stringstream ss;
            ss << setfill('0') << setw(3) << i;
            string varName = "mixf_" + ss.str();
            domn->v.push_back(new dv_mixf_hips( domn, varName, true, true ));
            domn->v[i]->ScHips = domn->io->scalarSc[i].as<double>();
        }
    }
    else
        domn->v.push_back(new dv_mixf_hips( domn, "mixf", true, true ));

    int ii = 0;
    domn->mixf   = domn->v.at(ii++);

    //-------------------- initialize profiles

    for(int i=0; i<domn->ngrd; i++){
        domn->mixf->d.at(i) = i<domn->ngrd/2 ? 0.0 : 1.0;
        for(int k=0; k<domn->v.size(); k++) {
            if(!domn->v[k]->L_transported)
                continue;
            domn->v[k]->d[i] = i<domn->ngrd/2 ? 0.0 : 1.0;
        }
    }

    //------------------- set minimial mesher

    vector<dv*> phi;
    domn->mesher->init(domn, phi);

    //------------------- for variable Sc, set the Batchelor level and i_plus level for each scalar

    double iEta = domn->io->params["nLevels"].as<int>() - 3;       // nLevels may have been increased in solver_hips initializer; want the original
    for(int k=0; k<domn->v.size(); k++) {
        if(domn->v[k]->L_transported){
            if(domn->v[k]->ScHips < 1.0)
                domn->v[k]->i_batchelor = iEta + 1.5*log(domn->v[k]->ScHips)/log(4);
            else if(domn->v[k]->ScHips > 1.0)
                domn->v[k]->i_batchelor = iEta +     log(domn->v[k]->ScHips)/log(4);
            else
                domn->v[k]->i_batchelor = iEta;
            domn->v[k]->i_plus = ceil(domn->v[k]->i_batchelor);
            cout << endl << "Scalar index, Sc, i_batchelor, i_plus: " << k << " " << domn->v[k]->ScHips << " " << domn->v[k]->i_batchelor << " " << domn->v[k]->i_plus << endl;
        }
    }

}

