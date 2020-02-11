/**
 * @file domaincase_hips_simpleRxn.cc
 * Header file for class domaincase_hips_simpleRxn
 */

#include "domaincase_hips_simpleRxn.h"
#include "domain.h"
#include "dv.h"
#include "dv_mixf_hips.h"
#include "dv_ygas_cold_hips.h"

////////////////////////////////////////////////////////////////////////////////
/** domaincase_hips_simpleRxn initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_hips_simpleRxn::init(domain *p_domn) {

    domn = p_domn;

    domn->v.push_back(new dv_mixf_hips(     domn, "mixf", true, true ));
    domn->v.push_back(new dv_ygas_cold_hips(domn, "y_A",  true, true ));
    domn->v.push_back(new dv_ygas_cold_hips(domn, "y_B",  true, true ));
    domn->v.push_back(new dv_ygas_cold_hips(domn, "y_R",  true, true ));
    domn->v.push_back(new dv_ygas_cold_hips(domn, "y_P",  true, true ));

    int ii = 0;
    domn->mixf   = domn->v.at(ii++);
    domn->ysp    = domn->v.begin()+ii;          // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += 4;

    //------------- set Sc of transported scalars: input file needs to be in the same order as above: mixf, enth, species

    if(domn->pram->LScHips){
        domn->v[0]->ScHips = domn->io->scalarSc[0].as<double>();
        domn->v[1]->ScHips = domn->io->scalarSc[1].as<double>();
        domn->v[2]->ScHips = domn->io->scalarSc[2].as<double>();
        domn->v[3]->ScHips = domn->io->scalarSc[3].as<double>();
        domn->v[4]->ScHips = domn->io->scalarSc[4].as<double>();
    }

    //-------------------- initialize profiles

    for(int i=0; i<domn->ngrd; i++){
        domn->mixf->d.at(i) = i<domn->ngrd/2  ? 0.0 : 1.0;
        domn->ysp[0]->d[i]  = i<domn->ngrd/2  ? 1.0 : 0.0;
        domn->ysp[1]->d[i]  = i<domn->ngrd/2  ? 0.0 : 1.0;
        domn->ysp[2]->d[i]  = 0.0;
        domn->ysp[3]->d[i]  = 0.0;
    }

    enforceMassFractions();

    //------------------- set minimial mesher

    vector<dv*> phi;
    domn->mesher->init(domn, phi);

}

