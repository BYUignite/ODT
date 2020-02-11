/**
 * @file domaincase.cc
 * Header file for class domaincase
 */

#include "domaincase.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/** Make sure mass fractions are normalized and bounded between 0 and 1
 */
void domaincase::enforceMassFractions() {

    double sum;
    for(int i=0; i<domn->ngrd; i++) {
        sum = 0.0;
        for(int k=0; k<domn->gas->nSpecies(); k++) {
            if(domn->ysp[k]->d.at(i) < 0.0) domn->ysp[k]->d.at(i) = 0.0;
            if(domn->ysp[k]->d.at(i) > 1.0) domn->ysp[k]->d.at(i) = 1.0;
            sum += domn->ysp[k]->d.at(i);
        }
        for(int k=0; k<domn->gas->nSpecies(); k++)
            domn->ysp[k]->d.at(i) /= sum;
    }
}

////////////////////////////////////////////////////////////////////////////////
/** Make sure mass fractions are normalized and bounded between 0 and 1
 */
void domaincase::enforceSootMom() {

    if (!domn->pram->Lsoot)
        return;

    double sum;
    for(int i=0; i<domn->ngrd; i++) {
        for(int k=0; k<domn->pram->nsvar; k++)
            if(domn->svar[k]->d[i] < 0.0)
                domn->svar[k]->d[i] = 0.0;
    }
}
