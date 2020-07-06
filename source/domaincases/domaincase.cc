/**
 * @file domaincase.cc
 * @brief Source file for class domaincase
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

