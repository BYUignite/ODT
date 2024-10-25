/**
 * @file radiation.cc
 * @brief Source file for class radiation
 */

#include "radiation.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Static members


///////////////////////////////////////////////////////////////////////////////
/** Contructor for radiation class object
 *
 *  @author Victoria B. Stephens 
 *
 *  @param p_domn \input pointer to domain object
 *
 */

radiation::radiation(domain *p_domn) {

    domn = p_domn;

    if(!domn->pram->Lrad)
        return;

    radProps = new radiationProperties();
    radProps->init(domn);

    sigmaSB = 5.670E-8;                     // W/m2*K4

}
