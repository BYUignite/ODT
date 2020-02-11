
#include "micromixer_flmlt.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer_flmlt::setNominalStepSize() {

    if(domn->pram->LisFlmltX) 
        micromixer::setNominalStepSize();
    else{
        domn->mesher->setGridDx(domn, dx);

        double coef = 0.0;
        double dmb;
        for (int i=0; i < domn->ngrd; i++) {
            dmb = 0.5*domn->chi->d.at(i)/dx[i]/dx[i];
            if (dmb > coef)
                coef = dmb;
        }
        dtStepNominal = domn->pram->diffCFL * 0.5 / coef;
    }
}


///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

bool micromixer_flmlt::adaptGridIfNeeded() {

    if(!domn->pram->LletFlmltAdpt)
        return false;

    double tperiod;
    if(domn->pram->LisFlmlt) 
        tperiod = 1.0/domn->pram->chi0 * 10;
    else{
        tperiod = domn->pram->domainLength*domn->pram->domainLength*domn->pram->rho0/domn->pram->kvisc0;
    }

    if(tNextAdapt == -1.0)   // just to initialize it
        tNextAdapt = tperiod;
    if(time > tNextAdapt) {
        domn->mesher->adaptGrid(0, domn->ngrd-1);
        tNextAdapt += tperiod;
        return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
/** Set the grid factor array
 */

void micromixer_flmlt::setGf() {
    if (domn->pram->LisFlmltX)
        micromixer::setGf();
    else
        return;
}

///////////////////////////////////////////////////////////////////////////////
/**Set the cell sizes vectors: dxc and dx */

void micromixer_flmlt::setGridDxcDx() {

    if (domn->pram->LisFlmltX)
        micromixer::setGridDxcDx();
    else
        return;

}
