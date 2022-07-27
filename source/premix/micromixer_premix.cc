
#include "micromixer_premix.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

bool micromixer_premix::adaptGridIfNeeded() {

    domn->mesher->adaptGrid(0, domn->ngrd-1);   //doldb
    return true;  //doldb

    double tperiod = domn->pram->domainLength*domn->pram->domainLength*domn->pram->rho0/domn->pram->kvisc0;

    if(tNextAdapt == -1.0)   // just to initialize it
        tNextAdapt = tperiod;
    if(time > tNextAdapt) {
        domn->mesher->adaptGrid(0, domn->ngrd-1);
        tNextAdapt += tperiod;
        return true;
    }
    else
        return false;
}

///////////////////////////////////////////////////////////////////////////////
/** Update grid for gas expansion
 */

void micromixer_premix::updateGrid() {

    /////////////// expand grid cells

    domn->rho->setVar();

    vector<double> dxc2(domn->ngrd);
    for(int i=0; i<domn->ngrd; i++)
        dxc2[i] = dxc.at(i)*(oldrho_or_rhov.at(i)/domn->rho->d.at(i));
    if(domn->pram->Lspatial)
        for(int i=0; i<domn->ngrd; i++)
            dxc2[i] /= domn->uvel->d.at(i);

    domn->mesher->setGridFromDxc(dxc2);

    /////////////// insert new cell to implement inlet velocity

    if(domn->pram->uInflow <=0)
        return;

    //---------- move all faces over (assuming vBClo > 0)

    double dx = domn->pram->uInflow * dt;
    for(int i=0; i<domn->ngrdf; i++)
        domn->posf->d[i] += dx;

    //---------- insert a cell at the beginning of domain

    for(int k=0; k < domn->v.size(); k++)
        domn->v[k]->d.insert(domn->v[k]->d.begin(), domn->domc->inlet_cell_dv_props[k]);
    domn->pos->d[0] = 0.5*(domn->posf->d[0] + domn->posf->d[1] );
    domn->ngrd++;
    domn->ngrdf++;

    domn->mesher->merge2cells(0, true);

    if((domn->posf->d[1]-domn->posf->d[0]) > 2.0*(domn->posf->d[2]-domn->posf->d[1])){
        vector<double> faces{domn->posf->d[0], 0.5*(domn->posf->d[0]+domn->posf->d[1]), domn->posf->d[1]};
        domn->mesher->splitCell(0,1,faces);
    }

    //---------- at the end of domain split cell and delete the overhang

    domn->mesher->enforceDomainSize();     // chop the domain (for strang, just at the final step)

    LforceSetNominalStepSize = true;       // we changed the grid, so indicate to recompute nominal timestep size

}