
#include "micromixer.h"
#include "domain.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric> //accumulate

#include "interp_linear.h"

///////////////////////////////////////////////////////////////////////////////
/** micromixer constructor function
 */

micromixer::micromixer() {
    cvode = new cvodeDriver();
    nsteps = 0;
}

///////////////////////////////////////////////////////////////////////////////
/** micromixer initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */

void micromixer::init(domain *p_domn) {
    domn = p_domn;

    bool LincludeRhsMix = (domn->pram->Lsolver=="SEMI-IMPLICIT") ? true : false;
    cvode->init(domn, LincludeRhsMix);
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
  */

void micromixer::advanceOdt(const double p_tstart, const double p_tend, const int iLevel) { // iLevel is for hips

    tstart = p_tstart;
    tend   = p_tend;

    setNominalStepSize();
    if(domn->pram->LdoDL) do_DL("init");

    for(time=tstart; time<tend; time+=dt, nsteps++) {

        if(adaptGridIfNeeded() || LforceSetNominalStepSize)
           setNominalStepSize();                             // resets LforceSetNominalStepSize
        setStepSize();
        if(domn->pram->Lsolver=="EXPLICIT")
            advanceOdtSingleStep_Explicit();
        else if(domn->pram->Lsolver=="SEMI-IMPLICIT")
            advanceOdtSingleStep_SemiImplicit();
        else if(domn->pram->Lsolver=="STRANG")
            advanceOdtSingleStep_StrangSplit();

        domn->io->dumpDomainIfNeeded();
        domn->io->outputFlmltProgress();

    }

    if(domn->pram->LdoDL) do_DL("calc a");

}
///////////////////////////////////////////////////////////////////////////////
/**Set the cell sizes vectors: dxc and dx */

void micromixer::setGridDxcDx() {

    domn->mesher->setGridDxc(domn, dxc, domn->pram->cCoord);
    domn->mesher->setGridDx(domn, dx);

}

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setNominalStepSize() {

    if (domn->pram->Lspatial) {
        double velMin = *min_element(domn->uvel->d.begin(), domn->uvel->d.end());
        if (velMin <= 0.0) {
            cout << "\nError micromixer::setNominalStepSize: velMin = " << velMin << ": neg or 0" << endl;
            exit(0);
        }
    }

    domn->mesher->setGridDx(domn, dx);

    double coef = 0.0;
    double dmb;
    for (int i=0; i < domn->ngrd; i++) {
        dmb = domn->dvisc->d.at(i) / domn->rho->d.at(i) / dx.at(i) / dx.at(i);
        if (domn->pram->Lspatial)
            dmb /= domn->uvel->d.at(i);
        if (dmb > coef)
            coef = dmb;
    }
    dtStepNominal = domn->pram->diffCFL * 0.5 / coef;

    LforceSetNominalStepSize = false;
}

///////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setStepSize() {

    dt = dtStepNominal;

    if (time+dt > domn->io->dumpTimes.at(domn->io->iNextDumpTime)) {
        dt = domn->io->dumpTimes.at(domn->io->iNextDumpTime) - time;
        domn->io->LdoDump = true;
    }

    if (time + dt > tend) {
        dt = (tend - time)*(1.0+1.0E-8);
        domn->io->LdoDump = false;
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
 */

void micromixer::advanceOdtSingleStep_Explicit(){

    setGridDxcDx();
    setGf();
    if(domn->pram->LdoDL) do_DL("set DL_1");

    domn->domc->setCaseSpecificVars();

    set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    domn->v[0]->resetSourceFlags();             // sets L_source_done = false for all transported vars
    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported) {
            domn->v.at(k)->getRhsMix(gf, dxc);
            domn->v.at(k)->getRhsSrc();
        }

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + dt*( domn->v.at(k)->rhsMix.at(i) + domn->v.at(k)->rhsSrc.at(i));

    updateGrid();            // update cell sizes due to rho or rho*v variations (continuity)

    if(domn->pram->LdoDL) do_DL("set DL_2");

    domn->mesher->enforceDomainSize();     // chop the domain

}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_SemiImplicit() {

    if(domn->pram->Lsolver!="SEMI-IMPLICIT")
        return;

    setGridDxcDx();
    setGf();
    domn->domc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());
    if(domn->pram->LdoDL) do_DL("set DL_1");

    //--------------- Set the explicit (mixing) terms

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    //--------------- Perform the implicit integration on each cell

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);

    //---------------

    updateGrid();            // update cell sizes due to density variations (continuity)

    if(domn->pram->LdoDL) do_DL("set DL_2");

    domn->mesher->enforceDomainSize();     // chop the domain

}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_StrangSplit() {

    if(domn->pram->Lsolver!="STRANG")
        return;

    //--------------- First step: phi_1 = phi_0 + 0.5*dt*D(phi_0)

    setGridDxcDx();
    setGf();
    domn->domc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //--------------- Second step: phi_2 = phi_1 + dt*S(phi_2)
    // (actually: dphi/dt = S(phi), with initial condition phi=phi_1. Implicit.)

    setGridDxcDx();
    //not needed: setGf();
    domn->domc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    domn->v[0]->resetSourceFlags();             // sets L_source_done = false for all transported vars
    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsSrc();

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //--------------- Third step: phi_3 = phi_2 + 0.5*dt*D(phi_2)

    setGridDxcDx();
    setGf();
    domn->domc->setCaseSpecificVars();
    set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);

    updateGrid();            // update cell sizes due to density variations (continuity)

    //-------------------------

    domn->mesher->enforceDomainSize();     // chop the domain (for strang, just at the final step)

}

///////////////////////////////////////////////////////////////////////////////
/** Set the grid factor array
 */

void micromixer::setGf(){

    gf.resize(domn->ngrdf, 0.0);

    for (int i=1, im=0; i<domn->ngrd; i++, im++) // interior
        gf.at(i) = 2.0 / (dx.at(im) + dx.at(i));

    gf.at(0)          = 2.0 / dx.at(0);                // lo boundary
    gf.at(domn->ngrd) = 2.0 / dx.at(domn->ngrd - 1);   // hi boundary
    if (domn->pram->bcType == "PERIODIC") {            // periodic
        gf.at(0)          = 2.0 / (dx.at(0) + dx.at(domn->ngrd - 1));
        gf.at(domn->ngrd) = gf.at(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Update grid for gas expansion
 */

void micromixer::updateGrid() {

    //-------------- return for fixed grid cases

    if( (domn->pram->bcType=="WALL" && (domn->pram->vBClo==0 || domn->pram->vBChi==0)) ||
        domn->pram->LisFlmlt || domn->pram->LisFlmltX ||
        domn->pram->LisHips)
        return;

    //-------------- handle wall suction/blowing case

    if(domn->pram->bcType=="WALL" && (domn->pram->vBClo !=0 || domn->pram->vBChi != 0)){

        //---------- move all faces over (assuming vBClo > 0)

        if(domn->pram->vBClo <=0) {
            cout << "\nError micromixer::updateGrid: vBClo is <= 0.0 for suction/blowing case." << endl;
            exit(0);
        }

        double dx = domn->pram->vBClo * dt;
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

    //-------------- general variable density/outflow case

    else {

        domn->rho->setVar();

        vector<double> dxc2(domn->ngrd);
        for(int i=0; i<domn->ngrd; i++)
            dxc2[i] = dxc.at(i)*(oldrho_or_rhov.at(i)/domn->rho->d.at(i));
        if(domn->pram->Lspatial)
            for(int i=0; i<domn->ngrd; i++)
                dxc2[i] /= domn->uvel->d.at(i);

        //-------------

        domn->mesher->setGridFromDxc(dxc2);
    }


}

///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

bool micromixer::adaptGridIfNeeded() {

    if (domn->pram->Lspatial && *min_element(dx.begin(), dx.end()) < 0.9*domn->pram->dxmin) {
#ifndef SILENT
        *domn->io->ostrm << endl << "#------- ADAPTING DURING DIFFUSION" << " " << domn->ngrd;
#endif
        domn->mesher->adaptGrid(0, domn->ngrd-1);
        return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
/** Processes the DL instability
 *  @param doWhat \input string that indicates what to do.
 * */

void micromixer::do_DL(string doWhat) {

    if(!domn->pram->LdoDL)
        return;

    //--------------------------------------------------
    if(doWhat == "init") {

        uDL_1     = vector<double>(domn->ngrd, 0.0);
        uDL_2     = vector<double>(domn->ngrd, 0.0);
        xDL_1     = domn->pos->d;
        xDL_2     = domn->pos->d;
        posDL_old = domn->pos->d;
        domn->aDL->d = vector<double>(domn->ngrd, 0.0);
    }
    //--------------------------------------------------
    else if(doWhat == "set DL_1") {

        xDL_1 = xDL_2;
        uDL_1 = uDL_2;
        posDL_old = domn->pos->d;

    }
    //--------------------------------------------------
    else if(doWhat == "set DL_2") {

        xDL_2.resize(domn->ngrd);
        uDL_2.resize(domn->ngrd);

        for(int i=0; i<domn->ngrd; i++) {
            uDL_2.at(i) = (domn->pos->d[i] - posDL_old.at(i)) / (domn->pram->Lspatial ? dt/domn->uvel->d[i] : dt);
            xDL_2.at(i) = 0.5*(posDL_old.at(i) + domn->pos->d[i]);
        }
    }
    //--------------------------------------------------
    else if(doWhat == "calc a") {

        vector<double> dmb;

        dmb = uDL_1;
        Linear_interp Linterp(xDL_1, dmb);
        uDL_1.resize(domn->ngrd);
        for(int i=0; i<domn->ngrd; i++)
            uDL_1.at(i)= Linterp.interp(domn->pos->d[i]);

        dmb = uDL_2;
        Linterp = Linear_interp(xDL_2, dmb);
        uDL_2.resize(domn->ngrd);
        for(int i=0; i<domn->ngrd; i++)
            uDL_2.at(i)= Linterp.interp(domn->pos->d[i]);

        for(int i=0; i<domn->ngrd; i++)
            domn->aDL->d.at(i) = (uDL_2.at(i) - uDL_1.at(i)) / (domn->pram->Lspatial ? dt/domn->uvel->d[i] : dt);
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Store the old density for continuity.
  * This is a function for generality on inheritance.
  */

void micromixer::set_oldrho_or_rhov() {
    oldrho_or_rhov = domn->rho->d;
}


#include <iomanip>
void micromixer::check_balance(int io) {

    setGridDxcDx();
    double mom = 0.0;

    for(int i=0; i<domn->ngrd; i++)
        mom += domn->rho->d[i] * domn->uvel->d[i] * (domn->uvel->d[i]-0.0) * dxc[i];

    cout << scientific;
    cout << setprecision(13);
    cout << endl << "check: " << io << " " << mom;
}

