/**
 * @file dv_posf.cc
 * Header file for class dv_posf
 */


#include "dv_posf.h"
#include "domain.h"
#include <iostream>
#include <cstdlib>
#include <numeric> //accumulate

////////////////////////////////////////////////////////////////////////////////
/*! dv_posf  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_posf::dv_posf(domain    *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrdf, 0.0);

    double dx = domn->pram->domainLength / domn->ngrd;
    d.at(0) = domn->pram->xDomainCenter - 0.5*domn->pram->domainLength;
    for(int i=1; i<domn->ngrdf; i++)
        d.at(i) = d.at(i-1) + dx;

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_posf splitCell function
 *
 * @param isplt  \input index of cell to split
 * @param nsplt  \input number of cells to split cell into
 * @param cellFaces \input original left edge, new interior faces, orig. right edge.
 *
 * note, the number of cells is changed in the calling function, not here
 */

void dv_posf::splitCell(const int isplt,
                        const int nsplt,
                        const vector<double> &cellFaces) {

    d.insert( d.begin() + isplt+1, nsplt, 0.0 );
    for(int i=isplt+1, j=1; i<=isplt+nsplt; i++, j++)
        d.at(i) = cellFaces.at(j);

}
////////////////////////////////////////////////////////////////////////////////
/*! dv_posf merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input Thermo says an adiabatic const P mixing will change volume,
 *            but sometimes we want to retain the old volume (e.g. when we merge a small cell
 *            on the edge of the domain when enforcing the boundaries. In that case, we
 *            will chop or extend the cell anyway, so there is no conservation issue.
 */

void dv_posf::merge2cells(const int    imrg,
                          const double m1,
                          const double m2,
                          const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

    double invC = 1.0/domn->pram->cCoord;
    double C    = domn->pram->cCoord;

    vector<double> dxc;
    domn->mesher->setGridDxc(domn, dxc, domn->pram->cCoord);
    dxc.at(imrg) = (m1+m2)/domn->rho->d.at(imrg);
    if(domn->pram->Lspatial)
        dxc.at(imrg) /= domn->uvel->d.at(imrg);


    domn->mesher->setGridFromDxc(dxc);       // does pos also --> pos is done for each merge and also at the end in meshManager::merge2cells

    ////----------- outflow boundary on both sides.

    //if(domn->pram->bcType=="OUTFLOW") {
    //
    //    double V2tot = accumulate(dxc.begin(), dxc.end(), 0.0);
    //    double dmb;
    //    d[0] = -pow(0.5*V2tot, invC);
    //    for(int ie=1, iw=0, i=0; ie<d.size(); ie++, iw++, i++) {
    //        if(d[iw] <= 0.0) {
    //            dmb = pow(abs(d[iw]),C) - dxc[i];
    //            if(dmb >=0) d[ie] = -pow( dmb, invC);
    //            else        d[ie] =  pow(-dmb, invC);
    //        }
    //        else {
    //            d[ie] = pow(pow(d[iw],C) + dxc[i], invC);
    //        }
    //    }
    //}

    ////------------------ Wall on left, outlet on right. Assume all posf values are >= 0.0, expand to the right.

    //else if(domn->pram->bcType == "Wall_OUT") {
    //    d[0] = 0.0;
    //    for(int ie=1, iw=0, i=0; ie<domn->ngrdf; ie++, iw++, i++)
    //        d[ie] = pow(pow(d[iw],C) + dxc[i], invC);
    //}

    ////------------------

    //else {
    //    cout << endl << "ERROR: dv_posf::merge2cells: not setup for given bcType" << endl;
    //    exit(0);
    //}

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_posf merger2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input Thermo says an adiabatic const P mixing will change volume,
 *            but sometimes we want to retain the old volume (e.g. when we merge a small cell
 *            on the edge of the domain when enforcing the boundaries. In that case, we
 *            will chop or extend the cell anyway, so there is no conservation issue.

void dv_posf::merge2cells(const int    imrg,
                          const double m1,
                          const double m2,
                          const bool   LconstVolume) {

    d.erase(d.begin() + imrg+1);

    if(LconstVolume || domn->pram->bcType=="WALL")    //todo: generalize this (works for constant density flows only, and not spatial (due to velocity)).
        return;

    vector<double> dxc;
    domn->mesher->setGridDxc(domn, dxc);

    double invC = 1.0/domn->pram->cCoord;
    double C    = domn->pram->cCoord;

    //----------- outflow boundary on both sides.

    if(domn->pram->bcType=="OUTFLOW") {

        double dmb;
        double pm1;

        //----------- do cell imrg faces

        double xc;                   // position in cell imrg with half cell vol on each side
        if(d.at(imrg) >= 0.0) {         // all on right (before expand)
            xc = pow(pow(d.at(imrg),C)+0.5*dxc.at(imrg),invC);      // xc is found using the original cell vol.
            dxc.at(imrg) = domn->pram->Lspatial ? (m1+m2)/domn->rho->d.at(imrg)/domn->uvel->d.at(imrg)  : (m1+m2)/domn->rho->d.at(imrg); // now we update the vol to compute faces from xc
            d.at(imrg+1) = pow(pow(xc,C)+0.5*dxc.at(imrg),invC);
            dmb       = pow(d.at(imrg+1),C)-dxc.at(imrg);
            pm1       = dmb<0.0 ? -1.0 : 1.0;   // (may cross after expand)
            d.at(imrg)   = pm1*pow(pm1*dmb,invC);
        }
        else if(d.at(imrg+1) <= 0.0) {  // all on left (before expand)
            xc = -pow(pow(abs(d.at(imrg+1)),C)+0.5*dxc.at(imrg),invC);
            dxc.at(imrg) = domn->pram->Lspatial ? (m1+m2)/domn->rho->d.at(imrg)/domn->uvel->d.at(imrg)  : (m1+m2)/domn->rho->d.at(imrg);
            d.at(imrg)   = -pow(pow(abs(xc),C)+0.5*dxc.at(imrg),invC);
            dmb       = pow(abs(d.at(imrg)),C) - dxc.at(imrg);
            pm1       = dmb<0.0 ? -1.0 : 1.0;   // (may cross after expand)
            d.at(imrg+1) = -pm1*pow(pm1*dmb,invC);

        }
        else {                       // cell splits center (before expand)
            dmb = pow(d.at(imrg+1),C) - 0.5*dxc.at(imrg);
            pm1 = dmb<0.0 ? -1.0 : 1.0;
            xc  = pm1*pow(pm1*dmb,invC);
            dxc.at(imrg) = domn->pram->Lspatial ? (m1+m2)/domn->rho->d.at(imrg)/domn->uvel->d.at(imrg)  : (m1+m2)/domn->rho->d.at(imrg);
            if(xc >= 0.0) {
                d.at(imrg+1) = pow(pow(xc,C)+0.5*dxc.at(imrg),invC);
                dmb = pow(d.at(imrg+1),C) - dxc.at(imrg);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(imrg) = pm1*pow(pm1*dmb,invC);
            }
            else {
                d.at(imrg) = -pow(pow(abs(xc),C)+0.5*dxc.at(imrg),invC);
                dmb = pow(abs(d.at(imrg)),C)-dxc.at(imrg);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(imrg+1) = -pm1*pow(pm1*dmb,invC);
            }
        }

        //----------- work right: imrg+1 to end

        for(int ie=imrg+2, iw=imrg+1, i=imrg+1; ie<d.size(); ie++, iw++, i++) {
            if(d.at(iw) > 0.0)
                d.at(ie) = pow(pow(d.at(iw),C)+dxc.at(i),invC);
            else {
                dmb = pow(abs(d.at(iw)),C)-dxc.at(i);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(ie) = -pm1*pow(pm1*dmb,invC);
            }
        }

        //----------- work left: imrg-1 to 0

        for(int iw=imrg-1, i=imrg-1, ie=imrg; iw>=0; iw--, i--, ie--) {
            if(d.at(ie) < 0.0)
                d.at(iw) = -pow(pow(abs(d.at(ie)),C)+dxc.at(i),invC);
            else {
                dmb = pow(d.at(ie),C) - dxc.at(i);
                pm1 = dmb<0.0 ? -1.0 : 1.0;
                d.at(iw) = pm1*pow(pm1*dmb,invC);
            }
        }
    }

    //------------------ Wall on left, outlet on right. Assume all posf values are >= 0.0, expand to the right.

    else if(domn->pram->bcType == "Wall_OUT") {
        dxc.at(imrg) = (m1+m2)/domn->rho->d.at(imrg);
        for(int iw=imrg, ie=imrg+1, i=imrg; ie<d.size(); ie++, iw++, i++)
            d.at(ie) = pow(pow(d.at(iw),C)+dxc.at(imrg),invC);
    }

    //------------------

    else {
        cout << endl << "ERROR: dv_posf::merge2cells: not setup for given bcType" << endl;
        exit(0);
    }

}
 */

////////////////////////////////////////////////////////////////////////////////
/*! Set data array from region of domain.
 *  @param i1 \input index of starting cell of domn to build from
 *  @param i2 \input index of ending cell of domn to build from
 *  See domain::setDomainFromRegion for additional details.
 */

void dv_posf::setDvFromRegion(const int i1, const int i2){

    // note, we are owned by the eddyline, so domn is eddl, so to get domn data, use domn->domn
    const vector<double> &domn_data = domn->domn->varMap.find(var_name)->second->d;

    if(i2 >= i1)
        d.assign(domn_data.begin()+i1, domn_data.begin()+i2+2  );
    else {           // wrap around (periodic assignment)
        d.assign(domn_data.begin()+i1, domn_data.end()-1);
        double idmb  = d.size();
        d.insert(d.end(), domn_data.begin(), domn_data.begin()+i2+2 );
        for(int i=idmb; i<d.size(); i++)
            d.at(i)+=domn->domn->Ldomain();
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Resize data
 */

void dv_posf::resize() {
    d.resize(domn->ngrdf);
}













