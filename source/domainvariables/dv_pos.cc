/**
 * @file dv_pos.cc
 * @brief Source file for class dv_pos
 */

#include "dv_pos.h"
#include "domain.h"
#include <cstdlib>

////////////////////////////////////////////////////////////////////////////////
/*! dv_pos  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_pos::dv_pos(domain    *line,
               const      string s,
               const bool Lt,
               const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    double dx = domn->pram->domainLength / domn->ngrd;
    d.at(0) = domn->pram->xDomainCenter + 0.5*(-domn->pram->domainLength + dx);
    for(int i=1; i<domn->ngrd; i++)
        d.at(i) = d.at(i-1) + dx;

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_pos splitCell function
 *
 * @param isplt  \input index of cell to split
 * @param nsplt  \input number of cells to split cell into
 * @param cellFaces \input original left edge, new interior faces, orig. right edge.
 *
 * note, the number of cells is changed in the calling function, not here
 * note, this is organized to be independent of whether we split posf before or after.
 * (that is, we often compute pos1 as simply 0.5*(posf1+posf2), but not here.
 */

void dv_pos::splitCell(const int isplt,
                       const int nsplt,
                       const vector<double> &cellFaces) {

    d.insert( d.begin() +isplt+1, nsplt, 0.0);
    for(int i=isplt, j=0; i<=isplt+nsplt; i++,j++)
        d.at(i) = 0.5*(cellFaces.at(j)+cellFaces.at(j+1));

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_pos merge2cells function
 *
 * Function presumes that the variable being merged is a quantity per unit mass.
 * Merging conservatively: (rho*phi*dx)_merged = (rho*phi*dx)_1 + (rho*phi*dx)_2
 * Then solve for phi_merged.
 *
 * NOTE: this should only be called when posf is correct.
 *
 * @param imrg    \input merge cells imrg and imrg+1
 * @param imrg \input merge cells imrg and imrg+1
 * @param m1   \input mass in cell imrg
 * @param m2   \input mass in cell imrg
 * @param LconstVolume \input (for posf, default is false)
 */

void dv_pos::merge2cells(const int    imrg,
                         const double m1,
                         const double m2,
                         const bool   LconstVolume) {

    setVar();
}

////////////////////////////////////////////////////////////////////////////////
/*! dv_pos setVar function
 *  @param ipt \input optional point to compute at
 *  Sets the grid position from posf
 *  So, you should make sure posf is consistent with pos before calling.
 *
 *  NOTE: this should only be called when posf is correct.
 */

void dv_pos::setVar(const int ipt){

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be = -1" << endl;
        exit(0);
    }

    d.resize(domn->posf->d.size()-1);

    for(int i=0; i<d.size(); i++)
        d.at(i) = 0.5*(domn->posf->d.at(i) + domn->posf->d.at(i+1));

}

////////////////////////////////////////////////////////////////////////////////
/*! Set data array from region of domain.
 *  @param i1 \input index of starting cell of domn to build from
 *  @param i2 \input index of ending cell of domn to build from
 *  See domain::setDomainFromRegion for additional details.
 */

void dv_pos::setDvFromRegion(const int i1, const int i2){

    // note, we are owned by the eddyline, so domn is eddl, so to get domn data, use domn->domn
    const vector<double> &domn_data = domn->domn->varMap.find(var_name)->second->d;

    if(i2 >= i1)
        d.assign(domn_data.begin()+i1, domn_data.begin()+i2+1  );
    else {           // wrap around (periodic assignment)
        d.assign(domn_data.begin()+i1, domn_data.end());
        d.insert(d.end(), domn_data.begin(), domn_data.begin()+i2+1 );
        int idmb = d.size();
        for(int i=idmb; i<d.size(); i++)
            d.at(i)+=domn->domn->Ldomain();
    }

}
