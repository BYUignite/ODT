/**
 * @file dv_uvw.cc
 * Header file for class dv_uvw
 */


#include "dv_uvw.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_uvw  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_uvw::dv_uvw(domain  *line,
               const      string s,
               const bool Lt,
               const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 */

void dv_uvw::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    //-------------------------

    if(ipt==-1) {
        if(var_name == "uvel" && domn->pram->cCoord != 3.0) {
            rhsSrc = vector<double>(domn->ngrd, -domn->pram->dPdx);
            for(int i=0; i<domn->ngrd; i++) {
                if(domn->pram->Lbuoyant)
                    rhsSrc.at(i) += (domn->rho->d.at(i) - domn->rho->d.at(domn->ngrd-1))*domn->pram->g;
                rhsSrc.at(i) /= domn->rho->d.at(i);
            }
        }

        if(domn->pram->Lspatial)
            for(int i=0; i<domn->ngrd; i++)
                rhsSrc.at(i) /= domn->uvel->d.at(i);
    }
    else {
        if(var_name == "uvel" && domn->pram->cCoord != 3.0) {
            rhsSrc.at(ipt) = -domn->pram->dPdx;
            if(domn->pram->Lbuoyant)
                rhsSrc.at(ipt) += (domn->rho->d.at(ipt) - domn->rho->d.at(domn->ngrd-1))*domn->pram->g;
            rhsSrc.at(ipt) /= domn->rho->d.at(ipt);
        }

        if(domn->pram->Lspatial)
            rhsSrc.at(ipt) /= domn->uvel->d.at(ipt);
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void dv_uvw::getRhsMix(const vector<double> &gf,
                       const vector<double> &dxc){

    rhsMix.resize(domn->ngrd, 0.0);

    //------------------ Set fluxes

    flux.resize(domn->ngrdf);
    vector<double> dvisc_f(domn->ngrdf);

    interpVarToFacesHarmonic(domn->dvisc->d, dvisc_f);

    //---------- Interior faces

    for (int i=1, im=0; i < domn->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * dvisc_f.at(i)*(d.at(i) - d.at(im));

    //---------- Boundary faces

    if(domn->pram->bcType=="OUTFLOW") {
        flux.at(0) = 0.0;
        flux.at(domn->ngrd) = 0.0;
    }
    else if(domn->pram->bcType=="WALL") {
        double bclo = var_name=="uvel" ? domn->pram->uBClo : var_name=="vvel" ? domn->pram->vBClo : domn->pram->wBClo;
        double bchi = var_name=="uvel" ? domn->pram->uBChi : var_name=="vvel" ? domn->pram->vBChi : domn->pram->wBChi;
        flux.at(0)          = -gf.at(0)          * dvisc_f.at(0)          * (d.at(0) - bclo);
        flux.at(domn->ngrd) = -gf.at(domn->ngrd) * dvisc_f.at(domn->ngrd) * (bchi - d.at(domn->ngrd-1));
    }
    else if(domn->pram->bcType=="WALL_OUT") {
        double bclo = var_name=="uvel" ? domn->pram->uBClo : var_name=="vvel" ? domn->pram->vBClo : domn->pram->wBClo;
        flux.at(0)          = -gf.at(0) * dvisc_f.at(0) * (d.at(0) - bclo);
        flux.at(domn->ngrd) = 0.0;
    }
    else if(domn->pram->bcType=="PERIODIC") {
        int im = domn->ngrd - 1;
        flux.at(0)          = -gf.at(0) * dvisc_f.at(0) * (d.at(0) - d.at(im));
        flux.at(domn->ngrd) = flux.at(0);
    }
    else {
        *domn->io->ostrm << endl << "ERROR: bcType not recognized in dv_uvw::getRhsMix" << endl;
        exit(0);
    }

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<domn->ngrd; i++, ip++)
       rhsMix.at(i) = -domn->pram->cCoord / (domn->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(domn->posf->d.at(ip)), domn->pram->cCoord-1) -
                     flux.at(i)  * pow(abs(domn->posf->d.at(i) ), domn->pram->cCoord-1));

    if(domn->pram->Lspatial)
        for(int i=0; i<domn->ngrd; i++)
            rhsMix.at(i) /= domn->uvel->d.at(i);

}


