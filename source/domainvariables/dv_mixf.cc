/**
 * @file dv_mixf.cc
 * @brief Source file for class dv_mixf
 */

#include "dv_mixf.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_mixf  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_mixf::dv_mixf(domain     *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    Dmf           = domn->io->streamProps["Dmf"] ? domn->io->streamProps["Dmf"].as<double>() : 0.0;
    if(L_transported && Dmf == 0.0) {
        cout << endl << "ERROR: if you are transporting mixture fraction, you need to set Dmf";
        exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Set mixture fraction from the gas state
 *  @param ipt \input optional point to compute at
 */

void dv_mixf::setVar(const int ipt){

    if(L_transported)
        return;        // don't set mixf from other quantities if its transported

    if(ipt != -1) {
        cout << endl << "ERROR in setVar: ipt must be -1" << endl;
        exit(0);
    }

    d.resize(domn->ngrd);

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    for(int i=0; i<domn->ngrd; i++) {
        for(int k=0; k<nsp; k++)
            yi.at(k) = domn->ysp[k]->d.at(i);
        d.at(i) = domn->strm->getMixtureFraction(&yi.at(0));
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void dv_mixf::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);
}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void dv_mixf::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(domn->ngrd, 0.0);

    setFlux(gf, dxc);

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<domn->ngrd; i++, ip++)
       rhsMix.at(i) = -domn->pram->cCoord / (domn->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(domn->posf->d.at(ip)), domn->pram->cCoord-1) -
                     flux.at(i)  * pow(abs(domn->posf->d.at(i) ), domn->pram->cCoord-1));

    if(domn->pram->Lspatial)
        for(int i=0; i<domn->ngrd; i++)
            rhsMix.at(i) /= domn->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv set face fluxes
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void dv_mixf::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    flux.resize(domn->ngrdf);
    vector<double> rhoD(domn->ngrd, Dmf);
    vector<double> rhoD_f(domn->ngrdf);

    for(int i=0; i<domn->ngrd; i++)
        rhoD.at(i) *= domn->rho->d.at(i);
    for(int i=0; i<domn->ngrdf; i++)
        rhoD_f.at(i) = linearInterpToFace(i, rhoD);

    //==========================

    //---------- Interior faces

    for (int i=1, im=0, sum=0.0; i < domn->ngrd; i++, im++)
        flux.at(i) = -rhoD_f.at(i)*gf.at(i)*(d.at(i) - d.at(im));

    //---------- Boundary faces

    flux.at(0)          = 0.0;
    flux.at(domn->ngrd) = 0.0;
}

