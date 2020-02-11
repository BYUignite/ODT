/**
 * @file dv_sca.cc
 * Header file for class dv_sca
 */


#include "dv_sca.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_sca  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_sca::dv_sca(domain  *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    constantSource = domn->io->dvParams["scalarSource"] ? domn->io->dvParams["scalarSource"].as<double>() : 0.0;

}

////////////////////////////////////////////////////////////////////////////////
/*! Set the passive scalar
 *  @param ipt \input optional point to set at
 *
 * NOTE: The passive scalar evolves by mixing and stirring through the
 *  velocity field. The scalar distribution is a result of the flow and
 *  therefore does not require a setVar() method. Future changes, e.g.
 *  for the point-wise array manipulation, may be inserted here.
 */

//dv_sca::setVar(const int ipt) {}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 */

void dv_sca::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    //-------------------------

    // Volume forcing like in Kawamura's heated channel here:

    rhsSrc.resize(domn->ngrd, constantSource);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void dv_sca::getRhsMix(const vector<double> &gf,
                       const vector<double> &dxc){

    rhsMix.resize(domn->ngrd, 0.0);

    //------------------ Set fluxes

    flux.resize(domn->ngrdf);
    vector<double> rho_f(domn->ngrdf);

    interpVarToFacesHarmonic(domn->rho->d, rho_f);

    //---------- Interior faces

    for (int i=1, im=0; i<domn->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * rho_f.at(i)*domn->pram->sdiff0 * (d.at(i) - d.at(im)); // rho*sdiff0 is the "dynamic" diffusion coefficient

    //---------- Boundary faces

    if(domn->pram->bcType=="OUTFLOW") {
        flux.at(0) = 0.0;
        flux.at(domn->ngrd) = 0.0;
    }
    else if(domn->pram->bcType=="WALL") {
        double bclo = domn->pram->sBClo;
        double bchi = domn->pram->sBChi;
        flux.at(0)          = -gf.at(0)          * rho_f.at(0)          * domn->pram->sdiff0          * (d.at(0) - bclo);
        flux.at(domn->ngrd) = -gf.at(domn->ngrd) * rho_f.at(domn->ngrd) * domn->pram->sdiff0          * (bchi    - d.at(domn->ngrd-1));
    }
    else if(domn->pram->bcType=="WALL_OUT") {
        double bclo = domn->pram->sBClo;
        flux.at(0)          = -gf.at(0) * rho_f.at(0) * domn->pram->sdiff0 * (d.at(0) - bclo);
        flux.at(domn->ngrd) = 0.0;
    }
    else if(domn->pram->bcType=="PERIODIC") {
        flux.at(0)          = -gf.at(0) * rho_f.at(0) * domn->pram->sdiff0 * (d.at(0) - d.at(domn->ngrd-1));
        flux.at(domn->ngrd) = flux.at(0);
    }
    else {
        *domn->io->ostrm << endl << "ERROR: bcType not recognized in dv_sca::getRhsMix" << endl;
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

