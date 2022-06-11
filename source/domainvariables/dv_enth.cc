/**
 * @file dv_enth.cc
 * @brief Source file for class dv_enth
 */


#include "dv_enth.h"
#include "domain.h"
#include "radiation/rad_opthin.h"
#include "radiation/rad_twoflux.h"
#include "radiation/rad_fvdom.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_enth  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_enth::dv_enth(domain    *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {


    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);
    nspc          = domn->gas->thermo()->nSpecies();

    LdoSpeciesFlux = domn->io->dvParams["LdoSpeciesFlux"] ? domn->io->dvParams["LdoSpeciesFlux"].as<bool>() : true;

    rad = 0;

    if(domn->pram->Lrad) {
        if (domn->pram->radSolveType == "OPTHIN")
            rad = new rad_opthin(domn);
        else if (domn->pram->radSolveType == "TWOFLUX")
            rad = new rad_twoflux(domn);
        else if (domn->pram->radSolveType == "FVDOM")
            rad = new rad_fvdom(domn);
        else {
            *domn->io->ostrm << endl << "ERROR: radSolveType not recognized or not consistent with cCoord" << endl;
            exit(0);
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void dv_enth::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    if(LagSrc)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    //----------- Compute the radiative source

    if(domn->pram->Lrad) {

        vector<vector<double> > xMoleSp(domn->ngrd, vector<double>(nspc));

        for(int i=0; i<domn->ngrd; i++){
            try {
                domn->domc->setGasStateAtPt(i);
            } catch (const odtCanteraError& e) {
                throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
            }
            domn->gas->thermo()->getMoleFractions(&xMoleSp.at(i).at(0));
        }

        rad->getRadHeatSource(xMoleSp, domn->temp->d, domn->pram->pres, rhsSrc);
    }

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);

    LagSrc = true;      // this is reset to false in setCaseSpecificVars


}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void dv_enth::getRhsMix(const vector<double> &gf,
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

void dv_enth::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){


    flux.resize(domn->ngrdf, 0.0);

    vector<double>          tcond(domn->ngrd, 0.0);
    vector<double>          tcond_f(domn->ngrdf, 0.0);
    vector<vector<double> > hsp(nspc, vector<double>(domn->ngrd));
    vector<double>          hh(nspc);

    for(int i=0; i<domn->ngrd; i++) {
        try {
            domn->domc->setGasStateAtPt(i);
        } catch (const odtCanteraError& e) {
            throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
        }
        tcond.at(i) = domn->gas->transport()->thermalConductivity();    // W/m*K
        if(LdoSpeciesFlux) {
            domn->gas->thermo()->getEnthalpy_RT(&hh.at(0));               // non-dimensional enthalpy
            for (int k=0; k<nspc; k++)
                hsp.at(k).at(i) = hh.at(k) * domn->temp->d.at(i) * GasConstant / domn->gas->thermo()->molecularWeight(k);        // J/kg
        }
    }

    interpVarToFacesHarmonic(tcond, tcond_f);

    //========== Do the thermal conduction portion of the flux

    //---------- Interior faces

    for (int i=1, im=0; i < domn->ngrd; i++, im++)
        flux.at(i) = -gf.at(i) * tcond_f.at(i)*(domn->temp->d.at(i) - domn->temp->d.at(im));

    //---------- Boundary faces

    if(domn->pram->bcType=="OUTFLOW") {
        flux.at(0)          = 0.0;
        flux.at(domn->ngrd) = 0.0;
    }
    else if(domn->pram->bcType=="WALL") {
        if(domn->pram->hWallBCtype=="ADIABATIC") {
            flux.at(0)          = 0.0;
            flux.at(domn->ngrd) = 0.0;
        }
        else if(domn->pram->hWallBCtype=="ISOTHERMAL") {
            flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (domn->temp->d.at(0)  - domn->pram->TBClo);
            flux.at(domn->ngrd) = -gf.at(domn->ngrd) * tcond_f.at(domn->ngrd) * (domn->pram->TBChi - domn->temp->d.at(domn->ngrd-1));
        }
        else {
            cout << endl << "ERROR: hWallBCtype unknown" << endl;
            exit(0);
        }
    }
    else if(domn->pram->bcType=="WALL_OUT") {
        if(domn->pram->hWallBCtype=="ADIABATIC") {
            flux.at(0)          = 0.0;
            flux.at(domn->ngrd) = 0.0;
        }
        else if(domn->pram->hWallBCtype=="ISOTHERMAL") {
            flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (domn->temp->d.at(0) - domn->pram->TBClo);
            flux.at(domn->ngrd) = 0.0;
        }
        else {
            cout << endl << "ERROR: hWallBCtype unknown" << endl;
            exit(0);
        }
    }
    else if(domn->pram->bcType=="PERIODIC") {
        int im = domn->ngrd - 1;
        flux.at(0)          = -gf.at(0)          * tcond_f.at(0)          * (domn->temp->d.at(0) - domn->temp->d.at(im));
        flux.at(domn->ngrd) = flux.at(0);
    }
    else {
        *domn->io->ostrm << endl << "ERROR: bcType not recognized in dv_enth::setFlux" << endl;
        exit(0);
    }

    //========== Add in the species flux portion.

    if(LdoSpeciesFlux) {

        double h_f;
        double hjsum;
        for(int i=0; i<domn->ngrdf; i++){
            hjsum = 0.0;
            for(int k=0; k<nspc; k++) {
                h_f   =  linearInterpToFace(i, hsp.at(k));
                hjsum += h_f * domn->ysp[k]->flux.at(i);
            }
            flux.at(i) += hjsum;
        }
    }

}