/**
 * @file dv_enth_flmlt.cc
 * Header file for class dv_enth_flmlt
 */

#include "dv_enth_flmlt.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv_enth_flmlt  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_enth_flmlt::dv_enth_flmlt(domain    *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, 0.0);

    if(domn->pram->heatloss != 0.0){      // nonadiabatic cases, set interpolators for enthalpy profile,
        vector<double> mixf, ha, hs;      //    see setVar, below
        domn->io->read_h_mixf_flmlt_profile(mixf, ha, hs);
        Linterp_s = Linear_interp(mixf, hs);
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *  todo: add in pressure term: unsteady, and nonuniform.
 */

void dv_enth_flmlt::getRhsSrc(const int ipt){

    if(!L_transported)
        return;

    if(LagSrc)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    LagSrc = true;      // this is reset to false in setCaseSpecificVars

}

////////////////////////////////////////////////////////////////////////////////
/*! dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void dv_enth_flmlt::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported)
        return;

    rhsMix.resize(domn->ngrd, 0.0);
}

////////////////////////////////////////////////////////////////////////////////
/*! dv_enth_flmlt setVar function
 *  @param ipt \input optional point to compute at
 *  heatLoss = (ha(Z) - h(Z))/hs(Z)
 *  h(Z) = ha(Z) - heatLoss*hs(Z)
 *  Z is mixture fraction
 *  ha is the adiabatic enthalpy profiles (linear between boundaries).
 *  hs is the sensible enthalpy profile = ha(Z) - h(Y_k,a(Z), Tmix(Z))
 *  Y_k,a(Z) is the adiabatic species profile.
 *  Tmix(Z) is a mixing temperature. This can be linear between the boundaries,
 *     but it should have the boundary temperature values so that hs = 0 on the boundaries.
 */

void dv_enth_flmlt::setVar(const int ipt){

    d.resize(domn->ngrd);

    //------------- adiabatic case

    if(domn->pram->heatloss == 0.0){   
        if(ipt == -1)
            for(int i=0; i<domn->ngrd; i++)
                d[i]= domn->strm->h0 * (1.0-domn->pos->d.at(i)) + domn->strm->h1 * (domn->pos->d.at(i));
        else
            d[ipt]= domn->strm->h0 * (1.0-domn->pos->d.at(ipt)) + domn->strm->h1 * (domn->pos->d.at(ipt));
    }
    //------------- nonadiabatic case
    else {            
        double ha, hs; 

        if(ipt == -1)
            for(int i=0; i<domn->ngrd; i++){
                ha = domn->strm->h0 * (1.0-domn->pos->d.at(i)) + domn->strm->h1 * (domn->pos->d.at(i));
                hs = Linterp_s.interp(domn->pos->d[i]);
                d[i] = ha - domn->pram->heatloss * hs;
            }
        else{
            ha = domn->strm->h0 * (1.0-domn->pos->d.at(ipt)) + domn->strm->h1 * (domn->pos->d.at(ipt));
            hs = Linterp_s.interp(domn->pos->d[ipt]);
            d[ipt] = ha - domn->pram->heatloss * hs;
        }

    }

}

////////////////////////////////////////////////////////////////////////////////
/* Set sensible enthalpy: ha(Ta, Ya) - h(Tmix, Ya)
 * Here, we treat ha as enthalpy. This function should be called by adiabatic 
 *   cases, otherwise it results in the sensible enthalpy between the current and "cold" states.
 * Cold state is taken as Tmix = T0*(1-mixf) + T1*(mixf), which is equal to hmix for const cp.
 * @param hsens \inout array of sensible enthalpies, for use in computing h from heatloss
 */

void dv_enth_flmlt::set_hsens(vector<double> &hsens){

    if(domn->pram->heatloss != 0.0){
        cout << endl << "ERROR: set_hsens should be called for adiabatic cases" << endl;
        exit(0);
    }

    hsens.resize(domn->ngrd);

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    double hhot, hcold;
    double Tmix;

    for(int i=0; i<domn->ngrd; i++){
        hhot = domn->enth->d[i];

        for(int k=0; k<nsp; k++) {
            yi.at(k) = domn->ysp[k]->d.at(i);
        }
        Tmix = domn->strm->T0*(1.0-domn->pos->d[i]) + domn->strm->T1*domn->pos->d[i];
        domn->gas->setState_TPY(Tmix, domn->pram->pres, &yi[0]);
        hcold = domn->gas->enthalpy_mass();

        hsens[i] = hhot - hcold;
    }

}
