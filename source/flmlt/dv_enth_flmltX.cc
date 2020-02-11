/**
 * @file dv_enth_flmltX.cc
 * Definitions file for class dv_enth_flmltX
 */

#include "dv_enth_flmltX.h"
#include "domain.h"
#include <cstdlib>              // exit

////////////////////////////////////////////////////////////////////////////////
/* Constructor function
 * 
 */
dv_enth_flmltX::dv_enth_flmltX(domain            *line,
                               const string       s,
                               const bool         Lt,
                               const bool         Lo) : dv_enth(line, s, Lt, Lo){


    if(domn->pram->heatloss != 0.0){      // nonadiabatic cases, set interpolators for enthalpy profile,
        vector<double> mixf, ha, hs;      //    see setVar, below
        domn->io->read_h_mixf_flmlt_profile(mixf, ha, hs);
        Linterp_a = Linear_interp(mixf, ha);
        Linterp_s = Linear_interp(mixf, hs);
    }

}

////////////////////////////////////////////////////////////////////////////////
/* Set sensible enthalpy: ha(Ta, Ya) - h(Tmix, Ya)
 * Here, we treat ha as enthalpy. This function should be called by adiabatic 
 *   cases, otherwise it results in the sensible enthalpy between the current and "cold" states.
 * Cold state is taken as Tmix = T0*(1-mixf) + T1*(mixf), which is equal to hmix for const cp.
 * @param hsens \inout array of sensible enthalpies, for use in computing h from heatloss
 */

void dv_enth_flmltX::set_hsens(vector<double> &hsens){

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
        Tmix = domn->strm->T0*(1.0-domn->mixf->d[i]) + domn->strm->T1*domn->mixf->d[i];
        domn->gas->setState_TPY(Tmix, domn->pram->pres, &yi[0]);
        hcold = domn->gas->enthalpy_mass();

        hsens[i] = hhot - hcold;
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! dv_enth_flmlt setVar function
 *  @param ipt \input optional point to compute at
 *  Note, if we are transporting enthalpy, then we return, which should be the 
 *  case for adiabatic flamelets. For non-adiabatic flamelets, calculate enthalpy
 *  based on adiabatic profile and heatloss paramater:
 *      heatLoss = (ha(Z) - h(Z))/hs(Z) -->
 *      h(Z) = ha(Z) - heatLoss*hs(Z)
 *      Z is mixture fraction
 */

void dv_enth_flmltX::setVar(const int ipt){

    if(L_transported)
        return;

    if(domn->pram->heatloss == 0.0){
        cout << endl << "ERROR: dv_enth_flmltX::setVar should not reach this point if adiabatic; \n"
             << "   the Linterp_a and Linterp_s wont be set, etc." << endl;
        exit(0);
    }

    d.resize(domn->ngrd);

    double ha, hs;

    if(ipt == -1)
        for(int i=0; i<domn->ngrd; i++) {
            ha = Linterp_a.interp(domn->mixf->d[i]);
            hs = Linterp_s.interp(domn->mixf->d[i]);
            d[i] = ha - domn->pram->heatloss * hs;
        }
    else {
        ha = Linterp_a.interp(domn->mixf->d[ipt]);
        hs = Linterp_s.interp(domn->mixf->d[ipt]);
        d[ipt] = ha - domn->pram->heatloss * hs;
    }

}


