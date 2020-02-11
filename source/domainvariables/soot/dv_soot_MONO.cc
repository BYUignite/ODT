/**
 * @file dv_soot_MONO.cc
 * Header file for class dv_soot_MONO
 * @author Victoria B. Lansinger
 */

#include "dv_soot_MONO.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *
 *  @param ipt \input optional point to compute source at.
 *
 *  adapted from dv_ygas.cc
 */

void dv_soot_MONO::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    if(L_source_done)                                    //  return if this routine was already called (by dv_ygas)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    static vector<vector<double> > rrSvar(nsvar);        // temp storage for moment rates [nsvar][ngrd]
    static vector<double>          yi(domn->strm->nspc); // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = (domn->ngrd)-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                                       // to save cost, compute needed terms for all dv_soot_MONO objects using this one.

        for(int k=0; k<nsvar; k++)                     // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);           // resize to match grid size

        for(int k=0; k<gasSootSources.size(); k++)         // resize gasSootSources for each gas speckes k
            gasSootSources[k].resize(domn->ngrd,0.0);

        for(int i=iS; i<=iE; i++) {                    // loop over grid points

             domn->domc->setGasStateAtPt(i);         // keep commented for decoupled soot
            domn->domc->enforceSootMom();

            double M0    = domn->svar[0]->d[i] * domn->rho->d[i];       // M0 = #/m3
            double M1    = domn->svar[1]->d[i] * domn->rho->d[i];       // M1 = rhoYs = kg/m3
                                                                                                // the abs should not be needed given enforceSootMom
            //---------- set weights and abscissas

            if (M0 <= 0.0) {
                wts[0] = 0.0;
                absc[0] = 0.0;
            }
            else {
                wts[0] = M0;                 // defined weights and abscissas for the monodisperse case
                absc[0] = M1/M0;
            }

            //--------- chemical soot rates

            for(int k=0; k<yi.size(); k++)
                yi[k] = domn->ysp[k]->d[i];
            set_gas_state_vars(domn->temp->d[i],
                               domn->pram->pres,
                               domn->rho->d[i],
                               domn->gas->meanMolecularWeight(),
                               domn->dvisc->d[i],
                               yi);


            double Jnuc  = getNucleationRate(absc, wts);         // #/m3*s
            double Kgrw  = getGrowthRate(M0, M1);                // kg/m2*s
            double Koxi  = getOxidationRate(M0, M1);             // kg/m2*s
            double Coag  = getCoagulationRate(absc[0], absc[0]);

            //--------- nucleation terms

            double N0 = Jnuc;                                    // #/m3*s
            double N1 = Jnuc*Cmin*MW_c/Na;                       // kg/m3*s

            //---------- PAH condensation terms

            double Cnd0 = 0.0;
            double Cnd1 = 0.0;

            if (nucleation_mech=="PAH") {                        // condense PAH if nucleate PAH
                Cnd1 = DIMER*m_dimer*getCoagulationRate(m_dimer, absc[0])*wts[0];
            }

            //--------- growth terms

            double Am2m3 = 0.0;                                  // m^2_soot / m^3_total
            if (M0 > 0.0)
                Am2m3 = M_PI * pow(abs(6/(M_PI*rhoSoot)*M1/M0),2.0/3.0) * abs(M0);    // m^2_soot / m^3_total = pi*di^2*M0

            double G0 = 0.0;                                     // zero by definition, #/m3*s
            double G1 = Kgrw*Am2m3;                              // kg/m3*s

            //--------- oxidation terms

            double X0 = 0.0;                                     // zero by definition, #/m3*s
            double X1 = -Koxi*Am2m3;                             // kg/m3*s

            ////--------- coagulation terms

            double C0 = -0.5*Coag*wts[0]*wts[0];                 // #/m3*s
            double C1 = 0.0;                                     // zero by definition, kg/m3*s

            //--------- combine to make source terms

            rrSvar[0][i] = (N0 + Cnd0 + G0 + X0 + C0) / rho;     // (#/m3*s)/rho = #/kg*s
            rrSvar[1][i] = (N1 + Cnd1 + G1 + X1 + C1) / rho;     // (kg-soot/m3*s)/rho = kg-soot/kg*s

            //---------- compute gas source terms

            set_gasSootSources(N1, Cnd1, G1, X1, i);

        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);

    L_source_done = true;

}
