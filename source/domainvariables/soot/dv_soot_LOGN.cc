/**
 * @file dv_soot_LOGN.cc
 * Header file for class dv_soot_LOGN
 * @author Victoria B. Lansinger
 */

#include "dv_soot_LOGN.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *
 *    copied from dv_ygas.cc
 */

void dv_soot_LOGN::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    if(L_source_done)                                    //  return if this routine was already called (by dv_ygas)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);                    // resize rhsSrc
    static vector<double>        yi(domn->strm->nspc); // [nspc]

    static vector<vector<double> > rrSvar(nsvar);      // temp storage for moment rates [nsvar][ngrd]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = (domn->ngrd)-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    for(int k=0; k<nsvar; k++)                         // for each mom/section,
        rrSvar.at(k).resize(domn->ngrd);               // resize to match grid size

    if(kMe==0) {                                       // to save cost, compute needed terms for all dv_soot_LOGN objects using this one.

        for(int k=0; k<nsvar; k++)                     // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);           // resize to match grid size

        for(int k=0; k<gasSootSources.size(); k++)         // resize gasSootSources for each gas speckes k
            gasSootSources[k].resize(domn->ngrd,0.0);

        for(int i=iS; i<=iE; i++) {                    // loop over grid points

            // domn->domc->setGasStateAtPt(i);         // keep commented for decoupled soot
            domn->domc->enforceSootMom();

            M0 = domn->svar[0]->d[i] * domn->rho->d[i];             // M0 = #/m3
            M1 = domn->svar[1]->d[i] * domn->rho->d[i];             // M1 = rhoYs = kg/m3
            M2 = domn->svar[2]->d[i] * domn->rho->d[i];             // M2 = kg2/m3

            for(int k=0; k<yi.size(); k++)
                yi[k] = domn->ysp[k]->d[i];
            set_gas_state_vars(domn->temp->d[i],
                               domn->pram->pres,
                               domn->rho->d[i],
                               domn->gas->meanMolecularWeight(),
                               domn->dvisc->d[i],
                               yi);

            //--------- nucleation and condensation terms

            double N0;                                              // #/m3*s
            double N1;                                              // kg/m3*s
            double N2;                                              // kg2/m3*s

            double Cnd0 = 0.0;                                      // by definition.
            double Cnd1;
            double Cnd2;

            double Kfm = get_Kfm();                // used in coagulation below
            double Kc  = get_Kc();                 // used below
            double Kcp = get_Kcp();                // used below

            double mmin = Cmin*MW_c/Na;
            double Jnuc;
            if(nucleation_mech !="PAH")
                Jnuc = getNucleationRate();
            else {

                //------ nucleation

                double wdotD = set_m_dimer();


                double mD  = m_dimer;
                double Ifm = Kfm*b_coag*( M0*pow(mD,1./6.) + 2*Mk(1./3.)*pow(mD,-1./6.) +
                                          Mk(2./3.)*pow(mD,-1./2.) + Mk(-1./2.)*pow(mD,2./3.) +
                                          2*Mk(-1./6.)*pow(mD,1./3.) + Mk(1./6.) );
                double Ic  = Kc*( 2*M0 + Mk(-1./3.)*pow(mD,1./3.) + Mk(1./3.)*pow(mD,-1./3.) +
                                  Kcp*( M0*pow(mD,-1./3.) + Mk(-1./3.) +
                                        Mk(1./3.)*pow(mD,-2./3.) + Mk(-2./3.)*pow(mD,1./3.)) );

                double I_beta_DS = Ic*Ifm/(Ic+Ifm);           // harmonic mean

                double beta_DD = getCoagulationRate(mD, mD);   // dimer self-collision rate
                double DIMER = 2.0*wdotD/(I_beta_DS + sqrt(I_beta_DS*I_beta_DS + 4*beta_DD*wdotD));   // #/m3

                Jnuc = 0.5*beta_DD*DIMER*DIMER;              // #/m3*s

                //------ PAH condensation

                double Ifm1 = Ifm;
                double Ifm2 = Kfm*b_coag*( M1*pow(mD,1./6.) + 2*Mk(4./3.)*pow(mD,-1./6.) +
                                           Mk(5./3.)*pow(mD,-1./2.) + Mk( 1./2.)*pow(mD,2./3.) +
                                           2*Mk( 5./6.)*pow(mD,1./3.) + Mk(7./6.) );
                double Ic1  = Ic;
                double Ic2  = Kc*( 2*M1 + Mk( 2./3.)*pow(mD,1./3.) + Mk(4./3.)*pow(mD,-1./3.) +
                                   Kcp*( M1*pow(mD,-1./3.) + Mk( 2./3.) +
                                         Mk(4./3.)*pow(mD,-2./3.) + Mk( 1./3.)*pow(mD,1./3.)) );

                Cnd1 =     mD*DIMER* (Ic1*Ifm1)/(Ic1+Ifm1);             // applying harmonic means
                Cnd2 = 2.0*mD*DIMER* (Ic2*Ifm2)/(Ic2+Ifm2);
            }
            //-----

            N0 = Jnuc;                                       // #/m3*s
            N1 = Jnuc*mmin;                                  // kg/m3*s
            N2 = Jnuc*mmin*mmin;                             // kg2/m3*s

            //--------- growth terms

            double Kgrw  = getGrowthRate(M0, M1);             // kg/m2*s

            double term = Kgrw * M_PI*pow(6.0/rhoSoot/M_PI,2.0/3.0);

            double G0 = 0.0;                                        // zero by definition, #/m3*s
            double G1 = term * Mk(2./3.);                           // kg/m3*s
            double G2 = term * Mk(5./3.) * 2;                       // kg2/m3*s

            //--------- oxidation terms

            double Koxi  = getOxidationRate(M0, M1);          // kg/m2*s

            double X0 = 0.0;                                        // zero by definition, #/m3*s
            double X1 = Koxi * M_PI*pow(6.0/rhoSoot/M_PI,2.0/3.0) * Mk(2./3.);      // kg/m3*s
            double X2 = Koxi * M_PI*pow(6.0/rhoSoot/M_PI,2.0/3.0) * Mk(5./3.) * 2;  // kg2/m3*s

            //--------- coagulation terms

            //---- free molecular
            double C0_fm = -Kfm * b_coag * (M0*Mk(1./6.) + 2.0*Mk(1./3.)*Mk(-1./6.) +
                                            Mk(2./3.)*Mk(-1./2.));       // #/m3*s
            double C1_fm = 0.0;                                                // zero by definition, kg/m3*s
            double C2_fm = 2*Kfm* b_coag * (M1*Mk(7./6.) + 2*Mk(4./3.)*Mk(5./6.) +
                                            Mk(5./3.)*Mk(1./2.));        // kg2/m3*s

            //---- continuum
            double C0_c = -Kc*( M0*M0 + Mk(1./3.)*Mk(-1./3.) + Kcp*(M0*Mk(-1./3.) + Mk(1./3.)*Mk(-2./3.)) );
            double C1_c = 0.0;
            double C2_c = 2*Kc*(M1*M1 + Mk(2./3.)*Mk(4./3.) + Kcp*(M1*Mk(2./3.) + Mk(1./3.)*Mk(4./3.)));

            //----- harmonic mean
            double C0 = C0_fm*C0_c/(C0_fm+C0_c);
            double C1 = 0.0;
            double C2 = C2_fm*C2_c/(C2_fm+C2_c);

            //--------- combinine to make source terms

            rrSvar[0][i] = (N0 + G0 + Cnd0 - X0 + C0) / domn->rho->d[i];   // (#/m3*s)/rho = #/kg*s
            rrSvar[1][i] = (N1 + G1 + Cnd1 - X1 + C1) / domn->rho->d[i];   // (kg-soot/m3*s)/rho = kg-soot/kg*s
            rrSvar[2][i] = (N2 + G2 + Cnd2 - X2 + C2) / domn->rho->d[i];   // (kg-soot/m3*s)ˆ2/rho = kg-sootˆ2/kg*s

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

////////////////////////////////////////////////////////////////////////////////
/*! Mk function
 *    Calculates fractional moments
 *
 *    @param k   \input  fractional moment to compute, corresponds to exponent
 *
 */

double dv_soot_LOGN::Mk(const double &k) {

    double M0_exp = 1 + 0.5*k*(k-3);
    double M1_exp = k*(2-k);
    double M2_exp = 0.5*k*(k-1);

    if (M2 == 0.0 && M2_exp < 0)
        M2_exp = 0;

    return pow(M0, M0_exp) * pow(M1, M1_exp) * pow(M2, M2_exp);
}

