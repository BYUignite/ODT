/**
 * @file dv_soot_QMOM.cc
 * Header file for class dv_soot_QMOM
 * @author Victoria B. Lansinger
 */

#include "dv_soot_QMOM.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>

void pdAlg(int nm, int np, vector<double> &mu, vector<double> &wts, vector<double> &absc );
void wheeler(const vector<double> &m, int N, vector<double> &w, vector<double> &x );
void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs, int &Nout, vector<double> &w, vector<double> &x );
void adaptiveWheelerAlgorithm(const std::vector<double>& moments, std::vector<double>& w, std::vector<double>& x,const double& rMin, const double& eAbs);

////////////////////////////////////////////////////////////////////////////////
/*! getRhsSrc function
 *
 *      dv source term part of the rhs function. Uses QMOM. Modeled off of
 *      dv_ygas class.
 *
 *      @param ipt \input optional point to compute source at.
 *
 */

void dv_soot_QMOM::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    if(L_source_done)                                    //  return if this routine was already called (by dv_ygas)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);                       // resize rhsSrc

    static vector<vector<double> > rrSvar(nsvar);         // temp storage for moment rates [nsvar][ngrd]
    static vector<double>          yi(domn->strm->nspc);  // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                                            // to save cost, compute needed terms for all dv_soot_QMOM objects using this one.

        for(int k=0; k<nsvar; k++)                          // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);                // resize to match grid size

        for(int k=0; k<gasSootSources.size(); k++)          // resize gasSootSources for each gas speckes k
            gasSootSources[k].resize(domn->ngrd,0.0);

        for(int i=iS; i<=iE; i++) {                         // loop over grid points

            domn->domc->setGasStateAtPt(i);              // set gas state; uncomment to couple gas/soot
            domn->domc->enforceSootMom();

            //---------- get moments from domain

            vector<double> M(nsvar,0.0);
            for (int k=0; k<nsvar; k++)
                M[k] = domn->svar[k]->d[i] * domn->rho->d[i];

            //---------- set weights and abscissas

            getWtsAbs(M, wts, absc);                        // PD and wheeler algorithms called in here
            for(int i=0; i<wts.size();i++){
                if(wts[i] < 0.0)  wts[i]  = 0.0;
                if(absc[i] < 0.0) absc[i] = 0.0;
            }

            //---------- get chemical soot rates

            for(int k=0; k<yi.size(); k++)
                yi[k] = domn->ysp[k]->d[i];
            set_gas_state_vars(domn->temp->d[i],
                               domn->pram->pres,
                               domn->rho->d[i],
                               domn->gas->meanMolecularWeight(),
                               domn->dvisc->d[i],
                               yi);

            double Jnuc = getNucleationRate(absc, wts);             // #/m3*s
            double Kgrw = getGrowthRate(M[0], M[1]);                // kg/m2*s
            double Koxi = getOxidationRate(M[0], M[1]);             // kg/m2*s

            //---------- nucleation terms

            vector<double> Mnuc(nsvar,0.0);                         // nucleation source terms for moments
            double m_nuc = Cmin*MW_c/Na;                            // mass of nucleated particle
            for (int k=0; k<nsvar; k++)
                Mnuc[k] = pow(m_nuc,k) * Jnuc;                      // Nr = m_min^r * Jnuc

            //---------- PAH condensation terms

            vector<double> Mcnd(nsvar,0.0);                             // initialize to 0.0
            if (nucleation_mech == "PAH") {                             // condense PAH if nucleate PAH
                for (int k=1; k<nsvar; k++) {                           // Mcnd[0] = 0.0 by definition
                    for (int ii=0; ii<absc.size(); ii++)
                        Mcnd[k] += getCoagulationRate(m_dimer, absc[ii])*pow(absc[ii],k-1)*wts[ii];
                    Mcnd[k] *= DIMER*m_dimer*k;
                }
            }

            //---------- growth terms

            vector<double> Mgrw(nsvar,0.0);                           // growth source terms for moments
            double Acoef = M_PI*pow(abs(6.0/M_PI/rhoSoot),2.0/3.0);   // Acoef = kmol^2/3 / kg^2/3
            for (int k=1; k<nsvar; k++)                               // Mgrw[0] = 0.0 by definition
                Mgrw[k] = Kgrw * Acoef * k * Mk(k-1.0/3.0);           // kg^k/m3*s

            //---------- oxidation terms

            vector<double> Moxi(nsvar,0.0);
            for (int k=1; k<nsvar; k++)                               // Moxi[0] = 0.0 by definition
                Moxi[k] = -Koxi * Acoef * k * Mk(k-1.0/3.0);           // kg^k/m3*s

            //---------- coagulation terms

            vector<double> Mcoa(nsvar,0.0);                           // coagulation source terms: initialize to zero!
            for(int k=0; k<nsvar; k++) {
                if(k==1) continue;
                for(int ii=1; ii<absc.size(); ii++)        // off-diagonal terms (looping half of them) with *2 incorporated
                    for(int j=0; j<ii; j++)
                        Mcoa[k] += getCoagulationRate(absc[ii], absc[j])*wts[ii]*wts[j]* (k==0 ? -1.0 : (pow(absc[ii]+absc[j],k))-pow(absc[ii],k)-pow(absc[j],k) );
                for(int ii=0; ii<absc.size(); ii++)              // diagonal terms
                    Mcoa[k] += getCoagulationRate(absc[ii], absc[ii])*wts[ii]*wts[ii]* (k==0 ? -0.5 : pow(absc[ii],k)*(pow(2,k-1)-1) );        //(pow(absc[ii]+absc[ii],k))-2.0*pow(absc[ii],k) );
            }

            //---------- combinine to make source terms

            for (int k=0; k<nsvar; k++)
                rrSvar[k][i] = (Mnuc[k] + Mcnd[k] + Mgrw[k] + Moxi[k] + Mcoa[k]) / rho;

            //---------- compute gas source terms

            set_gasSootSources(Mnuc[1], Mcnd[1], Mgrw[1], Moxi[1], i);

        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);                               // assign to domn's rhsSrc variable

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);

    L_source_done = true;
}

////////////////////////////////////////////////////////////////////////////////
/*! Mk function
 *
 *      Calculates fractional moments from weights and abscissas.
 *
 *      @param exp  \input  fractional moment to compute, corresponds to exponent
 *      @param wts  \input  weights
 *      @param absc \input  abscissas
 *      @param Mk   \output fractional moment value
 *
 */

double dv_soot_QMOM::Mk(double exp) {

    double Mk = 0;

    for(int k=0; k<nsvar/2; k++) {
        if (wts[k] == 0 || absc[k] == 0)
            return 0;
        else
            Mk += wts[k] * pow(absc[k],exp);
    }
    return Mk;

}

////////////////////////////////////////////////////////////////////////////////
/*! getWtsAbs function
 *
 *      Calculates weights and abscissas from moments using PD algorithm or
 *      wheeler algorithm.
 *
 *      @param M        \input  vector of moments
 *      @param wts      \input  weights
 *      @param absc     \input  abscissas
 *
 *      Notes:
 *      - Use wheeler over pdalg whenever possible.
 *      - wts and abs DO NOT change size; if we downselect to a smaller number
 *      of moments, the extra values are set at and stay zero
 *      - using w_temp and a_temp means we don't have to resize wts and absc,
 *      which is more convenient when wts and absc are used to reconstitute
 *      moment source terms.
 */

void dv_soot_QMOM::getWtsAbs(vector<double> M, vector<double> &wts, vector<double> &absc) {

    for (int k=0; k<nsvar; k++) {                  // if any moments are zero, return with zero wts and absc
        if (M[k] <= 0.0)
            return;
    }

    int N = nsvar;                                 // local nsvar
    bool negs = false;                             // flag for downselecting if there are negative wts/abs

    vector<double> w_temp(N/2,0.0);
    vector<double> a_temp(N/2,0.0);

    do {                                           // downselection loop

        negs = false;                              // reset flag

        for (int k=0; k<N/2; k++) {                // reinitialize wts and abs with zeros
            w_temp[k] = 0.0;
            a_temp[k] = 0.0;
        }

        if (N == 2) {                              // in 2 moment case, return monodisperse output
           wts[0]  = M[0];
           absc[0] = M[1]/M[0];
           return;
        }

        //pdAlg(N,N/2,ymom,wts,absc);              // PD algorithm
        wheeler(M, N/2, w_temp, a_temp);           // wheeler algorithm

        for (int k=0; k<N/2; k++) {
            if (w_temp[k] < 0.0 || a_temp[k] < 0.0 || a_temp[k] > 1.0)
                negs = true;
        }

        if (negs == true) {                        // if we found negative values
            N = N - 2;                             // downselect to two fewer moments and try again
            w_temp.resize(N/2);
            a_temp.resize(N/2);
        }

    } while (negs == true);                        // end of downselection loop

    for (int k = 0; k < w_temp.size(); k++) {      // assign temporary variables to output
        wts[k]  = w_temp[k];
        absc[k] = a_temp[k];
    }

}
