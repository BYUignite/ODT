/**
 * @file dv_soot_MOMIC.cc
 * Header file for class dv_soot_MOMIC
 * @author Victoria B. Lansinger
 */

#include "dv_soot_MOMIC.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>

////////////////////////////////////////////////////////////////////////////////
/*! getRhsSrc function
 *
 *      dv source term part of the rhs function. Uses MOMIC. Modeled off of
 *      dv_ygas class.
 *
 *      @param ipt \input optional point to compute source at.
 *
 */

void dv_soot_MOMIC::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    if(L_source_done)                                       //  return if this routine was already called (by dv_ygas)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);                         // resize rhsSrc

    static vector<vector<double> > rrSvar(nsvar);           // temp storage for moment rates [nsvar][ngrd]
    static vector<double>          yi(domn->strm->nspc);    // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    for(int k=0; k<nsvar; k++)                              // for each mom/section,
        rrSvar.at(k).resize(domn->ngrd);                    // resize to match grid size

    if(kMe==0) {                                            // to save cost, compute needed terms for all dv_soot_MOMIC objects using this one.

        for(int k=0; k<nsvar; k++)                          // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);                // resize to match grid size

        for(int k=0; k<gasSootSources.size(); k++)          // resize gasSootSources for each gas speckes k
            gasSootSources[k].resize(domn->ngrd,0.0);

        for(int i=iS; i<=iE; i++) {                         // loop over grid points

            // domn->domc->setGasStateAtPt(i);              // set gas state; uncomment to couple gas/soot
            domn->domc->enforceSootMom();                   // make sure moments are positive or zero

            //---------- get moments from domain

            vector <double> M(nsvar,0.0);
            for (int k=0; k<nsvar; k++) {
                M[k] = domn->svar[k]->d[i] * domn->rho->d[i];
            }

            //---------- determine how many moments to use

            int N = nsvar;                                 // local number of moments
            downselectIfNeeded(i, M, N);                   // downselect() will not change anything if all moment values >0

            //---------- get chemical soot rates

            for(int k=0; k<yi.size(); k++)
                yi[k] = domn->ysp[k]->d[i];
            set_gas_state_vars(domn->temp->d[i],
                               domn->pram->pres,
                               domn->rho->d[i],
                               domn->gas->meanMolecularWeight(),
                               domn->dvisc->d[i],
                               yi);

            double Jnuc = getNucleationRate();             // #/m3*s
            double Kgrw = getGrowthRate();                 // kg/m2*s
            double Koxi = getOxidationRate();              // kg/m2*s

            //---------- calculate MOMIC source terms

            //---------- nucleation terms

            vector<double> Mnuc(nsvar,0.0);                         // nucleation source terms for moments

            double m_nuc = Cmin*MW_c/Na;                            // mass of nucleated particle
            for (int k=0; k<N; k++)
                Mnuc[k] = pow(m_nuc,k) * Jnuc;                      // Nr = m_min^r * Jnuc

            //---------- PAH condensation terms

            vector<double> Mcnd(nsvar,0.0);                         // initialize to 0.0

            if (nucleation_mech == "PAH") {                         // condense PAH if nucleate PAH
                for (int k=1; k<N; k++) {                           // Mcnd[k] = 0.0 by definition
                    Mcnd[k] = getCoag(domn->temp->d[i],domn->pram->pres,domn->dvisc->d[i],M,k);
                    Mcnd[k] *= DIMER*m_dimer*k;
                }
            }

            //---------- growth terms

            vector<double> Mgrw(nsvar,0.0);                         // growth source terms for moments

            double Acoef = M_PI*pow(abs(6.0/M_PI/rhoSoot),2.0/3.0); // Acoef = kmol^2/3 / kg^2/3
            for (int k=1; k<N; k++)                                 // Mgrw[0] = 0.0 by definition
                Mgrw[k] = Kgrw * Acoef * k * MOMIC(k-1.0/3.0,M);    // kg^k/m3*s

            Mgrw[0] = 0.0;                                          // by definition

            //---------- oxidation terms

            vector<double> Moxi(nsvar,0.0);

            for (int k=1; k<N; k++)                                 // Moxi[0] = 0.0 by definition
                Moxi[k] = Koxi * Acoef * k * MOMIC(k-1.0/3.0,M);    // kg^k/m3*s

            Moxi[0] = 0.0;                                          // by definition

            //---------- coagulation terms

            vector<double> Mcoa(nsvar,0.0);

            if (domn->pram->coagulation_mech != "NONE") {
                for (int k=0; k<N; k++) {
                    Mcoa[k] = getCoag(domn->temp->d[i],domn->pram->pres,domn->dvisc->d[i],M,k);
                }
            }

            Mcoa[1] = 0.0;                                         // by definition; redundant

            //---------- combinine to make source terms

            for (int k=0; k<nsvar; k++)
                rrSvar[k][i] = (Mnuc[k] + Mcnd[k] + Mgrw[k] - Moxi[k] + Mcoa[k]) / domn->rho->d[i];

            //---------- compute gas source terms

            set_gasSootSources(Mnuc[1], Mcnd[1], Mgrw[1], Moxi[1], i);

        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);                                   // assign to domn's rhsSrc variable

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);

    L_source_done = true;
}

////////////////////////////////////////////////////////////////////////////////
/*! lagrangeInterp function
 *
 *      Calculates the Lagrange interpolated value from whole order moments.
 *
 *      @param x_i  \input      x value of desired interpolation
 *      @param x    \input      vector of x values to interpolate amongst
 *      @param y    \input      vector of y values to interpolate amongst
 *      @param y_i  \output     interpolated y value
 *
 */

double dv_soot_MOMIC::lagrangeInterp(double x_i, vector<double> x, vector<double> y) {

    double y_i = 0.0;

    for (int j = 0; j < x.size(); j++) {
        double L = 1.0;
        for (int m = 0; m < x.size(); m++) {
            if (m != j) {
                L = L * (x_i - x[m])/(x[j] - x[m]);
            }
        }
        y_i = y_i + y[j] * L;
    }

    return y_i;

}

////////////////////////////////////////////////////////////////////////////////
/*! MOMIC function
 *
 *      Calculates the desired fractional moment by lagrange interpolation
 *      between whole order moments. Because it uses log moments, it will crash
 *      if any moment is less than or equal to zero.
 *
 *      @param p     \input     desired interpolation value
 *      @param M     \input     vector of whole order moments
 *
 */

double dv_soot_MOMIC::MOMIC(double p, vector<double> M) {

    if (p == 0) {
        return M[0];
    }

    int size = M.size();
    if (p < 0 && M.size() != 2) {
        size = 3;
    }

    vector<double> log_mu(size,0.0);
    vector<double> x(size,0.0);

    for (int i = 0; i < size; i++) {
        log_mu[i] = log10(M[i] / M[0]);     // reduced moments
        x[i] = i;
    }

    double log_mu_p = lagrangeInterp(p,x,log_mu);

    return pow(10.0, log_mu_p) * M[0];

}

////////////////////////////////////////////////////////////////////////////////
/*! f_grid function
 *
 *      Calculates the grid function described in Frenklach 2002 MOMIC paper
 *      using lagrange interpolation between whole order moments
 *
 *      @param x     \input x grid point
 *      @param y     \input y grid point
 *      @param M     \input vector of whole order moments
 *
 */

double dv_soot_MOMIC::f_grid(int x, int y, vector<double> M) {

    double f1_0 = MOMIC(x-1.0/2.0,M)*MOMIC(y+1.0/6.0,M) + 2.0*MOMIC(x-1.0/6.0,M)*MOMIC(y-1.0/6.0,M) + MOMIC(x+1.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    double f1_1 = MOMIC(x-1.0/2.0,M)*MOMIC(y+7.0/6.0,M) + 2.0*MOMIC(x-1.0/6.0,M)*MOMIC(y+5.0/6.0,M) + MOMIC(x+1.0/6.0,M)*MOMIC(y+1.0/2.0,M) +
                  MOMIC(x+1.0/2.0,M)*MOMIC(y+1.0/6.0,M) + 2.0*MOMIC(x+5.0/6.0,M)*MOMIC(y-1.0/6.0,M) + MOMIC(x+7.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    if (y >= 4) {

        vector<double> temp_x(2,0.0);
        temp_x[0] = 0.0;
        temp_x[1] = 1.0;

        vector<double> temp_y(2,0.0);
        temp_y[0] = log10(f1_0);
        temp_y[1] = log10(f1_1);

        double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

        return pow(10.0, value);
    }

    double f1_2 =     MOMIC(x-1.0/2.0,M)*MOMIC(y+13.0/6.0,M) + 2.0*MOMIC(x-1.0 /6.0,M)*MOMIC(y+11.0/6.0,M) +     MOMIC(x+1.0 /6.0,M)*MOMIC(y+3.0/2.0,M) +
                  2.0*MOMIC(x+1.0/2.0,M)*MOMIC(y+7.0 /6.0,M) + 4.0*MOMIC(x+5.0 /6.0,M)*MOMIC(y+5.0 /6.0,M) + 2.0*MOMIC(x+7.0 /6.0,M)*MOMIC(y+1.0/2.0,M) +
                      MOMIC(x+3.0/2.0,M)*MOMIC(y+1.0 /6.0,M) + 2.0*MOMIC(x+11.0/6.0,M)*MOMIC(y-1.0 /6.0,M) +     MOMIC(x+13.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    if (y >= 3) {

        vector<double> temp_x(3,0.0);
        temp_x[0] = 0.0;
        temp_x[1] = 1.0;
        temp_x[2] = 2.0;

        vector<double> temp_y(3,0.0);
        temp_y[0] = log10(f1_0);
        temp_y[1] = log10(f1_1);
        temp_y[2] = log10(f1_2);

        double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

        return pow(10.0, value);
    }

    double f1_3 =     MOMIC(x-1.0/2.0,M)*MOMIC(y+19.0/6.0,M) + 2.0*MOMIC(x-1.0 /6.0,M)*MOMIC(y+17.0/6.0,M) +     MOMIC(x+1.0 /6.0,M)*MOMIC(y+5.0/2.0,M) +
                  3.0*MOMIC(x+1.0/2.0,M)*MOMIC(y+13.0/6.0,M) + 6.0*MOMIC(x+5.0 /6.0,M)*MOMIC(y+11.0/6.0,M) + 3.0*MOMIC(x+7.0 /6.0,M)*MOMIC(y+3.0/2.0,M) +
                  3.0*MOMIC(x+3.0/2.0,M)*MOMIC(y+7.0 /6.0,M) + 6.0*MOMIC(x+11.0/6.0,M)*MOMIC(y+5.0 /6.0,M) + 3.0*MOMIC(x+13.0/6.0,M)*MOMIC(y+1.0/2.0,M) +
                      MOMIC(x+5.0/2.0,M)*MOMIC(y+1.0 /6.0,M) + 2.0*MOMIC(x+17.0/6.0,M)*MOMIC(y-1.0 /6.0,M) +     MOMIC(x+19.0/6.0,M)*MOMIC(y-1.0/2.0,M);

    vector<double> temp_x(4,0.0);
    temp_x[0] = 0.0;
    temp_x[1] = 1.0;
    temp_x[2] = 2.0;
    temp_x[3] = 3.0;

    vector<double> temp_y(4,0.0);
    temp_y[0] = log10(f1_0);
    temp_y[1] = log10(f1_1);
    temp_y[2] = log10(f1_2);
    temp_y[3] = log10(f1_3);

    double value = lagrangeInterp(1.0/2.0, temp_x, temp_y);

    return pow(10.0, value);

}

////////////////////////////////////////////////////////////////////////////////
/*! getCoag function
 *
 *      Calculates coagulation rate for MOMIC based on a weighted average of
 *      continuum and free-molecular values. See Frenklach's 2002 MOMIC paper.
 *      Adapted from python code by Alex Josephson.
 *
 *      @param T    \input  local gas temperature (K)
 *      @param P    \input  local gas pressure (Pa)
 *      @param mu   \input  local dynamic/absolute gas viscosity (kg/s/m)
 *      @param M    \input  vector of previous moment values
 *      @param r    \input  number of the moment to be calculated
 *
 */

double dv_soot_MOMIC::getCoag(double T, double P, double mu, vector<double> M, int r) {

    if (r == 1) {                                               // coagulation does not affect M1
        return 0;                                               // this is a shortcut to save a little computation
    }
    double kb    = 1.38064852E-23;

    // Calculate Knudsen number to determine regime

    double mu_1     = M[1]/M[0];                                // average particle mass (kg)
    double d_g      = pow(6.0*kb*T/P/M_PI, 1.0/3.0);            // average gas molecular diameter (m)
    double d_p      = pow(6.0*mu_1/rhoSoot/M_PI, 1.0/3.0);     // average particle diameter (m)
    double lambda_g = kb*T/(pow(2.0,0.5)*M_PI*pow(d_g,2.0)*P);  // gas mean free path (m)
    double Kn       = lambda_g/d_p;                             // Knudsen number

    // Continuum regime

    double Rate_C = 0.0;

    double K_C = 2.0*kb*T/(3.0*mu);
    double K_Cprime = 1.257*lambda_g*pow(M_PI*rhoSoot/6.0,1.0/3.0);

    if (r == 0) {
        Rate_C = -K_C*(pow(M[0],2.0) + MOMIC(1.0/3.0,M)*MOMIC(-1.0/3.0,M) +
                  K_Cprime*(3.0*MOMIC(-1.0/3.0,M)*M[0] + MOMIC(2.0/3.0,M)*MOMIC(1.0/3.0,M)));
    }
    else {
        Rate_C = 0.0;
        for (int k=1; k<r; k++) {
            if (k > r-k) {
                Rate_C = Rate_C;
            }
            else {
                Rate_C = Rate_C + boost::math::binomial_coefficient<double>(r,k)*
                         (2.0*M[k]*M[r-k] + MOMIC(k+1.0/3.0,M)*MOMIC(r-k-1.0/3.0,M) + MOMIC(k-1.0/3.0,M)*MOMIC(r-k+1.0/3.0,M) +
                          2.0*K_Cprime* (2.0*MOMIC(k-1.0/3.0,M)*M[r-k] + M[k]*MOMIC(r-k-1.0/3.0,M) + MOMIC(k-2.0/3.0,M)*MOMIC(r-k+1.0/3.0,M)));
            }
        }
        Rate_C = 0.5*K_C*Rate_C;
    }

    // Free-molecular regime

    double Rate_F = 0.0;

    double K_f = 2.2*pow(3.0/(4.0*M_PI*rhoSoot),2.0/3.0) * pow(8.0*M_PI*kb*T,1.0/2.0);

    if (r == 0) {
        Rate_F = -0.5*K_f*f_grid(0,0,M);
    }
    else {
        Rate_F = 0.0;
        for (int k=1; k<r; k++) {
            if (k > r-k) {
                Rate_F = Rate_F;
            }
            else {
                Rate_F = Rate_F + boost::math::binomial_coefficient<double>(r,k)*f_grid(k,r-k,M);
            }
        }
        Rate_F = 0.5*K_f*Rate_F;
    }

    return Rate_F/(1+1/Kn) + Rate_C/(1+Kn);

}

////////////////////////////////////////////////////////////////////////////////
/*! downselectIfNeeded function
 *
 *      Checks the moment set can be used with MOMIC. Downselects if needed.
 *      Zero and negative moment values will trip the flag for downselection.
 *      The floor value is N = 2.
 *
 *      For the case where M0 <= 0, we set N = 0, which means that the source
 *      terms will be initalized but not calculated such that the rates are all
 *      zero. If M1 <= 0, we assign it a value based on M0 and a lognormal
 *      distribution that matches the initialized profile in the domc. This
 *      is a workaround since the lagrange interpolation can't handle M1 = 0.
 *
 *      @param ipt  \input  grid point to evaluate at
 *      @param &M   \input  moment set
 *      @param &N   \input  number of downselected moments
 *
 */

void dv_soot_MOMIC::downselectIfNeeded(int ipt, vector<double> &M, int &N) {

    // CHECK: M0 <= 0.0

    if (M[0] <= 0.0) {
        N = 0;
        return;
    }

    // CHECK: M1 <= 0.0

    if (M[1] <= 0.0) {
        double M0 = 1.0;
        double sigL = 3.0;                              // sigL and mavg should be same as in domaincase
        double mavg = 1.0E-21;
        M[1] = M0 * mavg * exp(0.5 * pow(sigL,2.0));    // give M1 a value based on M0 and lognormal dist.
        domn->svar[1]->d[ipt] = M[1];
    }

    // CHECK: all remaining moments

    bool zeros = false;                                 // flag for breaking do/while loop

    do {

        if (N <= 2) {                                   // will downselect to 2 moments, but not further
            break;
        }

        zeros = false;                                  // reset flag

        for (int i=0; i<N; i++) {
           if (M[i] <= 0.0) {
               zeros = true;
           }
        }

        if (zeros == true) {                            // if flagged, downselect by two moments
            N = N - 1;
        }

    } while (zeros != 0);                               // keep going until we get no flag

    M.resize(N);                                        // resize M based on downselected N

}
