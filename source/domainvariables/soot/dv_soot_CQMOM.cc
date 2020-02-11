/**
 * @file dv_soot_CQMOM.cc
 * Header file for class dv_soot_CQMOM
 * @author Victoria B. Lansinger
 */

#include "dv_soot_CQMOM.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>
#include <cfloat>      // DBL_EPSILON
#include <vector>
#include <algorithm>

#include "eispack.h"

void adaptive_CQMOM(int N_v, int N_s, const vector<vector<double> > &m,
                    const double &eabs_v, const double &eabs_s,
                    const vector<double> &rmin_v, const vector<double> &rmin_s,
                    int &Nuse_v, vector<int> &Nuse_s, vector<vector<double> > &W,
                    vector<vector<double> > &X_v, vector<vector<double> > &X_s );

//extern "C" void dstev_(char *JOBZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 *
 *  copied from dv_ygas.cc
 */

void dv_soot_CQMOM::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);                     // resize rhsSrc

    static vector<vector<double> > rrSvar(nsvar);        // temp storage for moment rates [nsvar][ngrd]
    int nsvar_v = domn->pram->nsvar_v;
    int nsvar_s = domn->pram->nsvar_s;

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    for(int k=0; k<nsvar; k++)                                        // for each mom/section,
        rrSvar.at(k).resize(domn->ngrd);                            // resize to match grid size

    if(kMe==0) {                                                    // to save cost, compute needed terms for all dv_soot_CQMOM objects using this one.

        for(int k=0; k<nsvar; k++)                                    // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);                        // resize to match grid size

        for(int i=iS; i<=iE; i++) {                                    // loop over grid points

            domn->domc->setGasStateAtPt(i);                            // set gas state

            // prepare moment values for CQMOM solver
            vector<vector<double> > mom_matrix(nsvar_v);       // define matrix storage
            for (int k=0; k<nsvar_v; k++) {
                mom_matrix[k].resize(nsvar_s, 0.0);                 // resize and initialize with zeros
            }

            int counter = 0;                                                            // set counter for number of transported moments
            for (int k=0; k<nsvar_v; k++) {
                mom_matrix[k][0] = domn->svar[k]->d[i] * domn->rho->d[i];               // fills first s column (longer than the rest)
                counter++;                                                              // counter increase
            }
            for (int j=1; j<nsvar_s; j++) {
                for (int k=0; k<nsvar_v/2; k++) {
                    mom_matrix[k][j] = domn->svar[counter]->d[i] * domn->rho->d[i];     // fills the rest of the s columns
                }
            }

            // set variables needed for adaptive wheeler
            double eabs_v = domn->io->sootParams["eabs_v"].as<double>();
            double eabs_s = domn->io->sootParams["eabs_s"].as<double>();

            vector<double> rmin_v(nsvar_v/2, 1e-5);
            vector<double> rmin_s(nsvar_s/2, 1e-5);

            nused_v = 0;
            nused_s.resize(nsvar_s/2,0);

            // do adaptive CQMOM calculation
            adaptive_CQMOM(nsvar_v/2, nsvar_s/2, mom_matrix, eabs_v, eabs_s, rmin_v, rmin_s,
                            nused_v, nused_s, wts, absc_v, absc_s);

            // get chemical soot rates
            double Jnuc = chem->getNucleationRate(i);                // #/m3*s
            double Kgrw = chem->getGrowthRate(i);                    // kg/m2*s
            double Koxi = chem->getOxidationRate(i);                // kg/m2*s
            //double Kfc  = chem->getCoagulationRate(i);  //              // #/m3*s

            // nucleation terms
            vector<vector<double> > Mnuc(nsvar_v);
            for (int k=0; k<nsvar_v; k++) {
               Mnuc[k].resize(nsvar_s,0.0);
            }

            double m_min = domn->pram->Cmin*chem->MWc/chem->Na;
            double s_min = pow(6*pow(M_PI,2.0)*chem->MWc*domn->pram->Cmin/domn->pram->rho_soot/chem->Na,2.0/3.0);

            for (int k=0; k<nused_v*2; k++) {
                for (int j=0; j<nused_s[k]*2; j++) {
                    Mnuc[k][j] = pow(m_min,k) * pow(s_min,j) * Jnuc;
                }
            }

            // growth terms
            vector<vector<double> > Mgrw(nsvar_v);
            for (int k=0; k<nsvar_v; k++) {
               Mgrw[k].resize(nsvar_s,0.0);
            }

            for (int k=0; k<nused_v*2; k++) {
                for (int j=0; j<nused_s[k]*2; j++) {
                    Mgrw[k][j]    = -(k + 2.0/3.0*j) * Mkj(k-1,j+1) * Kgrw;            // Gk = r* M_(k-1/3)*Ap*ks = r*M(k-1/3)*Rgrw*MWc
                }
            }

            // oxidation terms
            vector<vector<double> > Moxi(nsvar_v);
            for (int k=0; k<nsvar_v; k++) {
               Moxi[k].resize(nsvar_s,0.0);
            }

            for (int k=0; k<nused_v*2; k++) {
                for (int j=0; j<nused_s[k]*2; j++) {
                    Moxi[k][j]    = (k + 2.0/3.0*j) * Mkj(k-1,j+1) * Koxi;            // Xk = r* M_(k-1/3)*Ap*ks = r*M(k-1/3)*Rgrw*MWc
                }
            }

            // coagulation terms
            vector<vector<double> > Mcoa(nsvar_v);
            for (int k=0; k<nsvar_v; k++) {
               Mcoa[k].resize(nsvar_s,0.0);
            }

            // getCoag(Mcoa,i);

            // combinine to make source terms
            vector<vector<double> > rrSvar_temp(nsvar_v);       // size [nsvar_v][nsvar_s]
            for (int k=0; k<nsvar_v; k++) {
               rrSvar_temp[k].resize(nsvar_s,0.0);
            }

            for (int k=0; k<nused_v*2; k++) {
                for (int j=0; j<nused_s[k]*2; j++) {
                    rrSvar_temp[k][j] = (Mnuc[k][j] + Mgrw[k][j] - Moxi[k][j] + Mcoa[k][j]) / domn->rho->d[i];
                }
            }

            // convert back to vector form for rrSvar
            int counter2 = 0;                                                            // set counter for number of transported moments
            for (int k=0; k<nsvar_v; k++) {
                rrSvar[k][i] = rrSvar_temp[k][0];
                counter2++;                                                              // counter increase
            }
            for (int j=1; j<nsvar_s; j++) {
                for (int k=0; k<nsvar_v/2; k++) {
                    rrSvar[counter2][i] = rrSvar_temp[k][j];
                }
            }


        }
    }

     for(int i=iS; i<=iE; i++)
         rhsSrc.at(i) = rrSvar.at(kMe).at(i);

     if(domn->pram->Lspatial)
         for(int i=iS; i<=iE; i++)
             rhsSrc.at(i) /= domn->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! Mk function
 *    Calculates fractional moments from weights and abscissas
 *  CQMOM moment calculation
 *
 *    @param wts     \input     weights; combined conditional*regular
 *    @param absc    \input     abscissas
 *    @param Mk    \output fractional moment value
 *
 */

double dv_soot_CQMOM::Mkj(double v, double s) {

    double Mkj = 0;

    for (int k=0; k<nused_v; k++) {
        for(int j=0; k<nused_s[k]; j++) {
            if (wts[k][j] == 0 || absc_v[k][j] == 0 || absc_s[k][j] == 0)
                Mkj += 0;
            else
                Mkj += wts[k][j] * pow(absc_v[k][j],v) * pow(absc_s[k][j],s);
        }
    }

    return Mkj;

}

