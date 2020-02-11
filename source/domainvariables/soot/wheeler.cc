#include<vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "eispack.h"

using namespace std;

//extern "C" void dstev_(char *JOBZ, int *N, double *D, double *E, double *Z, int *LDZ, double *WORK, int *INFO);

/**--------------------------------------------------------------------------------
 * David Lignell
 * Adaptive Wheeler algorithm for computing weights and abscissas from moments.
 * From Yuan and Fox (2011) Journal of Computational Physics 230(10) P. 8216-8246.
 * LApack's dstev function to compute eigenvalues and eigenvectors of symm tridiag matrix.

 * @param m    \input vector of moments (size = 2N)
 * @param N    \input number of quadrature nodes (abscissas)
 * @param rmin \input vector of ratios wmin/wmax for different number of nodes used.
                      rmin[0]    = 'vacuum' state.
                      rmin[nl-1] = minimum ratio for nl nodes: nl=2, 3, ...N
 * @param eabs \input ratio min(abs dx)/max(abs dx), where dx is the spacing between any two abscissas.
                      That is, the minimium relative distance between distinct absicass.
                      eabs = 1E-3 is quoted by Yuan, but that would seem to require too high a spread in the abscissas for soot problems (?)
 * @param Nout \output number of nodes used
 * @param w    \output vector weights
 * @param x    \output vector abscissas
 *
 * From the paper (note: converted matlab indexing to C++):
 *     The input parameters rmin and eabs are problem dependent.
 *     Generally, the optimal values for rmin depend on how accurately
 *     the moments mom are reproduced by the finite-volume spatial fluxes
 *     (see Appendix C for details), with less accu- rate moments requiring larger values.
 *     (Note that setting rmin[N-1]=1 will eliminate the possibility of having N nodes.)
 *     In most problems, rmin[0] (which controls the ‘vacuum’ state) can be set very small.
 *     The parameter rmin[1], which controls the switch from one-to two-node quadrature,
 *     can be relatively small since lower-order moments are typically more accurate.
 *     In contrast, rmin[N-1] >= 0 for N >= 3 should be set more conservatively.
 *     The parameter eabs controls the condition number of the Vandermonde matrices V1 and V2_a1 .
 *     Values for eabs around 1E-3 usually provide acceptable accuracy for multi-variable quadrature.

 */

void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs,
                      int &Nout, vector<double> &w, vector<double> &x ) {

    int                     N2 = N*2;
    vector<vector<double> > sig(N+1,vector<double>(N2+1, 0.0));
    vector<double>          a(N,0.0);
    vector<double>          b(N,0.0);

    vector<double>          eval(N);
    vector<double>          evec(N*N);
    vector<double>          j_diag(N);
    vector<double>          j_ldiag(N);

    double cutoff = 0.0;
    //int    werror = 0;

    w = vector<double>(N, 0.0);
    x = vector<double>(N, 0.0);

    if(m[0] <= 0)
        return;       // with w = x = 0

    if(N==1 || m[0] < rmin[0]) {
        Nout = 1;
        w[0] = m[0];
        x[0] = m[1]/m[0];
        return;
    }

    //---------- construct recurrence matrix

    for(int i=1; i<N2+1; i++)
        sig[1][i] = m[i-1];
    a[0] = m[1]/m[0];
    for(int k=2;k<N+1; k++) {
        for(int l=k; l<N2-k+2; l++)
            sig[k][l] = sig[k-1][l+1]-a[k-2]*sig[k-1][l]-b[k-2]*sig[k-2][l];
        a[k-1] = sig[k][k+1]/sig[k][k]-sig[k-1][k]/sig[k-1][k-1];
        b[k-1] = sig[k][k]/sig[k-1][k-1];
    }

    //---------- get max N using diag elements of sig
    int Norig = N;
    for(int k=N; k >=2; k--)
        if(sig[k][k] <= cutoff) {
            N = k-2;
            if(N==1) {
                Nout = 1;
                w[0] = m[0];
                x[0] = m[1]/m[0];
                return;
            }
        }

    //---------- compute quadrature using maximum n

    N2  = N*2;
    a   = vector<double>(N, 0.0);
    b   = vector<double>(N, 0.0);
    sig = vector<vector<double> >(N+1,vector<double>(N2+1, 0.0));

    for(int i=1; i<N2+1; i++)
        sig[1][i] = m[i-1];
    a[0] = m[1]/m[0];
    for(int k=2;k<N+1; k++) {
        for(int l=k; l<N2-k+2; l++)
            sig[k][l] = sig[k-1][l+1]-a[k-2]*sig[k-1][l]-b[k-2]*sig[k-2][l];
        a[k-1] = sig[k][k+1]/sig[k][k]-sig[k-1][k]/sig[k-1][k-1];
        b[k-1] = sig[k][k]/sig[k-1][k-1];
    }

    //---------- check if moments are not realizable (should never happen)

    double bmin = *min_element(b.begin(), b.end());
    if(bmin < 0) {
        //werror = 1;
        cout << endl << "Moments in adaptive_wheeler are not realizable!" << endl;
        return;
    }

    //---------- setup Jacobi matrix for N-point quadrature, adapt N using rmax and eabs

    for(int Nl=N; Nl>=1; Nl--) {

        if(Nl==1) {
            Nout = 1;
            w[0] = m[0];
            x[0] = m[1]/m[0];
            return;
        }

        vector<double> j_diag(Nl);
        vector<double> j_ldiag(Nl);
        vector<double> eval(Nl);
        vector<double> evec(Nl*Nl);
        x = vector<double>(Norig, 0.0);
        w = vector<double>(Norig, 0.0);

        for(int i=1; i<=Nl-1; i++) {
            j_diag[i-1] = a[i-1];
            j_ldiag[i] = -sqrt(abs(b[i]));
        }
        j_diag[Nl-1] = a[Nl-1];

        for(int i=0; i<Nl; i++)
            evec[i+Nl*i] = 1.0;

        int flag = tql2(Nl, &j_diag[0], &j_ldiag[0], &evec[0]);       // for eispack

        //char VorN = 'V';
        //vector<double> work(2*Nl-2);
        //int info;
        //dstev_( &VorN, &Nl, &j_diag[0], &j_ldiag[1], &evec[0], &Nl, &work[0], &info);

        for(int j=0; j<Nl; j++) {
            x[j] = j_diag[j];               // j_diag are now the vector of eigenvalues.
            w[j] = m[0]*pow(evec[0+j*Nl],2.0);
        }

        vector<double> dab(Nl,0.0);
        vector<double> mab(Nl,0.0);
        double dmb;
        for(int i=Nl-1; i>=1; i--) {
            dab[i] = abs(x[i] - x[0]);
            mab[i] = dab[i];
            for(int k=1; k<i; k++) {
                dmb = abs(x[i] - x[k]);
                if(dmb < dab[i])
                    dab[i] = dmb;
                if(dmb > mab[i])
                    mab[i] = dmb;
            }
        }
        double mindab = *min_element(dab.begin()+1, dab.end());  // min abs of spacing between any two abscissas.
        double maxmab = *max_element(mab.begin()+1, mab.end());  // max abs of spacing between any two abscissas.
        if(Nl==2)
            maxmab = 1;

        double wmin = *min_element(w.begin(), w.begin()+Nl);
        double wmax = *max_element(w.begin(), w.begin()+Nl);
        if( wmin/wmax > rmin[Nl-1] && mindab/maxmab > eabs ) {   // make eabs smaller to be less stringent and allow closer abscissas
            Nout = Nl;                                           // make all rmin very small to be less stringent and allow small relative weights.
            return;
        }
    }
}

//#include"eispack.hpp"      // for eispack
//#include <string>          // for eispack
/**--------------------------------------------------------------------------------
 * David Lignell
 * Wheeler algorithm for computing weights and abscissas from moments.
 * From Marchisio and Fox (2013) Computational Models for Polydisperse and Multiphase systems.
 * // Uses eispack function tql2 for eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
 * // If eispack version tql2 is desired, download eispack.hpp and eispack.cpp.
 * LApack's dstev function to compute eigenvalues and eigenvectors of symm tridiag matrix.

 * @param m vector of moments (size = 2N)
 * @param N input number of quadrature nodes (abscissas)
 * @param w output weights
 * @param x output abscissas
 */

void wheeler(const vector<double> &m, int N, vector<double> &w, vector<double> &x ) {

    int                     N2 = N*2;
    vector<vector<double> > sigma(N+1,vector<double>(N2, 0.0));
    vector<double>          a(N,0.0);
    vector<double>          b(N,0.0);
    vector<double>          eval(N);
    vector<double>          evec(N*N);
    vector<double>          j_diag(N);
    vector<double>          j_ldiag(N);

    for(int l=0; l<=N2-1; l++)
        sigma[1][l] = m[l];


    a[0] = m[1]/m[0];

    for(int k=1; k<=N-1; k++) {
        for(int l=k; l<=N2-k-1; l++)
            sigma[k+1][l] = sigma[k][l+1]-a[k-1]*sigma[k][l]-b[k-1]*sigma[k-1][l];
        a[k] = -sigma[k][k]/sigma[k][k-1]+ sigma[k+1][k+1]/sigma[k+1][k];
        b[k] = sigma[k+1][k]/sigma[k][k-1];
    }


    j_diag = a;
    for(int i=1; i<=N-1; i++)
        j_ldiag[i] = -sqrt(abs(b[i]));

    for(int i=0; i<N; i++)
        evec[i+N*i] = 1.0;

    int flag = tql2(N, &j_diag[0], &j_ldiag[0], &evec[0]);       // for eispack

    //char VorN = 'V';
    //vector<double> work(2*N-2);
    //int info;
    //dstev_( &VorN, &N, &j_diag[0], &j_ldiag[1], &evec[0], &N, &work[0], &info);

    x = j_diag;      // j_diag are now the vector of eigenvalues.

    for(int j=0; j<N; j++)
        w[j] = pow(evec[0+j*N],2.0)*m[0];

}
