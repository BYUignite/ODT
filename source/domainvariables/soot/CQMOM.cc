#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <cfloat>      // DBL_EPSILON

using namespace std;

//-------------------------------------------------------------------------------------------

void vandermonde_solver(const vector<double> &x, vector<double> &w, const vector<double> &q, const int &n);
void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs,
                      int &Nout, vector<double> &w, vector<double> &x );

//-------------------------------------------------------------------------------------------

/*!
 * From Yuan and Fox (2011) Journal of Computational Physics 230(10) P. 8216-8246.
 * Two internal coordinates: v, s.
 *
 * @param N_v    \input  max number of nodes (environments, abscissas) in the v direction.
 * @param N_s    \input  max number of nodes (environments, abscissas) in the s direction.
 * @param m      \input  matrix of moments m[2][3] is m_2,3 is m_iv,is
 * @param eabs_v \input  for adaptive Wheeler algorithm for v direction
 * @param eabs_s \input  for adaptive Wheeler algorithm for s direction
 * @param rmin_v \input  vector for adaptive Wheeler algorithm for v direction
 * @param rmin_s \input  vector for adaptive Wheeler algorithm for s direction
 * @param Nuse_v \output number of nodes USED/TO USE in v direction
 * @param Nuse_s \output vector number of nodes USED/TO USE in s direction at each v node used.
 * @param W      \output matrix of weights at each node
 * @param X_v    \output matrix of v-locations of all nodes (this could be done with a vector)
 * @param X_s    \output matrix of s-locations of all nodes
 */

void adaptive_CQMOM(int N_v, int N_s, const vector<vector<double> > &m,
                    const double &eabs_v, const double &eabs_s,
                    const vector<double> &rmin_v, const vector<double> &rmin_s,
                    int &Nuse_v, vector<int> &Nuse_s, vector<vector<double> > &W,
                    vector<vector<double> > &X_v, vector<vector<double> > &X_s ) {

    //-----------------

    X_v    = vector<vector<double> >(N_v, vector<double>(N_s, 0.0));
    X_s    = vector<vector<double> >(N_v, vector<double>(N_s, 0.0));
    W      = vector<vector<double> >(N_v, vector<double>(N_s, 0.0));
    Nuse_v = N_v;
    Nuse_s = vector<int>(N_v, 0);

    double small = 10*DBL_EPSILON;

    //-----------------

    if(N_v==0 || N_s==0){
        cout << endl << "ERROR in adaptive_CQMOM, nodex and nodey must be > 0" << endl;
        exit(0);
    }
    //-----------------

    if(m[0][0] <= 0)
        return;
    else if(m[0][0] < rmin_v[0] ) {
        W[0][0] = m[0][0];
        X_v[0][0] = m[1][0]/m[0][0];
        X_s[0][0] = m[0][1]/m[0][0];
        return;
    }
    else if(m[0][0] < small*100) {
        N_v = (N_v < 2) ? N_v : 2;          //todo check N_v, N_s changes...
        N_s = (N_s < 2) ? N_s : 2;
    }

    //----------------- 1D quadrature in the v direction

    vector<double> x_v(N_v);
    vector<double> w_v(N_v);
    vector<double> dmb(N_v);

    for(int i=0; i<N_v; i++)
        dmb[i] = m[i][0];

    adaptive_wheeler(dmb, N_v, rmin_v, eabs_v, Nuse_v, w_v, x_v);

    //----------------- calculate the conditional moments

    vector<vector<double> > m_c(Nuse_v, vector<double> (N_s, 0.0));
    for(int i=0; i<Nuse_v; i++)                // initialize
        m_c[i][0] = 1.0;

    vector<vector<double> > V(Nuse_v, vector<double>(Nuse_v));
    for(int i=0; i<Nuse_v; i++)
        for(int j=0; j<Nuse_v; j++)
            V[i][j] = pow(x_v[j], i);         // Vandermonde matrix, used directly for iterative improvement of linear solve

    vector<double> b_dx(Nuse_v);              // rhs of linear solve, and dx for iterative improvement
    vector<double> sol(Nuse_v);

    for(int j=1; j<2*N_s; j++) {              // solve [V][R](mc_ij) = (m_ij) system for j=1 to 2*N_s-1
                                              // [V], [R] are Nuse_v x Nuse_v; (mc_ij), (m_ij) are Nuse_v
        for(int i=0; i<Nuse_v; i++)           // [R] is diagonal with w_v elements
            b_dx[i] = m[i][j];

        vandermonde_solver(x_v, sol, b_dx, Nuse_v);      // solving [V]([R](mu)) = m for ([R](mu)) = sol
                                                         // sol is the unknown, b_dx is the RHS, x_v vector defines the Vandermond "A" matrix
        //--------- iterative improvement (1 iter), optional
        {
            //--- 1: solve A(x+dx) = b for (x+dx). This was just done above, and (x+dx) = sol.
            //--- 2: form the error: ( A*(x+dx) - b ) = (err)
            vector<double> err(Nuse_v, 0.0);             // form err
            for(int ii=0; ii<Nuse_v; ii++) {
                for(int jj=0; jj<Nuse_v; jj++)
                    err[ii] += V[ii][jj]*sol[ii];
                err[ii] -= b_dx[ii];
            }
            //--- 3: solve A(dx) = (err) for (dx)
            vandermonde_solver(x_v, b_dx, err, Nuse_v);    // solve for dx. b = dx here
            //--- 4: solve for x = (x+dx) - dx.
            for(int i=0; i<Nuse_v; i++)
                sol[i] -= b_dx[i];
        }

        //-------- Convert from ([R](mu)) to (mu) by "/w_v[i]".

        for(int i=0; i<Nuse_v; i++)
            m_c[i][j] = sol[i] / w_v[i];
    }

    //----------------- calculate the conditional weights and abscissas

    for(int i=0; i<Nuse_v; i++)
        adaptive_wheeler(m_c[i], N_s, rmin_s, eabs_s, Nuse_s[i], W[i], X_s[i]);

    //----------------- finalize: reset W and X_v

    for(int i=0; i<Nuse_v; i++) {
        for(int j=0; j<Nuse_s[i]; j++) {
            X_v[i][j] = x_v[i];
            W[i][j]   *= w_v[i];
        }
    }

}

////-------------------------------------------------------------------------------------------

//void wheeler(const vector<double> &m, int N, vector<double> &w, vector<double> &x );

//extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info);
//extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

////-------------------------------------------------------------------------------------------
//
///*!
// * From Marchisio and Fox (2013) Computational Models for Polydisperse and Multiphase systems.
// * Two internal coordinates: v, s.
// */
//
//void CQMOM(int N_v, int N_s, vector<vector<double> > &m,
//           vector<double> &W, vector<vector<double> > &X ) {
//
//    double cutoff = 1.0E-6;
//    int    flag = 0;
//
//    //--------- construct the quadrature in the first direction
//
//    vector<double> x_v(N_v);
//    vector<double> w_v(N_v);
//    vector<double> dmb(N_v);
//
//    for(int i=0; i<N_v; i++)
//        dmb[i] = m[i][0];
//
//    wheeler(dmb, N_v, w_v, x_v);
//
//    //--------- check it out...
//
//    for(int i=0; i<N_v; i++)
//        if(abs(w_v[i])/m[0][0] < cutoff){
//            cout << endl << "One of the weights is null! Reduce the number of nodes in direction 0" << endl;
//            flag = 1;
//        }
//
//    //---------
//
//    if(flag == 0) {
//
//        //-------- define the Vandermonde matrix (* the R matrix) for conditional moments
//
//        vector<double> VR(N_v*N_v);
//
//        for(int j=0; j<N_v; j++)
//            for(int i=0; i<N_v; i++)
//                VR[i+N_v*j] = pow(x_v[j], i) * w_v[j];    // column-major format (fortran style) for LApack solver
//
//        //--------- calculate the conditional moments
//
//        int  ipiv[N_v];
//        int  info;
//        char trans = 'N';
//        int  nrhs = 1;
//
//        vector<vector<double> > m_c(N_v, vector<double> (N_s));
//
//        for(int i=0; i<2*N_s; i++) {
//
//            for(int k=0; k<N_v; k++)
//                dmb[k] = m[k][i];
//
//            dgetrf_(&N_v, &N_v, &VR[0], &N_v, ipiv, &info);
//            dgetrs_(&trans, &N_v, &nrhs, &VR[0], &N_v, ipiv, &dmb[0], &N_v, &info);   // dmb is b --> x in Ax=b
//
//            for(int k=0; k<N_v; k++)       // store solution in m_c (conditional moment)
//                m_c[k][i] = dmb[k];
//        }
//
//        //--------- calculate the conditional weights and abscissas
//
//        vector<vector<double> > x_sc(N_s, vector<double>(N_v));
//        vector<vector<double> > w_sc(N_s, vector<double>(N_v));
//
//        dmb.resize(N_s);
//        vector<double> dmb1(N_s);
//
//        for(int i=0; i<N_v; i++) {
//
//            wheeler(m_c[i], N_s, dmb, dmb1);
//
//            for(int k=0; k<N_s; k++) {
//                w_sc[k][i] = dmb[k];
//                x_sc[k][i] = dmb1[k];
//            }
//        }
//
//        //--------- finish up
//
//        int indx = 0;
//        for(int i=0; i<N_v; i++) {
//            for(int j=0; i<N_s; j++){
//                W[indx]    = w_v[i]*w_sc[j][i];
//                X[0][indx] = x_v[i];
//                X[1][indx] = x_sc[j][i];
//                indx++;
//            }
//        }
//
//    } // end if(flag==0)
//
//    //---------
//
//    else {
//
//        int indx = 0;
//        for(int i=0; i<N_v; i++) {
//            for(int j=0; i<N_s; j++){
//                W[indx]   = 0.0;
//                X[0][indx] = 0.0;
//                X[1][indx] = 0.0;
//                indx++;
//            }
//        }
//
//    } // end else
//
//}
