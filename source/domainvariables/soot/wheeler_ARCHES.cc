#include<vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "eispack.h"

using namespace std;

// declare lapack eigenvalue solver
extern "C"{
  void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
              double* w, double* work, int* lwork, int* info );

  void dgesv(int *n, int *nrhs, double *a, int *lda,
             int *ipiv, double *b, int *ldb, int *info);
}

//uncomment to debug matricies
//#define cqmom_dbg

/***************************
 N-D Adaptive Wheeler alogrithm
 ****************************/
void adaptiveWheelerAlgorithm(const std::vector<double>& moments, std::vector<double>& w,
                              std::vector<double>& x,const double& rMin, const double& eAbs)
{
  using std::vector;
  using std::cout;
  using std::endl;

  int nEnvMax = moments.size()/2;
  int nEnvOut = moments.size()/2;

#ifdef cqmom_dbg
  int nMom = moments.size();
  cout << "Wheeler Adaptive Start" << endl;
  for (int i = 0; i < nMom; i++) {
    cout << "m[" << i << "]=" << moments[i] << endl;
  }
#endif
  bool isRealizable = false;

  while ( !isRealizable ) {
    if (nEnvOut == 1) {
      w[0] = moments[0];
      x[0] = moments[1]/moments[0];
      double d_small = 1.0e-10;
      if ( fabs(x[0]) < d_small) {
        //prevent very small values from propagating junk later
        x[0] = 0.0;
      }
#ifdef cqmom_dbg
      cout << "Singular Point" << endl;
#endif
      isRealizable = true;
      continue;
    }

    vector<double> a (nEnvOut, 0.0);
    vector<double> b (nEnvOut, 0.0);
    vector< vector<double> > sig (2*nEnvOut+1, vector<double> (2*nEnvOut+1,0.0) );

    for ( int i = 1; i<(2*nEnvOut+1); i++) {
      sig[1][i] = moments[i-1];
    }
    a[0] = moments[1]/moments[0];

    for (int k = 2; k<(nEnvOut+1); k++) {
      for (int j = k; j<(2*nEnvOut-k+2); j++) {
        sig[k][j] = sig[k-1][j+1]-a[k-2]*sig[k-1][j]-b[k-2]*sig[k-2][j];
      }
      a[k-1] = sig[k][k+1]/sig[k][k] - sig[k-1][k]/sig[k-1][k-1];
      b[k-1] = sig[k][k]/sig[k-1][k-1];
    }

#ifdef cqmom_dbg
    for (int i=0; i<nEnvOut; i++) {
      cout << "a[" << i << "] = " << a[i] <<endl;
    }
    for (int i=0; i<nEnvOut; i++) {
      cout << "b[" << i << "] = " << b[i] <<endl;
    }
    for (int i=0; i<2*nEnvOut+1; i++) {
      for (int j = 0; j<2*nEnvOut+1; j++) {
        cout << "S[" << i << "][" << j << "] = " << sig[i][j] << "  ";
      }
      cout << endl;
    }
#endif

    bool nonrealCheck = false;
    //check a vector for a nan - occurs in point distribution
    for ( int i = 0; i<nEnvOut; i++ ) {
      if ( std::isnan(a[i]) || std::isinf(a[i]) ) {
#ifdef cqmom_dbg
        cout << "WARNING: Arches: CQMOMInversion: not-a-number in a vector encountered. " << endl;
#endif
        nEnvOut--;
        nonrealCheck = true;
        break;
      }
    }

    if (nonrealCheck)
      continue;

    double d_small = 1.0e-14;
    //check the b vector for realizable space
    for ( int i = 0; i<nEnvOut; i++ ) {
      if ( (b[i] != 0.0 && b[i]<d_small) || std::isnan(b[i]) ) { //clip if b is very close to zero
#ifdef cqmom_dbg
        cout << "WARNING: Arches: CQMOMInversion: Negative b vector encountered." << endl;
#endif
        nEnvOut--;
        nonrealCheck = true;
        break;
      }
    }

    if (nonrealCheck)
      continue;

#ifdef cqmom_dbg
    vector< vector<double> > z (nEnvOut, vector<double> (nEnvOut, 0.0) );
    for ( int i = 0; i<(nEnvOut-1); i++) {
      z[i][i] = a[i];
      z[i][i+1] = sqrt( b[i+1] );
      z[i+1][i] = z[i][i+1];
    }
    z[nEnvOut-1][nEnvOut-1] = a[nEnvOut-1];
    cout << "made z matrix" << endl;
    for (int i=0; i<nEnvOut; i++) {
      for (int j = 0; j<nEnvOut; j++) {
        cout << "z[" << i << "][" << j << "] = " << z[i][j] << "  ";
      }
      cout << endl;
    }
#endif

    vector<double> z_(nEnvOut*nEnvOut,0.0);
    for (int i = 0; i<(nEnvOut-1); i++) {
      z_[i*nEnvOut + i] = a[i];
      z_[i*nEnvOut + i + 1] = sqrt( b[i+1] );
      z_[(i+1)*nEnvOut + i] = z_[i*nEnvOut + i + 1];
    }
    z_[nEnvOut*nEnvOut-1] = a[nEnvOut-1];

    vector<double> eigenVal (nEnvOut, 0.0);
    //_____________
    //solve eigenvalue problem with external package
    int  lda = nEnvOut, info, lwork;
    double wkopt;
    double* work;
    lwork = -1;
    char jobz='V';
    char matType = 'U';
    dsyev( &jobz, &matType, &nEnvOut, &z_[0], &lda, &eigenVal[0], &wkopt, &lwork, &info ); //with -1 this finds work size
    lwork = (int)wkopt;
    work = new double[lwork];
    // Solve eigenproblem. eigenvectors are stored in the z_ matrix, columnwise
    dsyev( &jobz, &matType, &nEnvOut, &z_[0], &lda, &eigenVal[0], work, &lwork, &info );
    bool status = ( info>0 || info<0 )? false : true;
    delete [] work;

    if (!status) {
#ifdef cqmom_dbg
      cout << "WARNING: Arches: CQMOMInversion: Solving Eigenvalue problem failed. Moment set: " << endl;
#endif
      nEnvOut--;
      continue;
    }
    //_____________
    //Solve actual weights and abscissas
    for( int i = 0; i < nEnvOut; i++) {
      w[i] = moments[0] * z_[i*nEnvOut] * z_[i*nEnvOut];
      x[i] = eigenVal[i];
    }

    //____________
    //Check that the minimum spacing and weight ratios are met
    vector<double> dab (nEnvOut);
    vector<double> mab (nEnvOut);
    double mindab, maxmab;

    double minTest;
    for (int i = nEnvOut-1; i>=1; i--) {
      dab[i] = fabs(x[i] - x[0]);
      for (int j = 1; j<i; j++) {
        minTest = fabs(x[i]-x[j]);
        if (minTest < dab[i]) {
          dab[i] = minTest;
        }
      }
    }
    double maxTest;
    for (int i = nEnvOut-1; i>=1; i--) {
      mab[i] = fabs(x[i] - x[0]);
      for (int j = 1; j<i; j++) {
        maxTest = fabs(x[i]-x[j]);
        if (maxTest > mab[i]) {
          mab[i] = maxTest;
        }
      }
    }

    mindab = dab[1]; maxmab = mab[1];
    for (int i=1; i<nEnvOut; i++) {
      if (dab[i] < mindab) {
        mindab = dab[i];
      }
      if (mab[i] > maxmab) {
        maxmab = mab[i];
      }
    }

    //check that prescribed condtions are met
    double maxW, minW;
    maxW = w[0];
    minW = w[0];
    for (int i = 0; i <nEnvOut; i++) {
      if (w[i] < minW) {
        minW = w[i];
      }
      if (w[i] > maxW) {
        maxW = w[i];
      }
    }
    if (minW/maxW > rMin && mindab/maxmab > eAbs) {
      isRealizable = true;
    } else {
      nEnvOut--;
#ifdef cqmom_dbg
      cout << "Weight Ratio: " << minW/maxW << " Abscissa Ratio: " << mindab/maxmab << endl;
      cout << "Decreasing environments " << nEnvOut+1 << " to " << nEnvOut << endl;
#endif
    }

  }
  //fill in any extra weight/abscissa values as zero
  if ( nEnvOut < nEnvMax ) {
    for (int i = nEnvOut; i < nEnvMax; i++ ) {
      w[i] = 0.0;
      x[i] = 0.0;
    }
  }
}
