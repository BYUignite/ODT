#include<vector>
#include <iostream>
#include<iomanip>
#include <cmath>
//#include"lapack.h"
using namespace std;

/**--------------------------------------------------------------------------------
 * David Lignell
 * Gordon's product difference algorithm
 * From Gordon, 1968, and Rober McGraw 1984.  See also Numerical Recipies.
 * Originally coded in matlab: pd.m
 * There is also pd2pt.f90 which is a hardcoded 2pt (4 mom) quadrature
 * Uses imtql2 from eispack for eigenvals, vecs of symmetric tridiag mat
 * @param nm input number of moments
 * @param np input number of points (=nm/2)
 * @param mu input
 * @param wts output
 * @param absc output
 */

void imtql2(int                               nm, //nmd
        int                               n,  //np
        std::vector<double>               &d,  //absc
        std::vector<double>               &e,  //b
        std::vector<std::vector<double> > &z,  //P
        int                               &ierr);

double pythago(double a, double b);

void pdAlg(int nm, int np, vector<double> &mu, vector<double> &wts, vector<double> &absc ) {

    //for (int i=0; i<4; i++)   //Finite moments are fed to PD but it crashes =/
    //cout<<setprecision(17)<< mu[i]<< "  ";
    //cout << "\n";
    int                        nmd=nm+1;    // increase for nm>8
    vector<double>             M(nm);
    vector<vector<double> >    P(nmd, vector<double>(nmd, 0.0));
    vector<double>             alph(nm);
    vector<double>             a(np);
    vector<double>             b(np); // Eigen Vector?
    int                        ierr;
    double                     m0, m1;    // scale factors

    //--------------------------------------------------------------------------



    M = mu;   // This assignment operator works for vectors
    m0 = M[0];    // scale all mom by m0

    for(int i=0; i<nm; i++)
        M[i]  = M[i]/m0;              //loop inserted for vector math


    m1 = M[1];                        // then rescale for the absc (m1,2,3 have units)
    for(int i=1; i<nm; i++)
        M[i] = M[i] / pow(m1, i);


    P[0][0] = 1.0;
    for(int i=0; i<nm; i++)
        P[i][1] = M[i];

    for(int i=1; i<=nm; i+=1)
        P[i-1][1] = P[i-1][1]*pow(-1.0,i-1.0);



    for(int j=3; j<=nm+1; j++)
        for(int i=1; i<=nm-j+3-1; i++)  // 2 to 1
            P[i-1][j-1] = P[0][j-2]*P[i][j-3] - P[0][j-3]*P[i][j-2];

    alph[0] = 0.0;
    for(int i=1; i<nm; i++){
        alph[i] = P[0][i+1]/(P[0][i]*P[0][i-1]);
    }

    for(int i=1; i<=np-1; i++){
        absc[i-1] = alph[2*i-1] + alph[2*i-2];
        b[i] = sqrt(abs(alph[2*i] * alph[2*i-1]));

    }
    absc[np-1] = alph[2*np-1] + alph[2*np-2];

    P = vector<vector<double> >(nmd, vector<double>(nmd,0));  // What is this?
    for(int i=0; i<np; i++)
        P[i][i] = 1.0;

    imtql2( nmd, np, absc, b, P, ierr);           // a --> eigenvalues = absc

    // reuse P --> eigenvec
    if(ierr != 0)
        cout << endl << "Error in dstev for pd alg" << endl;

    for(int i=0; i<np; i++) {
        absc[i] = absc[i] * m1;
        wts[i] = P[0][i]*P[0][i]*m0;
    }

}

//================================================================================

void  imtql2(int nm,          //nmd
        int  n,         //np
        vector<double> &d, //absc
        vector<double> &e, //b
        vector<vector<double> > &z, //P
        int &ierr){ // ierr

    int j,mml;
    double b,c,f,g,p,r,s,tst1,tst2;
    // Assume initial values are all equal to zero
    j=0;   mml =0;
    b= 0.0; c=0.0; f=0.0; g=0.0; p=0.0;r=0.0;s=0.0;tst2=0.0;
    //cout<< e[0] <<" " << e[1]<<" "<<e[2] <<"\n  ";
    /*  !
        !     this subroutine is a translation of the algol procedure imtql2,
        !     num. math. 12, 377-383(1968) by martin and wilkinson,
        !     as modified in num. math. 15, 450(1970) by dubrulle.
        !     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
        !
        !     this subroutine finds the eigenvalues and eigenvectors
        !     of a symmetric tridiagonal matrix by the implicit ql method.
        !     the eigenvectors of a full symmetric matrix can also
        !     be found if  tred2  has been used to reduce this
        !     full matrix to tridiagonal form.
        !
        !     on input
        !
        !        nm must be set to the row dimension of two-dimensional
        !          array parameters as declared in the calling program
        !          dimension statement.
        !
        !        n is the order of the matrix.
        !
        !        d contains the diagonal elements of the input matrix.
        !
        !        e contains the subdiagonal elements of the input matrix
        !          in its last n-1 positions.  e(1) is arbitrary.
        !
        !        z contains the transformation matrix produced in the
        !          reduction by  tred2, if performed.  if the eigenvectors
        !          of the tridiagonal matrix are desired, z must contain
        !          the identity matrix.
        !
        !      on output
        !
        !        d contains the eigenvalues in ascending order.  if an
        !          error exit is made, the eigenvalues are correct but
        !          unordered for indices 1,2,...,ierr-1.
        !
        !        e has been destroyed.
        !
        !        z contains orthonormal eigenvectors of the symmetric
        !          tridiagonal (or full) matrix.  if an error exit is made,
        !          z contains the eigenvectors associated with the stored
        !          eigenvalues.
        !
        !        ierr is set to
        !          zero       for normal return,
        !          j          if the j-th eigenvalue has not been
        !                     determined after 30 iterations.
        !
        !     calls pythago for  dsqrt(a*a + b*b) .
        !
        !     questions and comments should be directed to burton s. garbow,
        !     mathematics and computer science div, argonne national laboratory
        !
        !     this version dated august 1983.
        !
        !     ------------------------------------------------------------------
        */

    int i=0;
    int m=0 ;
    ierr = 0;
    if (n == 1)
        return;
    for( int i=2; i<=n; i++){// loop to 100   // LSF
        e[i-2] = e[i-1];                // 100   conveinent to renumber the elements of e
    }
    e[n-1] = 0.0;  //0.0e0                  // LSF

    for(int l=1; l<=n;  l++){ // go to 240
        j = 0;

        //     .......... look for small sub-diagonal element ..........
m105:
        for( m = l; m<=n ; m++){   // marked as 105 go to 110 m to n
            if (m == n){  // go to 120   //reduncant with for loop
                break;

            }
            tst1 = abs(d[m-1]) + abs(d[m]); //using cmath's abs in place of dabs
            tst2 = tst1 + abs(e[m-1]);
            if (tst2 == tst1)  //  go to 120
                break;
        }         //marked as 110
        p = d[l-1]; //marked as 120


        if (m == l)   //go to 240
            continue;    // This is the only way to exit the loop correctly, without throwing an error

        if (j == 30){  // go to 1000
            ierr=l;
            cout<< " j==30\n ";
            return;
        }
        j = j + 1;
        //     .......... form shift ..........
        g = (d[l] - p) / (2.0 * e[l-1]); // 2.0d
        r = pythago(g,1.0);           //what!?  huh!?
        if (g<0)
            g = d[m-1] - p + e[l-1] / (g -abs(r));   //uses dsign(r,g) here instead of if statement
        else
            g = d[m-1] - p + e[l-1] / (g +abs(r));   //uses dsign(r,g) here instead of if statement

        s = 1.0;   // removed d0 s on each of these 3 lines
        c = 1.0;
        p = 0.0;
        mml = m - l;
        //     .......... for i=m-1 step -1 until l do -- ..........

        for (int ii = 1; ii<=mml; ii++) {  // go to 200 1 to mml steps
            i = m - ii;
            f = s * e[i-1];
            b = c * e[i-1];
            r = pythago(f,g);
            e[i] = r;
            if (r == 0.0){    // go to 210

                //     .......... recover from underflow ..........
                d[i] = d[i] - p;  // marked as 210
                e[m-1] = 0.0;
                goto m105;
            }

            s = f / r;
            c = g / r;
            g = d[i] - p;
            r = (d[i-1] - g) * s + 2.0 * c * b;  // 2.0d0
            p = s * r;
            d[i] = g + p;
            g = c * r - b;
            //    .......... form vector ..........
            for (int  k = 1; k<=n; k++){    // loop 1 to n  to 180
                f = z[k-1][i];
                z[k-1][i] = s * z[k-1][i-1] + c * f;// fix row, colomns to columns, rows
                z[k-1][i-1] = c * z[k-1][i-1] - s * f;
            } // marked as 180

        }//marked as 200
        d[l-1] = d[l-1] - p;
        e[l-1] = g;
        e[m-1] = 0.0; //d0
        goto m105;



    }  // marked as 240

    //     .......... order eigenvalues and eigenvectors ..........
    int k=0;// prviously k was only used as a counter in a loop

    for (int ii=2;ii<=n; ii++){   // go from 2 to n to 300

        i = ii - 1;
        k = i;
        p = d[i-1];
        //
        for ( j = ii; j<=n; j++ ){  //  260 loop from ii to n to 260
            if (d[j-1] >= p)// go to 260
                continue;
            k = j;
            p = d[j-1];
        } // marked as  260
        //
        if (k == i)    //if k ==i goto 300
            continue ;

        d[k-1] = d[i-1];
        d[i-1] = p;
        //
        for ( j = 1;j<=n;j++) {    // 1 to n to 280
            p = z[j-1][i-1];
            z[j-1][i-1] = z[j-1][k-1];
            z[j-1][k-1] = p;
        }// 280    continue

    }// Marked as 300
    //
    //     .......... set error -- no convergence to an
    //                eigenvalue after 30 iterations ..........
    // cout<< ierr<<" normal exit\n";
    return;
}
//================================================================================

double pythago(double a, double b) {
    double p=0,r=0,s=0,t=0,u=0 ;
    p = (fabs(a) > fabs(b)) ? a : b;//abs to fabs
    if(p==0)
        return p;

    r = (abs(a) < abs(b)) ? a : b;
    r = r*r/(p*p);
    for(;;) {
        t = 4.0+r;
        if(t==4.0)
            return (p);
        s = r/t;
        u = 1.0 + 2.0*s;
        p = u*p;
        r = s*s*r/(u*u);
        if(std::isnan(t)){
            cout<<"Nans in PDalg" << t<< "  \n";
            int delme;
            cin >> delme;}
    }
}

