
#include <vector>

using namespace std;

//-------------------------------------------------------------------------------------------
/*
 * From Numerical Recipes in C page 92. G.B. Rybicki routine.
 * solve [x matrix] * w = q, where [x matrix] is formed from input vector x.
 * @param x \input  vector defining the Vandermonde matrix.
 * @param w \output vector containing the solution
 * @param q \input  vector containing the right had side
 * @param n \input  size of the system
 *
 */

void vandermonde_solver(const vector<double> &x, vector<double> &w, const vector<double> &q, const int &n) {

    int i,j,k;
    double b,s,t,xx;
    vector<double> c(n, 0.0);       // initialize to 0

    if(n==1)
        w[0] = q[0];

    else {
        c[n-1] = -x[0];
        for(i=1; i<n; i++) {
            xx = -x[i];
            for(j=(n-1-i); j<(n-1); j++)
                c[j] += xx*c[j+1];
            c[n-1] += xx;
        }
        for(i=0; i<n; i++) {
            xx = x[i];
            t = 1.0;
            b = 1.0;
            s = q[n-1];
            for(k=n-1; k>=1; k--) {
                b = c[k] + xx*b;
                s += q[k-1]*b;
                t = xx*t+b;
            }
            w[i] = s/t;
        }
    }
}
