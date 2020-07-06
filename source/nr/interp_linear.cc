/**
 * @file interp_linear.cc
 * @brief Source file for class interp_linear
 */

#include "interp_linear.h"

/////////////////////////////////////////////////////////////////////

double Linear_interp::interp(double x){

    int ilo = get_location(x);
    int ihi = ilo+1;
    if(ilo==nxy-1)     // if x coincides with the last point
        ihi = ilo-1;   //    then interpolate using the last two points.

    return Y[ilo] + (x-X[ilo])*(Y[ihi]-Y[ilo])/(X[ihi]-X[ilo]);
}

/////////////////////////////////////////////////////////////////////

/* Desired point x will normally be between two grid points.
   Find the lower grid point index.
   If x coincides with a grid point, then return that index.
*/
   
int Linear_interp::get_location(double x){

    int ilo = 0;
    int ihi = nxy-1;

    if(x <= X[ilo])
        return ilo;
    if(x >= X[ihi])
        return ihi;

    int im, ip;

    for(int iter=0; iter<10000; iter++){
        //im = static_cast<int>((ilo + ihi)/2);                            // midpoint reduction
        //im = static_cast<int>(ilo + (ihi-ilo)*(x-X[ilo])/(X[ihi]-X[ilo])); // uniform grid assumption
        im = ilo + static_cast<int>((ihi-ilo)*(x-X[ilo])/(X[ihi]-X[ilo])); // uniform grid assumption
        ip = im + 1;
        if(X[im] <= x && x < X[ip])
            return im;
        else if(x==X[ip])
            return ip;
        else if (x > X[ip])
            ilo = ip;
        else
            ihi = im;
        if(iter==9999){
            cout << endl << "Error in interp_linear::get_location: no convergence" << endl; 
            exit(0);
        }
    }
    return im;
}
