/*
 * @file interp_linear.h
 * @brief Header file for class Linear_interp
 */

#ifndef INTERP_LINEAR_H
#define INTERP_LINEAR_H

#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class Linear_interp {

    public:

        vector<double> *X;
        vector<double> *Y;

        int nxy;

    private:

        int ilo;
        int ihi;

        //////////////////////////////////////////////////////////////////

    public:

        double interp(double x){
            set_bounding_indicies(x);
            return (*Y)[ilo] + (x-(*X)[ilo])*((*Y)[ihi]-(*Y)[ilo])/((*X)[ihi]-(*X)[ilo]);
        }

    private: 

        void set_bounding_indicies(double x){
            if(x <= (*X)[0])
                ihi = 1;
            else if(x >= (*X).back())
                ihi = nxy-1;
            else {
                vector<double>::iterator itHi = lower_bound((*X).begin(), (*X).end(), x); // lower_bound gives values >= x
                ihi = itHi - (*X).begin();
            }
            ilo = ihi-1;
        }

        //////////////////////////////////////////////////////////////////

    public:

        //------------- constructors

        Linear_interp(){}

        Linear_interp(vector<double> &X_p, vector<double> &Y_p){
            X = &X_p;
            Y = &Y_p;
            nxy = (*X).size();
            if((*X).size()!=(*Y).size()){
                cout << endl << "Error in interp_linear: X, Y need same size" << endl;
                exit(0);
            }

            for(int i=1; i<nxy-1; i++) 
                if( ((*X)[i] == (*X)[i-1]) || 
                    ((*X)[i] == (*X)[i+1]) || 
                    ((*X)[i]<(*X)[i-1] && (*X)[i]<(*X)[i+1]) ||
                    ((*X)[i]>(*X)[i-1] && (*X)[i]>(*X)[i+1]) ){
                    cout << endl << "Error in interp_linear: X should be monotonic" << endl;
                    exit(0);
                }
        }

};

#endif
