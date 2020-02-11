#ifndef INTERP_LINEAR_H
#define INTERP_LINEAR_H

#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class Linear_interp {

    public:

        vector<double> X;
        vector<double> Y;

        int nxy;

        //////////////////////////////////////////////////////////////////

        double interp(double x);

    private: 

        int get_location(double x);

        //////////////////////////////////////////////////////////////////

    public:

        Linear_interp(){}

        Linear_interp(vector<double> &X_p, vector<double> &Y_p){
            X = X_p;
            Y = Y_p;
            nxy = X.size();
            if(X.size()!=Y.size()){
                cout << endl << "Error in interp_linear: X, Y need same size" << endl;
                exit(0);
            }

            for(int i=1; i<nxy-1; i++) 
                if( (X[i] == X[i-1]) || 
                    (X[i] == X[i+1]) || 
                    (X[i]<X[i-1] && X[i]<X[i+1]) ||
                    (X[i]>X[i-1] && X[i]>X[i+1]) ){
                    cout << endl << "Error in interp_linear: X should be monotonic" << endl;
                    exit(0);
                }
        }

};

#endif
