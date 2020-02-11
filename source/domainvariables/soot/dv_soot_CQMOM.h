/**
 * @file dv_soot_CQMOM.h
 * Header file for class dv_soot_CQMOM
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_CQMOM of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_CQMOM : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<vector<double> >     wts;        ///< weights
        vector<vector<double> >     absc_v;     ///< abscissas, v direction
        vector<vector<double> >     absc_s;     ///< abscissas, s direction

        int                         nused_v;    ///<
        vector<int>                 nused_s;    ///<

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double Mkj(double v, double s);
//        void adaptive_CQMOM(int N_v, int N_s, const vector<vector<double> > &m,
//                    const double &eabs_v, const double &eabs_s,
//                    const vector<double> &rmin_v, const vector<double> &rmin_s,
//                    int &Nuse_v, vector<int> &Nuse_s, vector<vector<double> > &W,
//                    vector<vector<double> > &X_v, vector<vector<double> > &X_s );
//        void vandermonde_solver(const vector<double> &x, vector<double> &w, const vector<double> &q, const int &n);
//        void adaptive_wheeler(const vector<double> &m, int N, const vector<double> &rmin, const double &eabs,
//                      int &Nout, vector<double> &w, vector<double> &x );

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_CQMOM(domain      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {}

        virtual ~dv_soot_CQMOM(){}

};


////////////////////////////////////////////////////////////////////////////////



