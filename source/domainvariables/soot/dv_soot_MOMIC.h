/**
 * @file dv_soot_MOMIC.h
 * Header file for class dv_soot_MOMIC
 */

#pragma once

#include "dv_soot.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_MOMIC of parent dv object.
 *
 *  @author Victoria B. Lansinger
 */

class dv_soot_MOMIC : virtual public dv_soot {

    //////////////////// DATA MEMBERS //////////////////////

    private:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1);

    private:

        double  lagrangeInterp(double x_i, vector<double> x, vector<double> y);
        double  MOMIC(double p, vector<double> M);
        double  f_grid(int x, int y, vector<double> M);
        double  beta(int p, int q, int ipt);
        double  getCoag(double T, double P, double mu, vector<double> M, int r);
        void    downselectIfNeeded(int ipt, vector<double> &M, int &N);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_MOMIC(domain      *line,
                     const string  s,
                     const bool    Lt,
                     const bool    Lo=true) : dv_soot(line, s, Lt, Lo) {}

        virtual ~dv_soot_MOMIC(){}

};


////////////////////////////////////////////////////////////////////////////////



