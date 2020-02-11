/**
 * @file solver_hips.h
 * Header file for class solver_hips
 */

#pragma once

#include <vector>
#include <utility>    // pair
#include "solver.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver_hips object
 *
 *  @author David O. Lignell
 *
 */

class solver_hips : public solver {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        vector<double> levelRates;     ///< list of eddy event rates at each level
        vector<pair<double,int> > eTL; ///< list of eddy times and levels
        int iEta;                      ///< Kolmogorov level (needed for variable Sc scalars)
        double totalEddyRate;          ///< total rate of all eddies 0 through nLevels-3
        vector<double> eddyLevelCDF;   ///< sample from this to get eddy level

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();

    private:

        void setEddyEventTimes();
        void selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap);
        void reset_rates_for_Sc(const vector<double> &levelTaus);
        void sample_hips_eddy(double &dt, double &iLevel);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver_hips(){}
        virtual void init(domain *p_domn);

        virtual ~solver_hips(){}

};


////////////////////////////////////////////////////////////////////////////////


