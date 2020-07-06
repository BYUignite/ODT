/**
 * @file eddy.h
 * @brief Header file for class eddy
 */

#pragma once

#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing eddy object
 *
 *  @author David O. Lignell
 */

class eddy {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain              *domn;              ///< pointer to domain object
        domain              *eddl;              ///< pointer to eddy line object

        double              eddySize;           ///< size of eddy
        double              leftEdge;           ///< left edge location of eddy
        double              rightEdge;          ///< right edge location of eddy
        double              invTauEddy;         ///< inverse eddy timescale
        double              Pa;                 ///< eddy acceptance probability
        bool                LperiodicEddy;      ///< a wrap-around eddy
        vector<double>      cCoef;              ///< coefficient of K kernel
        vector<double>      bCoef;              ///< coefficient of J kernel
        vector<double>      K;                  ///< eddy kernel K
        vector<double>      dxc;                ///< \delta(x^cCoord) is prop. to cell "volume"
        vector<double>      pos0;               ///< initial eddy cell locations, for kernel

        double              esdp1;              ///< eddy size distribution parameters.
        double              esdp2;
        double              esdp3;
        double              esdp4;

        vector<double> cca,ccb,ccc,ccd;                 ///< polynomial coefficient arrays for cylindricalAnomalyHack

    //////////////////// MEMBER FUNCTIONS /////////////////

        void   sampleEddySize();
        void   sampleEddyPosition();
        void   tripMap(domain *line,    const int iS, int iE, const double C, const bool LsplitAtEddy=false);
        bool   eddyTau(const double Z_value, const double C);
        void   computeEddyAcceptanceProb(const double dtSample);
        void   applyVelocityKernels(domain *line, const int iS, const int iE);

    private:

        void   fillKernel();
        void   fillKernel_planarAnalytic();
        double eddyFavreAvgVelocity(const vector<double> &dxc);
        void   set_kernel_coefficients();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        eddy(){}
        void init(domain *p_domn, domain *p_eddl);
        ~eddy(){}

};


////////////////////////////////////////////////////////////////////////////////


