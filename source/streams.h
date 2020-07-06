/**
 * @file streams.h
 * @brief Header file for class streams
 */

#pragma once

#include <vector>

using namespace std;

class domain;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing streams for use in mixing and or reaction problems.
 *  This is writting in terms of mixture fraction with streams defined in an
 *  input file.  The class can implement products of complete combustion, or
 *  equilibrium (through the Cantera IdealGasMix object, if desired).
 *  This class holds a pointer to a Cantera IdealGasMix object (defined up front
 *  in main) that computes thermodynamic, kinetic, and transport data.
 *
 *  @author David O. Lignell
 */

class streams {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain         *domn;           ///< pointer to domain object

        double         T0;              ///< stream mixf=0 temperature
        double         T1;              ///< stream mixf=1 temperature
        double         h0;              ///< stream mixf=0 enthalpy
        double         h1;              ///< stream mixf=1 enthalpy
        vector<double> y0;              ///< stream mixf=0 composition vector
        vector<double> y1;              ///< stream mixf=1 composition vector
        double         M0;              ///< stream mixf=0 mean molecular weight
        double         M1;              ///< stream mixf=1 mean molecular weight
        double         rho0;            ///< stream mixf=0 density
        double         rho1;            ///< stream mixf=1 density
        vector<double> D0;              ///< stream mixf=0 diffusivities 
        vector<double> D1;              ///< stream mixf=1 diffusivities
        vector<double> hsp0;            ///< stream mixf=0 species enthalpies
        vector<double> hsp1;            ///< stream mixf=1 species enthalpies

        double         mixfStoic;       ///< stoichiometric mixture fraction

        int            nspc;            ///< number of species in gas mechanism

                                        /// \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double         beta0;           ///< mixf = (beta-beta0) / (beta1-beta0)
                                        ///< \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double         beta1;           ///< mixf = (beta-beta0) / (beta1-beta0)

        vector<double> gCHON;           ///< gammas, as in beta = sum_i (y_i*gamma_i)

        int            comp2;           ///< for premixed combustion to distinguish between different input possibilities for the mixture

    //////////////////// MEMBER FUNCTIONS /////////////////

        void getProdOfCompleteComb(const double mixf,
                                   vector<double> &ypcc,
                                   double &hpcc,
                                   double &Tpcc);

        void getMixingState(const double mixf,
                            vector<double> &ymix,
                            double &hmix,
                            double &Tmix);

        double getMixtureFraction(const double *y,
                                  const bool doBeta01=false);

    private:

        void setStoicMixf();
        vector<double> setElementMassFracs(const double *y);
        vector<double> setElementMoleFracs(const double *y);
        vector<double> getElementMoles(const double *x,
                                       double &nOnotFromO2,
                                       double &nHnotFromH2O,
                                       double &nCnotFromCO2);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        streams(){}
        void init(domain *p_domn, const vector<double> &gammas);

        ~streams(){}




};



