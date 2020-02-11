/**
 * @file dv_soot.h
 * Header file for class dv_soot
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot of parent dv object.
 *  This is a virtual base class.
 *  @author Victoria B. Lansinger
 */

class dv_soot : public dv {

    //////////////////// DATA MEMBERS //////////////////////


    protected:

        static int              N;                      ///< iterator for assigning kMe values

        int                     kMe;                    ///< kMe
        static int              nsvar;                  ///< number of soot variables

        //-----------

        static constexpr double Na    = 6.02214086E26;  ///< Avogadro's constant: #/kmol
        static constexpr double kb    = 1.38064852E-23; ///< Boltzmann constant = Rg/Na: J/#*K
        static constexpr double Rg    = 8314.46;        ///< Universal gas constant
        static constexpr double eps_c = 2.2;            ///< coagulation constant
        static constexpr double Df    = 1.8;            ///< soot fractal dimension

        //-----------

        static double           Cmin;                   ///< number of carbons in soot nucleation size
        static double           MW_c;                   ///< mw of carbon 12.011;
        static double           MW_h;                   ///< mw of hydrogen 1.0079;
        static double           rhoSoot;                ///< soot density kg/m3
        static double           b_coag;                 ///< coagulation rate parameter (see Lignell thesis page 58.)
        static string           nucleation_mech;        ///< soot nucleation chemistry flag
        static string           growth_mech;            ///< soot growth chemistry flag
        static string           oxidation_mech;         ///< soot oxidation chemistry flag
        static string           coagulation_mech;       ///< soot coagulation mechanism flag

        static double           DIMER;                  ///< dimer concentration
        static double           m_dimer;                ///< dimer mass

        //-----------

        static double           rC2H2_rSoot_n;          ///<
        static double           rH2_rSoot_ncnd;         ///< mass rate ratio: for gasSootSources from nucleation
        static vector<double>   rPAH_rSoot_ncnd;        ///< mass rate ratio: for gasSootSources from nucleation/condensation
        static double           rO2_rSoot_go;           ///< mass rate ratio: for gasSootSources from growth/oxidation
        static double           rOH_rSoot_go;           ///<
        static double           rH_rSoot_go;            ///<
        static double           rCO_rSoot_go;           ///<
        static double           rH2_rSoot_go;           ///<
        static double           rC2H2_rSoot_go;         ///<

        //----------- gas state variables

        static double           T;                      ///< K
        static double           P;                      ///< Pa
        static double           rho;                    ///< kg/m3
        static double           MW;                     ///< kg/kmol mean molecular weight
        static double           mu;                     ///< kg/m*s
        static vector<double>  *yi;                     ///< pointer to species mass fractions

        //-----------

        static int              i_c2h2;                 ///< soot gas species indices
        static int              i_o2;
        static int              i_h;
        static int              i_h2;
        static int              i_oh;
        static int              i_h2o;
        static int              i_co;
        static int              i_elem_c;               ///< index if element C
        static int              i_elem_h;               ///< index if element H
        static vector<int>      i_pah;                  ///< vector of PAH species indicies
        static vector<double>   MW_sp;                  ///< vector of molecular weights


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void getRhsSrc(const int ipt=-1) = 0;
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    protected:

        void   set_gas_state_vars (const double &T_p, const double &P_p, const double &rho_p, const double &MW_p, const double &mu_p, vector<double> &y_p);

        double getNucleationRate  (const vector<double> &mi=vector<double>(0), const vector<double> &wi=vector<double>(0));
        double getGrowthRate      (const double &M0=-1, const double &M1=-1);
        double getOxidationRate   (const double &M0=-1, const double &M1=-1);
        double getCoagulationRate (const double &m1,    const double &m2);

        double get_gas_mean_free_path();
        double get_Kc();
        double get_Kcp();
        double get_Kfm();

        double set_m_dimer();
        void   set_Ndimer(const vector<double> &mi, const vector<double> &wi);

        void   set_gasSootSources(const double &N1, const double &Cnd1, const double &G1, const double &X1, const int igrd);

    private:

        double nucleation_LL      ();
        double nucleation_Linstedt();
        double nucleation_PAH     (const vector<double> &mi, const vector<double> &wi);

        double growth_Lindstedt   ();
        double growth_LL          (const double &M0, const double &M1);
        double growth_HACA        (const double &M0, const double &M1);

        double oxidation_LL       ();
        double oxidation_Lee_Neoh ();
        double oxidation_NSC_Neoh ();
        double oxidation_HACA     (const double &M0, const double &M1);

        double coagulation_LL     (const double &m1, const double &m2);
        double coagulation_Fuchs  (const double &m1, const double &m2);
        double coagulation_Frenk  (const double &m1, const double &m2);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot(domain      *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_soot(){}

};

////////////////////////////////////////////////////////////////////////////////



