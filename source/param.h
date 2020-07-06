/**
 * @file param.h
 * @brief Header file for class param
 */

#pragma once

#include "inputoutput.h"
#include <string>
#include <cstdlib>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class param {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain  *domn;           ///< pointer to domain object
        inputoutput *io;         ///< pointer to io object (has the input file)

        int     seed;            ///<  random number generator seed (negative to randomize it)
        double  tEnd;            ///<  ending time of realization
        double  domainLength;    ///<  length of domain (m)
        int     ngrd0;           ///<  initial grid points
        double  rho0;            ///<  initial uniform density (kg/m^3)
        double  kvisc0;          ///<  initial uniform kinematic viscosity (m^2/s)
        double  sdiff0;          ///<  initial uniform scalar diffusivity (m^2/s)
        double  dPdx;            ///<  initial pressure gradient (Pa/m)
        double  pres;            ///<  initial pressure (Pa)
        string  chemMechFile;    ///<  name of chemical mechanism file
        string  probType;        ///<  problem type: CHANNEL, CHANNEL_SCALAR, JETMIXL_RXN, COUETTE

        double  Z_param;         ///<  Viscous penalty parameter
        double  A_param;         ///<  Energy Distribution parameter alpha
        double  C_param;         ///<  Eddy frequency parameter
        string  LES_type;        ///<  NONE, THIRDS, ELAPSEDTIME, FRACDOMAIN, INTEGRALSCALE
        double  Z_LES;           ///<  large eddy suppression (nonpositive prevents les test)
        double  diffCFL;         ///<  multiplies min diffusion timestep
        double  cvode_atol;      ///<  absolute tolerace atol for cvode
        double  cvode_rtol;      ///<  relative tolerace rtol for cvode
        double  x0virtual;       ///<  LES virtual origin

        bool    LdoDL;           ///<  flag to do the DL energy from the DL instability
        bool    Lrad;            ///<  radiation flag
        bool    Lbuoyant;        ///<  flag to turn on bouyancy (horizontal domain)
        bool    LPeEddy;         ///<  flag to turn on potential energy for eddies (vertical domain)
        bool    LplanarExpCent0; ///<  flag: for planar cases (C=1) set the expansion center at 0 for outflow cases (normally expand about the expansion center.
        double  g;               ///<  gravity (default -9.81)
        string  Lsolver;         ///<  EXPLICIT, SEMI-IMPLICIT, or STRANG
        bool    Lperiodic;       ///<  periodic if true
        bool    Lspatial;        ///<  spatial formulation if true
        bool    LTMA;            ///<  true for the triplet map TMA: 3 = vol segments; false for TMB: 3 equal length segments
        bool    LplanarTau;      ///<  true for computing cylindrical/spherical tau_eddy using a planar formulation. If accepted, a cylindrical eddy is implemented
        bool    Lignition;        ///<  true if starting with unreacted mixing profile to allow ignition

        string  bcType;          ///<  OUTFLOW, PERIODIC, WALL, WALL_OUT
        int     cCoord;          ///<  1 = planar, 2 = cylindrical, 3 = spherical
        double  xDomainCenter;   ///<  position of the center of the domain

        double  gDens;           ///<  grid density for mesher
        double  dxmin;           ///<  min grid spacing: = dxmin / domain length
        double  dxmax;           ///<  max grid spacing = dxmax / domain length

        double  Pmax;            ///<  maximum eddy acceptance probability
        double  Pav;             ///<  Average acceptance probability
        double  dtfac;           ///<  maximum factor to increase dtSmean
        int     nDtSmeanWait;    ///<  number of eddy samples before increase dtSmean
        int     eddyMinCells;    ///<  eddy must overlap at least this many cells
        double  DAtimeFac;       ///<  time until catch-up adaption is DAtimeFac * dtCUmax
        double  tdfac;           ///<  factor between dtCUmax and dtCFL for temporal flows; DEFAULT = 1.0
        int     sLastDA;         ///<  size of the lastDA vector for timing adaptmesh after diff
        double  Lp;              ///<  Most probable eddy size frac of domainLength
        double  Lmax;            ///<  Max eddy size frac of domainLength
        double  Lmin;            ///<  Min eddy size frac of domainLength

        int     modDump;         ///<  accepted eddies before output file
        int     modDisp;         ///<  frequency to display results (# eddies)
        bool    Ltecplot;        ///<  set TRUE for tecplot friendly output

        bool    LmultiPhase;     ///<  true if domain has more than one phase (particles don't count)
        double  eSurfTens;       ///<  surface tension, J/m2 for liquid phases

        double  uBClo;           ///<  Dirichlet velocity boundary condition.
        double  uBChi;           ///<  Dirichlet velocity boundary condition.
        double  vBClo;           ///<  Dirichlet velocity boundary condition.
        double  vBChi;           ///<  Dirichlet velocity boundary condition.
        double  wBClo;           ///<  Dirichlet velocity boundary condition.
        double  wBChi;           ///<  Dirichlet velocity boundary condition.
        double  sBClo;           ///<  Dirichlet scalar boundary condition.
        double  sBChi;           ///<  Dirichlet scalar boundary condition.
        string  hWallBCtype;     ///<  ADIABATIC or ISOTHERMAL
        double  TBClo;           ///<  Required if hWallBCtype = ISOTHERMAL
        double  TBChi;           ///<  Required if hWallBCtype = ISOTHERMAL

        bool    Lrestart;        ///<  true to restart from file, else false
        string  rstType;         ///<  "single" or "multiple"
        double  trst;            ///<  restart time (from restart file), default is 0.0;

        double  umin_spatial;    ///< min u for spatial flows; used when kernels pull velocity

        //----------------- Radiation variables

        string                      radSolveType;   ///< OPTHIN, TWOFLUX, FVDOM
        string                      radCoefType;    ///< PLANCKMEAN, WSGG, RCSLW
        int                         npsi;
        int                         ntheta;

    //////////////////// MEMBER FUNCTIONS /////////////////

    private:

        template <class T>
        T errMsg(const string param) {
            *io->ostrm << endl << "ERROR: missing parameter: " + param << endl;
            exit(0);
            T dummy = static_cast<T> (0);
            return dummy;
        }

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        param(inputoutput *p_io);
        void init(domain *p_domn);
        ~param(){}

};


////////////////////////////////////////////////////////////////////////////////


