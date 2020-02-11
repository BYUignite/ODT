/**
 * @file micromixer.h
 * Header file for class micromixer
 */

#pragma once

#include <vector>
#include <string>
#include "cvodeDriver.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer object
 *
 *  @author David O. Lignell
 */

class micromixer {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain         *domn;          ///< pointer to domain object
        cvodeDriver    *cvode;         ///< pointer to cvode driver object for implicit ODE integration (stiff)

        double         tstart;
        double         time;           ///< current time
        double         tend;
        double         dtStepNominal;  ///< nominal step size
        double         dt;             ///< actual step size (shortened based on output or tend)

        vector<double> dxc;            ///< abs(\Delta(x^c))
        vector<double> dx;             ///< abs(\Delta(x))
        vector<double> gf;             ///< grid factor for derivatives: (df/dx) = gf * (f - f)

        bool           LdoDump;        ///<

        vector<double> uDL_1;          ///< for DL instability: old velocity
        vector<double> uDL_2;          ///< for DL instability: new velocity
        vector<double> xDL_1;          ///< for DL instability: = "old" cell center positions
        vector<double> xDL_2;          ///< for DL instability: = "new" cell center positions
        vector<double> posDL_old;      ///< for DL instability: = "new" cell center positions

        vector<double> oldrho_or_rhov; ///< store the old density for continuity

        int nsteps;                    ///< total number of timesteps taken during simulation

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void advanceOdt(const double p_tstart, const double p_tend, const int iLevel = -1);    // iLevel is for hips

        void check_balance(int io);

    protected:

        virtual void setGf();                ///< sets the gf array
        virtual void setGridDxcDx();         ///< sets the dxc array
        virtual void set_oldrho_or_rhov();   ///< record old rho (or rho*u) for continuity
        virtual bool adaptGridIfNeeded();   ///< expansion or contraction --> adapt
        virtual void setNominalStepSize();   ///< sets a nominal dt for the whole period

        void setStepSize();                  ///< set a local dt for interruptions (dump or tend)
        void updateGrid();                   ///< enforce the continuity condition: (e.g., rho*dx = const).
        void do_DL(string doWhat);

        void advanceOdtSingleStep_Explicit();
        void advanceOdtSingleStep_SemiImplicit();
        void advanceOdtSingleStep_StrangSplit();

        bool LforceSetNominalStepSize;       ///< used in updateGrid when splitting cells to indicate to reset timestep size later




    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer();
        void init(domain *p_domn);
        virtual ~micromixer(){ delete cvode; }

};


////////////////////////////////////////////////////////////////////////////////


