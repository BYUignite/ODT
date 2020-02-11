/**
 * @file cvodeDriver.h
 * Header file for class cvodeDriver.h
 */

#pragma once

#include "dv.h"

class domain;

#include "cvode/cvode.h"
#include "cvode/nvector_serial.h"
#include "cvode/cvode_dense.h"
#include "cvode/sundials_dense.h"
#include "cvode/sundials_types.h"

#include <map>
#include <vector>
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////

/** Class integrates a single cell at a time.
 *  Intended for stiff systems (e.g. combustion chemistry).
 *
 *  @author David O. Lignell
 */

class cvodeDriver {

    public :

        ////////////////////// DATA MEMBERS /////////////////////

        domain              *domn;           ///< pointer to main domain object
        int                 neq;             ///< number of eqns solved
        N_Vector            var;             ///< vector of variables being solved for CVode
        map<int,dv*>        tVarMap;         ///< map to transported variables. (Domain vars in any order, but here, solve transported).
        int                 iC;              ///< which cell are we integrating
        double              atol;            ///< CVODE tol
        double              rtol;            ///< CVODE tol

        bool                LincludeRhsMix;  ///< if true, mixing term is included in integration

    private :

        void                *cvode_mem;      ///< CVode memory
        bool                Ldestruct;       ///< true if we setup cvode and can therefore destruct
        vector<double>      vard;            ///< variable array dummy
        vector<double>      Sd;              ///< variable source dummy


        ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public :

        void integrateCell(int p_iC, double tres);

    private :

        void testCVflag(int flag, string func);

    public :

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

        cvodeDriver(){Ldestruct = false;}                         // constructor
        void init(domain *p_domn, const bool p_LincludeRhsMix);  // initializer
        ~cvodeDriver();                                           // destructor

};


