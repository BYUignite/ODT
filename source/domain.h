/**
 * @file domain.h
 * @brief Header file for class \ref domain
 */

#pragma once

#include "dv.h"
#include "domaincase.h"
#include "inputoutput.h"
#include "param.h"
#include "streams.h"
#include "micromixer.h"
#include "eddy.h"
#include "meshManager.h"
#include "solver.h"
#include "randomGenerator.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"
#include <vector>
#include <string>
#include <map>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing domain object
 *
 *  @author David O. Lignell
 */

class domain {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain                  *domn;     ///< (for one domain to point to another (eddl))

        int                     ngrd;      ///< number of grid cells
        int                     ngrdf;     ///< number of grid cell faces = ngrd+1

        vector<dv*>             v;         ///< All domain variables are stored in here.

        dv*                     pos;       ///< pointers to gas properties
        dv*                     posf;      ///< access as: posf->d[i], or posf->var_name, etc.
        dv*                     rho;
        dv*                     dvisc;
        dv*                     uvel;
        dv*                     vvel;
        dv*                     wvel;
        dv*                     sdiff;
        dv*                     sca;
        dv*                     phase;
        dv*                     enth;
        dv*                     temp;
        dv*                     mixf;
        dv*                     chi;
        dv*                     hr;
        dv*                     aDL;
        vector<dv*>::iterator   ysp;       ///< access as: ysp=v.begin(), (*ysp)->d[i] or (*(ysp+k))->d[i], or ysp[k]->d[i].
        vector<dv*>::iterator   svar;      ///< iterator for increment to go through moments (*(ysp+k))->d[i];)
        vector<dv*>::iterator   eta;       ///< iterator for increment to go through species etc. (*(ysp+k))->d[i];)

        map<string,dv*>         varMap;

        IdealGasPhase           *gas;        ///< pointer to cantera thermochemistry object (reaction rates, Cp, etc.)
        Transport               *tran;       ///< pointer to cantera transport object (viscosity, diffusivity, etc.)
        streams                 *strm;       ///< pointer to gas stream properties
        inputoutput             *io;         ///< pointer to input/output object
        param                   *pram;       ///< pointer to the parameters object
        micromixer              *mimx;       ///< pointer to micromixer for diffusion, reaction, domain evolution.
        eddy                    *ed;         ///< pointer to object for eddy operations
        domain                  *eddl;       ///< pointer to eddyline object
        solver                  *solv;       ///< pointer to solver object
        meshManager             *mesher;     ///< pointer to mesh manager object

        randomGenerator         *rand;

        int                     nTrans;      ///< number of transported variables on the domain.

        domaincase              *domc;       ///< domaincase class: set specific vars...


    //////////////////// MEMBER FUNCTIONS /////////////////

        int    domainPositionToIndex(double position, const bool LowSide, int dbg);
        void   setDomainFromRegion(const int i1, const int i2);
        double cyclePeriodicDomain(const int icycle);
        void   backCyclePeriodicDomain(const double backCycleDistance);
        double Ldomain();

    private:

        void initEddyDomain();


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        void init(inputoutput     *p_io,
                  meshManager     *p_mesher,
                  streams         *p_strm,
                  IdealGasPhase   *p_gas,
                  Transport       *p_tran,
                  micromixer      *p_mimx,
                  eddy            *p_ed,
                  domain          *p_eddl,
                  solver          *p_solv,
                  randomGenerator *p_rand,
                  bool             LisEddyDomain=false);
        domain(domain *p_domn, param *p_pram);
        virtual ~domain() {
            for(int k=0; k<v.size(); k++)
                delete v.at(k);
            delete domc;     
        }

};


////////////////////////////////////////////////////////////////////////////////


