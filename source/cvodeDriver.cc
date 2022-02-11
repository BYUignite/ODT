/**
 * @file cvodeDriver.cc
 * @brief Source file for class cvodeDriver
 */

#include <cstdlib>

#include "cvodeDriver.h"
#include "domain.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototypes to be used in this source file

static int RHSF(double t, N_Vector y, N_Vector ydot, void* f_data);

///////////////////////////////////////////////////////////////////////////////
/** Constructor function.
 *
 *  @param p_domn \input pointer to domain object
 */
void cvodeDriver::init(domain *p_domn, bool p_LincludeRhsMix) {

    //---------- set pointers

    domn           = p_domn;
    LincludeRhsMix = p_LincludeRhsMix;

    //----------

    if(domn->pram->Lsolver=="EXPLICIT")
        return;

    Ldestruct = true;

    //---------- set some CVode params

    atol = domn->pram->cvode_atol; // 1.0E-10;
    rtol = domn->pram->cvode_rtol; // 1.0E-4;

    //---------- set the number of equations to solve, and the variable map

    neq = domn->nTrans;
    for(int i=0, k=0; i<domn->v.size(); i++)
        if(domn->v.at(i)->L_transported)
            tVarMap[k++] = domn->v.at(i);

    //---------- initialize the dependent variable

    var = N_VNew_Serial(neq);

    //---------- set the CVode object

    int flag;

    cvode_mem = CVodeCreate(CV_BDF);

    if(!cvode_mem) {
        cout << endl << "ERROR INITIALIZING CVODE MEMORY" << endl;
        exit(0);
    }

    flag = CVodeInit(cvode_mem, RHSF, 0.0, var);
        testCVflag(flag, "CVodeInit");

    flag = CVodeSStolerances(cvode_mem, rtol, atol);
        testCVflag(flag, "CVodeSStolerances");

    SUNMatrix A = SUNDenseMatrix(neq, neq);
    SUNLinearSolver LS = SUNLinSol_Dense(var, A);

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
        testCVflag(flag, "CVodeSetLinearSolver");

    flag = CVodeSetUserData(cvode_mem, this);
        testCVflag(flag, "CVodeSetUserData");

    flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
        testCVflag(flag, "CVodeSetMaxNumSteps");
}

///////////////////////////////////////////////////////////////////////////////
/** Destructor */
cvodeDriver::~cvodeDriver() {
    if(Ldestruct) {
        N_VDestroy_Serial(var);
        CVodeFree(&cvode_mem);
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Public interface: Integrates a given cell (iC) for time tres.
 *  @param iC   \input integrate this cell.
 *  @param tres \input time to integrate for.
 */
void cvodeDriver::integrateCell(int p_iC, double tres) {

    iC = p_iC;
    int flag;

    //---------- Initialize the Dependent Variable

    for(int k=0; k<neq; k++)
        NV_Ith_S(var, k) = tVarMap[k]->d.at(iC);

    //---------- Reset CVode for the new cell

    flag = CVodeReInit(cvode_mem, 0.0, var);
        testCVflag(flag, "CVodeReInit");

    //---------- Integrate the solution

    double t;
    flag = CVodeSetStopTime(cvode_mem, tres);
        testCVflag(flag, "CVodeSetStopTime");
    flag = CVode(cvode_mem, tres, var, &t, CV_NORMAL);
        testCVflag(flag, "CVode");

    //----------- Recover the solution

    for(int k=0; k<neq; k++)
        tVarMap[k]->d.at(iC) = NV_Ith_S(var,k);

}

///////////////////////////////////////////////////////////////////////////////
/** Right hand side function for implicit integration
 *  f_data is cast to the cvodeDriver that calls cvode (so cvode is passed pointer "this").
 *  "static" means only defined in this file (just following cvode standard use).
 *
 *  @param t      \input, current time (from start of integration).
 *  @param y      \input, current dependent variable state during integration.
 *  @param ydot   \output, rhs of dy/dt = ydot = explict+implict terms. (Rate for the t integration.)
 *  @param f_data \input, pointer to be cast to user data for computing ydot from y.
 */
static int RHSF(double t, N_Vector y, N_Vector ydot, void* f_data) {

    cvodeDriver *cvd;
    cvd = static_cast<cvodeDriver*>(f_data);

    for(int k=0; k<cvd->neq; k++)                         // set the transported vars on the domain
        cvd->tVarMap[k]->d.at(cvd->iC) = NV_Ith_S(y,k);

    cvd->domn->domc->setCaseSpecificVars_cvode(cvd->iC);

    for(int k=0; k<cvd->neq; k++)                    // set the transported var souce terms
        cvd->tVarMap[k]->getRhsSrc(cvd->iC);

    for(int k=0; k < cvd->neq; k++) {
        NV_Ith_S(ydot,k) = cvd->tVarMap[k]->rhsSrc.at(cvd->iC);
        if(cvd->LincludeRhsMix)    // do this for SEMI-IMPLICIT, not STRANG
            NV_Ith_S(ydot,k) += cvd->tVarMap.at(k)->rhsMix.at(cvd->iC);
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
/** Message for bad cvode flag
 *  @param flag \input value flag
 *  @param func \input function with error
 */

void cvodeDriver::testCVflag(int flag, std::string func) {
    if(flag != CV_SUCCESS) {
        cout << endl << "ERROR in " << func << endl;
        exit(0);
    }
}