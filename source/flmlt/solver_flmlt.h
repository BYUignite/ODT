/**
 * @file solver_flmlt.h
 * Header file for class solver_flmlt
 */

#pragma once

#include "solver.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver_flmlt object
 *
 *  @author David O. Lignell
 */

class solver_flmlt : public solver {


    //////////////////// DATA MEMBERS //////////////////////

    public:


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver_flmlt(){}
        virtual void init(domain *p_domn);
        virtual ~solver_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////


