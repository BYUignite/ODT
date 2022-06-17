/**
 * @file solver_premix.h
 * Header file for class solver_premix
 */

#pragma once

#include "solver.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver_premix object
 *
 *  @author David O. Lignell
 */

class solver_premix : public solver {


    //////////////////// DATA MEMBERS //////////////////////

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

        virtual void calculateSolution();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver_premix(domain *p_domn);
        virtual ~solver_premix(){}

};


////////////////////////////////////////////////////////////////////////////////


