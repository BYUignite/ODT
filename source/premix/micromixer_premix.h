/**
 * @file micromixer_premix.h
 * Header file for class micromixer_premix
 */

#pragma once

#include "micromixer.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer_premix object
 *
 *  @author David O. Lignell
 */

class micromixer_premix : public micromixer {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        double tNextAdapt;


    //////////////////// MEMBER FUNCTIONS /////////////////

    protected:

        virtual bool adaptGridIfNeeded();
        virtual void updateGrid();           ///< enforce the continuity condition: (e.g., rho*dx = const).

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer_premix() : micromixer() { tNextAdapt = -1.0; }
        virtual ~micromixer_premix(){ }

};


////////////////////////////////////////////////////////////////////////////////


