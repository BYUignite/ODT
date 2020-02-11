/**
 * @file micromixer_flmlt.h
 * Header file for class micromixer_flmlt
 */

#pragma once

#include "micromixer.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer_flmlt object
 *
 *  @author David O. Lignell
 */

class micromixer_flmlt : public micromixer {

    //////////////////// DATA MEMBERS //////////////////////

    private:

        double tNextAdapt;


    //////////////////// MEMBER FUNCTIONS /////////////////

    protected:

        virtual void set_oldrho_or_rhov() {return;}
        virtual void updateGrid()         {return;}
        virtual void setNominalStepSize();   ///< sets a nominal dt for the whole period
        virtual void setGf();
        virtual void setGridDxcDx();

        virtual bool adaptGridIfNeeded();

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer_flmlt() : micromixer() { tNextAdapt = -1.0; }
        virtual ~micromixer_flmlt(){ }

};


////////////////////////////////////////////////////////////////////////////////


