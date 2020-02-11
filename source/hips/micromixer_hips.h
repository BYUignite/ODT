/**
 * @file micromixer_hips.h
 * Header file for class micromixer_hips
 */

#pragma once

#include "micromixer.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing micromixer_hips object
 *
 *  @author David O. Lignell
 */

class micromixer_hips : public micromixer {

    //////////////////// DATA MEMBERS //////////////////////

    private:


    //////////////////// MEMBER FUNCTIONS /////////////////

    private:

    void mixAcrossLevelTree(const int kVar, const int iMixLevel, const int iTree); //, const int iTree=-1);
    void forceProfile();

    public:

    virtual void advanceOdt(const double p_tstart, const double p_tend, const int iLevel = -1);    // iLevel is for hips

    protected:

        virtual void setGf()              {return;}
        virtual void setGridDxcDx()       {return;}
        virtual void set_oldrho_or_rhov() {return;}
        virtual void setNominalStepSize();     ///< sets a nominal dt for the whole period

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        micromixer_hips() : micromixer() { }
        virtual ~micromixer_hips(){ }

};


////////////////////////////////////////////////////////////////////////////////


