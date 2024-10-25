/**
 * @file dv_posf.h
 * @brief Header file for class dv_posf
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_posf of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_posf : public dv {

    public:

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

        virtual void splitCell(const int            isplt,
                               const int            nsplt,
                               const vector<double> &cellFaces);

        virtual void setDvFromRegion(const int i1, const int i2);
        virtual void resize();

    private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_posf(){}
        dv_posf(domain      *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_posf(){}

};


////////////////////////////////////////////////////////////////////////////////

