/**
 * @file dv.h
 * @brief Header file for class \ref dv
 */

#pragma once

#include <string>
#include <vector>
#include "chemMech.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing base dv object.
 *  Derived dv will be u,v,w,Yi,etc.
 *
 *  @author David O. Lignell
 */

class dv {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        string                        var_name;               ///< name of variable
        vector<double>                d;                      ///< the data
        bool                          L_transported;          ///< flag true if var is transported
        bool                          L_output;               ///< flag true if included in output
        bool                          LagSrc;                 ///< flag to lag source term in implicit solve (initially put in for enthalpy radiation)

        domain                        *domn;                  ///< pointer to domain object (parent)

        vector<double>                rhsSrc;                 ///< the data
        vector<double>                rhsMix;                 ///< the data

        vector<double>                flux;

        //---------- for gas-soot coupling

        bool                           L_source_done;         ///< flag set true when source updated: (for gas-soot sources)
        static vector<vector<double> > gasSootSources;        ///< [nspc][ngrd] source terms for gas from soot rxns

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void   setVar(const int ipt=-1){}

        virtual void   merge2cells(const int    imrg,
                                   const double m2,
                                   const double m1,
                                   const bool   LconstVolume=false);

        virtual void   splitCell(const int            isplt,
                                 const int            nsplt,
                                 const vector<double> &cellFaces);

        virtual void   getRhsSrc(const int ipt=-1){if(!L_transported) return;}
        virtual void   getRhsMix(const vector<double> &gf,
                                 const vector<double> &dxc){if(!L_transported) return;}

        virtual void   interpVarToFacesHarmonic(const vector<double> &cvar, vector<double> &fvar);
        virtual double linearInterpToFace(const int &iface, const vector<double> &vec);
        virtual void   setDvFromRegion(const int i1, const int i2);
        virtual void   resize();

        void           resetSourceFlags();


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv(){}
        dv(domain       *line,
           const string s,
           const bool   Lt,
           const bool   Lo=true);

        virtual ~dv(){}

};


////////////////////////////////////////////////////////////////////////////////