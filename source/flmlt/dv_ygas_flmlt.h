/**
 * @file dv_ygas_flmlt.h
 * Header file for class dv_ygas_flmlt
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_ygas_flmlt of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_ygas_flmlt : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int                     nspc;        ///< number of gas species

        int                            kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_ygas_flmlt(){}
        dv_ygas_flmlt(domain      *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~dv_ygas_flmlt(){}

};


////////////////////////////////////////////////////////////////////////////////

