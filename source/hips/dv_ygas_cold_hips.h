/**
 * @file dv_ygas_cold_hips.h
 * Header file for class dv_ygas_cold_hips
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_ygas_hips of parent dv object.
 *
 *  @author David O. Lignell
 */

class dv_ygas_cold_hips : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int   nspc;        ///< number of gas species
        int          kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

        double Da1;
        double Da2;

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_ygas_cold_hips(){}
        dv_ygas_cold_hips(domain   *line,
                const string       s,
                const bool         Lt,
                const bool         Lo=true);

        virtual ~dv_ygas_cold_hips(){}

};


////////////////////////////////////////////////////////////////////////////////

