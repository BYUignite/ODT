/**
 * @file dv_ygas.h
 * @brief Header file for class dv_ygas
 */

#pragma once

#include "dv.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_ygas of parent lv object.
 *
 *  @author David O. Lignell
 */

class dv_ygas : public dv {


    //////////////////// DATA MEMBERS //////////////////////

    private:

        static int                     nspc;        ///< number of gas species

        int                            kMe;         ///< index of this spc in list: 0 to nspc-1; set from var_name

        double                         aP;
        double                         aP_x;

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////


        virtual void getRhsSrc(const int ipt=-1);
        virtual void getRhsMix(const vector<double> &gf,
                               const vector<double> &dxc);

    private:

        void setFlux(const vector<double> &gf,
                     const vector<double> &dxc);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_ygas(){}
        dv_ygas(domain     *line,
                const string s,
                const bool   Lt,
                const bool   Lo=true);

        virtual ~dv_ygas(){}

};


////////////////////////////////////////////////////////////////////////////////

