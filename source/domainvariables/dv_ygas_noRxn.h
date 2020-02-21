/**
 * @file dv_ygas_noRxn.h
 * Header file for class dv_ygas_noRxn
 */

#pragma once

#include "dv_ygas.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_ygas_noRxn of parent dv_ygas object.
 *
 *  @author David O. Lignell
 */

class dv_ygas_noRxn : public dv_ygas {


    //////////////////// DATA MEMBERS //////////////////////

    private:

    public:

    //////////////////// MEMBER FUNCTIONS /////////////////

        virtual void getRhsSrc(const int ipt=-1);

    private:

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_ygas_noRxn(){}
        dv_ygas_noRxn(domain     *line,
                      const string s,
                      const bool   Lt,
                      const bool   Lo=true) : dv_ygas(line, s, Lt, Lo) {}

        virtual ~dv_ygas_noRxn(){}

};


////////////////////////////////////////////////////////////////////////////////

