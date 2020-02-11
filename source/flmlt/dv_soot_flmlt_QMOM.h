/**
 * @file dv_soot_flmlt_QMOM.h
 * Header file for class dv_soot_flmlt_QMOM
 *
 * https://en.wikipedia.org/wiki/Dominance_(C%2B%2B)
 */

#pragma once

#include "dv_soot_QMOM.h"
#include "dv_soot_flmlt.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_QMOM of parent dv object.
 *
 *  @author Victoria B. Lansinger
 *
 *  NOTE: getRhsMix will come from dv_soot_flmlt, and getRhsSrc will come from dv_soot_QMOM
 */

class dv_soot_flmlt_QMOM : public dv_soot_flmlt, public dv_soot_QMOM {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_flmlt_QMOM(domain       *line,
                           const string  s,
                           const bool    Lt,
                           const bool    Lo=true) :
                 dv_soot(line, s, Lt, Lo),
                 dv_soot_flmlt(line, s, Lt, Lo),
                 dv_soot_QMOM(line, s, Lt, Lo) {}

        virtual ~dv_soot_flmlt_QMOM(){}

};


////////////////////////////////////////////////////////////////////////////////



