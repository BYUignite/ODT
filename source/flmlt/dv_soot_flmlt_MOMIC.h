/**
 * @file dv_soot_flmlt_MOMIC.h
 * Header file for class dv_soot_flmlt_MOMIC
 *
 * https://en.wikipedia.org/wiki/Dominance_(C%2B%2B)
 */

#pragma once

#include "dv_soot_MOMIC.h"
#include "dv_soot_flmlt.h"
#include <string>
#include <vector>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing child dv_soot_MOMIC of parent dv object.
 *
 *  @author Victoria B. Lansinger
 *
 *  NOTE: getRhsMix will come from dv_soot_flmlt, and getRhsSrc will come from dv_soot_MOMIC
 */

class dv_soot_flmlt_MOMIC : public dv_soot_flmlt, public dv_soot_MOMIC {

    //////////////////// DATA MEMBERS //////////////////////


    //////////////////// MEMBER FUNCTIONS /////////////////


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        dv_soot_flmlt_MOMIC(domain       *line,
                           const string  s,
                           const bool    Lt,
                           const bool    Lo=true) :
                 dv_soot(line, s, Lt, Lo),
                 dv_soot_flmlt(line, s, Lt, Lo),
                 dv_soot_MOMIC(line, s, Lt, Lo) {}

        virtual ~dv_soot_flmlt_MOMIC(){}

};


////////////////////////////////////////////////////////////////////////////////



