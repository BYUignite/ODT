/**
 * @file dv_temp_profile.cc
 * Header file for class dv_temp_profile
 */


#include "dv_temp_profile.h"
#include "domain.h"
#include <iostream>
#include <cstdlib>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/*! Constructor
 */

dv_temp_profile::dv_temp_profile(domain      *line,
                             const string s,
                             const bool   Lt,
                             const bool   Lo) :
    dv(line, s, Lt, Lo) {
    
    for(int i=0; i<domn->io->dvParams["xT"].size(); ++i) {
        x.push_back(domn->io->dvParams["xT"][i][0].as<double>());
        T.push_back(domn->io->dvParams["xT"][i][1].as<double>());
    }

    Linterp = Linear_interp(x, T);

    if(domn->pram->domainLength != x.back()){
        cout << endl << "ERROR: when using dv_temp_profile, domainLength must equal x profile" << endl;
        exit(0);
    }

}

////////////////////////////////////////////////////////////////////////////////
/*! Set temperature from the gas state
 *  @param ipt \input optional point to compute at
 */

void dv_temp_profile::setVar(const int ipt){

    d.resize(domn->ngrd);
    if(ipt == -1)
        for(int i=0; i<domn->ngrd; i++)
            d.at(i) = Linterp.interp(domn->pos->d.at(i));
    else
        d.at(ipt) = Linterp.interp(domn->pos->d.at(ipt));
}