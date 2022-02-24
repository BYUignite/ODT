/**
 * @file radiation.cc
 * @brief Source file for class radiation
 */

#include "radiation.h"
#include "domain.h"
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////
// Static members


///////////////////////////////////////////////////////////////////////////////
/** Constructor for radiation class object
 *
 *  @param p_domn \input pointer to domain object
 *
 */

radiation::radiation(domain *p_domn) {

    domn = p_domn;

    if(!domn->pram->Lrad)
        return;

    // --------------- populate list of gas species indices

    bool fmissing = false;

    nRadSp = 4;                     // CH4, CO2, H2O, CO
    iRadIndx.resize(nRadSp);

    int isp;

    isp = domn->gas->thermo()->speciesIndex("CH4");
    isp = (isp > 0) ? isp : domn->gas->thermo()->speciesIndex("ch4");
    iRadIndx[0] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->thermo()->speciesIndex("CO2");
    isp = (isp > 0) ? isp : domn->gas->thermo()->speciesIndex("co2");
    iRadIndx[1] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->thermo()->speciesIndex("H2O");
    isp = (isp > 0) ? isp : domn->gas->thermo()->speciesIndex("h2o");
    iRadIndx[2] = isp;
    if(isp < 0) fmissing = true;

    isp = domn->gas->thermo()->speciesIndex("CO");
    isp = (isp > 0) ? isp : domn->gas->thermo()->speciesIndex("co");
    iRadIndx[3] = isp;
    if(isp < 0) fmissing = true;

    if(fmissing)
        cout << endl << "Warning: one or more radiating species missing from gas mechanism" << endl;

    // --------------- create radProps object

    if (domn->pram->radCoefType == "PLANCKMEAN") {
        radProps = new rad_planck_mean();
    }
    else if (domn->pram->radCoefType == "WSGG") {
        radProps = new rad_wsgg();
    }
    else if (domn->pram->radCoefType == "RCSLW") {
        radProps = new rad_rcslw(4, 1800, 101325, 0.0, 0.189416, 0.0941293, 0.000859649);       // init w/ stoich methane/air at 1800K
    }
    else {
        *domn->io->ostrm << endl << "ERROR: radCoefType not recognized" << endl;
        exit(0);
    }

}