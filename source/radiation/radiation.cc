/**
 * @file radiation.cc
 * @brief Source file for class radiation
 */

#include "radiation.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Static members


///////////////////////////////////////////////////////////////////////////////
/** Constructor for radiation class object
 *
 *  @author Victoria B. Stephens 
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
//        init_planck_mean_coefs();
        radProps = new rad_planck_mean();
    }
    else if (domn->pram->radCoefType == "WSGG") {
//        init_WSGG_coefs();
//        nGG = 5;                        // 1 clear gas and 4 gray gases
        radProps = new rad_wsgg();
    }
    else if (domn->pram->radCoefType == "RCSLW") {
//        init_RCSLW_coefs();
        radProps = new rad_rcslw(4);
    }
    else {
        *domn->io->ostrm << endl << "ERROR: radCoefType not recognized" << endl;
        exit(0);
    }

//    radProps = new radiationProperties();
//    radProps->init(domn);

    sigmaSB = 5.670E-8;                     // W/m2*K4

}