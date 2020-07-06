/**
 * @file dv_rho_mf.cc
 * @brief Source file for class dv_rho_mf
 */


#include "dv_rho_mf.h"
#include "domain.h"

////////////////////////////////////////////////////////////////////////////////
/*! dv_rho_mf  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_rho_mf::dv_rho_mf(domain    *line,
        const      string s,
        const bool Lt,
        const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;
    d             = vector<double>(domn->ngrd, domn->pram->rho0);

    rho0      = domn->io->streamProps["rho0"]     ?   domn->io->streamProps["rho0"].as<double>()  :  1.0;    //errMsg<double>("rho0");
    rho1      = domn->io->streamProps["rho1"]     ?   domn->io->streamProps["rho1"].as<double>()  :  1.0;    //errMsg<double>("rho1");

    temp0     = domn->io->streamProps["temp0"]     ?  domn->io->streamProps["temp0"].as<double>()  :  1.0;    //errMsg<double>("temp0");
    temp1     = domn->io->streamProps["temp1"]     ?  domn->io->streamProps["temp1"].as<double>()  :  1.0;    //errMsg<double>("temp1");
    tempFlame = domn->io->streamProps["tempFlame"] ?  domn->io->streamProps["tempFlame"].as<double>() : 1.0;  //errMsg<double>("tempFlame");
    Zst       = domn->io->streamProps["Zst"]       ?  domn->io->streamProps["Zst"].as<double>()    :  0.5;    //errMsg<double>("Zst");

    if ( !(domn->io->streamProps["tempFlame"]) && ( domn->io->streamProps["temp1"] || domn->io->streamProps["temp1"] ) ) {
        //    if ( ( ( temp0 > 1.01 ) || (temp1 > 1.01 ) ) && ( tempFlame < 1.01 ) ) {
        cout << endl << "ERROR: If tempFlame is not specified, do not specify different temp0 and temp1." << endl;
        cout << endl << "If you want different temp0 and temp1, specify a colinear tempFlame with Zst." << endl;
        exit(2);
    }
    if ( !( ( Zst >= 0.0 ) && (Zst <= 1.000001 ) ) ) {
        cout << endl << "ERROR: Need Zst between 0 and 1." << endl;
        exit(2);
    }


    }

    ////////////////////////////////////////////////////////////////////////////////
    /*! dv_rho_mf merger2cells function
     *
     * @param imrg    \input merge cells imrg and imrg+1
     * @param imrg \input merge cells imrg and imrg+1
     * @param m1   \input mass in cell imrg
     * @param m2   \input mass in cell imrg
     * @param LconstVolume \input (for posf, default is false)
     */

    void dv_rho_mf::merge2cells(const int    imrg,
            const double m1,
            const double m2,
            const bool   LconstVolume) {

        setVar(imrg);
        d.erase(d.begin() + imrg+1);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /*! dv_rho_mf setVar function
     *  @param ipt \input optional point to compute at
     */

    void dv_rho_mf::setVar(const int ipt) {

        double temp;
        d.resize(domn->ngrd,domn->pram->rho0);
        if(ipt == -1)
            for(int i=0; i<domn->ngrd; i++) {
                if ( domn->mixf->d.at(i) < Zst )
                    temp = temp0 * ( Zst - domn->mixf->d.at(i) ) / Zst + tempFlame * domn->mixf->d.at(i) / Zst;
                else
                    temp = temp1 * ( domn->mixf->d.at(i) - Zst ) / ( 1 - Zst ) + tempFlame * ( 1 - domn->mixf->d.at(i) ) / ( 1 - Zst );
                d.at(i) = 1.0 / temp / ( (1-domn->mixf->d.at(i))/(rho0*temp0) + domn->mixf->d.at(i)/(rho1*temp1) );
            }
        else {
            //d.at(ipt) = 1.0/( (1-domn->mixf->d.at(ipt))/rho0 + domn->mixf->d.at(ipt)/rho1 );
            if ( domn->mixf->d.at(ipt) < Zst )
                temp = temp0 * ( Zst - domn->mixf->d.at(ipt) ) / Zst + tempFlame * domn->mixf->d.at(ipt) / Zst;
            else
                temp = temp1 * ( domn->mixf->d.at(ipt) - Zst ) / ( 1 - Zst ) + tempFlame * ( 1 - domn->mixf->d.at(ipt) ) / ( 1 - Zst );
            d.at(ipt) = 1.0 / temp / ( (1-domn->mixf->d.at(ipt))/(rho0*temp0) + domn->mixf->d.at(ipt)/(rho1*temp1) );

        }
    }

