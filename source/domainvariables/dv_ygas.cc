/**
 * @file dv_ygas.cc
 * @brief Source file for class dv_ygas
 */


#include "dv_ygas.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////
// Setup static members

int dv_ygas::nspc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

//void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

////////////////////////////////////////////////////////////////////////////////
/*! dv_ygas  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_ygas::dv_ygas(domain  *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(domn->ngrd, 0.0);

    nspc = domn->gas->thermo()->nSpecies();

    string spName(var_name.begin()+2, var_name.end());        // var_name is like y_O2. Just want the O2 part.
    kMe                = domn->gas->thermo()->speciesIndex(spName);

    aP   = domn->io->params["aP"]   ? domn->io->params["aP"].as<double>()       : 1.0;
    //aP_x = domn->io->params["aP_x"] ? domn->io->params["aP_x"].as<double>()/4.6 : 1.0E-10;
    aP_x = domn->io->params["aP_x"] ? domn->io->params["aP_x"].as<double>() : 1.0E-10;

}

////////////////////////////////////////////////////////////////////////////////
/*! lv source term part of the rhs function
 *  @param ipt \input optional point to compute source at.
 * Gas temperature needs to be set to use problem specific RR
 */

void dv_ygas::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    static vector<vector<double> > rrSpc(nspc);    // [nspc][ngrd]
    static vector<double>          yi(nspc,0.0);       // [nspc]
    static vector<double>          rr(nspc, 0.0);       // [nspc]

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = domn->ngrd-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                 // to save cost, compute needed terms for all dv_ygas objects using this one.

        //-------- Soot source terms
        //         Soot sources need to be done first: check with L_source_done flag
        //         If not done, do them here.
        //         Then when soot sources are done, the L_source_done is true and we'll avoid double computing.
        //         This allows arbitrary ordering of the domain variables. If soot always comes before ygas in the list we wouldn't need this.

        if(domn->pram->Lsoot && !domn->svar[0]->L_source_done)
            for(int k=0; k<domn->pram->nsvar; k++)
                domn->svar[k]->getRhsSrc(ipt);            // this will set L_source_done flag true

        for(int k=0; k<nspc; k++) {        // resize and reset rrSpc for each gas species k
            rrSpc[k].resize(domn->ngrd);
            for (int j=0; j<domn->ngrd; j++)
                rrSpc[k][j] = 0;
        }

        for(int i=iS; i<=iE; i++) {

            // make sure rho and T are set first (though it should be for the diffuser at least).
            for(int k=0; k<nspc; k++)
                yi.at(k) = domn->ysp[k]->d.at(i);
            domn->chem->getProblemSpecificRR(domn->rho->d.at(i), domn->temp->d.at(i), domn->pram->pres, &yi.at(0), &rr.at(0));

            for(int k=0; k<nspc; k++) {
                rrSpc.at(k).at(i) = rr.at(k) * domn->gas->thermo()->molecularWeight(k) / domn->rho->d.at(i);   // kmol/(m³ s)*(kg/kmol)*(kg/m3) = 1/s
                if (domn->pram->Lsoot)
                    rrSpc.at(k).at(i) += gasSootSources.at(k).at(i);    // gasSootSources set in soot source terms
            }
        }

    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSpc.at(kMe).at(i) *
             (domn->mimx->time < aP_x ? aP : 1.0 ); //doldb this domain
             //(1.0+(aP-1.0)*exp(-domn->mimx->time/aP_x)); //doldb this domain

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);
}

////////////////////////////////////////////////////////////////////////////////
/*! lv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void dv_ygas::getRhsMix(const vector<double> &gf,
                        const vector<double> &dxc){

    if(!L_transported) return;

    rhsMix.resize(domn->ngrd, 0.0);

    setFlux(gf, dxc);

    //------------------ Compute the mixing term

    for(int i=0,ip=1; i<domn->ngrd; i++, ip++)
       rhsMix.at(i) = -domn->pram->cCoord / (domn->rho->d.at(i) * dxc.at(i)) *
                    (flux.at(ip) * pow(abs(domn->posf->d.at(ip)), domn->pram->cCoord-1) -
                     flux.at(i)  * pow(abs(domn->posf->d.at(i) ), domn->pram->cCoord-1));

    if(domn->pram->Lspatial)
        for(int i=0; i<domn->ngrd; i++)
            rhsMix.at(i) /= domn->uvel->d.at(i);

}

////////////////////////////////////////////////////////////////////////////////
/*! lv set face fluxes
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void dv_ygas::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    static vector<vector<double> > rhoD(nspc);
    static vector<vector<double> > rhoD_f(nspc);
    static vector<vector<double> > rhoDYinvM(nspc);
    static vector<vector<double> > rhoDYinvM_f(nspc);
    static vector<vector<double> > ysp_f(nspc);
    static vector<double>          Di(nspc);
    static vector<double>          dMdx;
    static vector<double>          MMw;

    if(kMe==0) {                 // to save cost, compute needed terms for all dv_ygas objects using this one.

        MMw.resize(domn->ngrd);
        dMdx.resize(domn->ngrdf);
        for(int k=0; k<nspc; k++) {
            rhoD.at(k).resize(domn->ngrd);
            rhoD_f.at(k).resize(domn->ngrdf);
            rhoDYinvM.at(k).resize(domn->ngrd);
            rhoDYinvM_f.at(k).resize(domn->ngrdf);
            ysp_f.at(k).resize(domn->ngrdf);
            domn->ysp[k]->flux.resize(domn->ngrdf);
        }
        for(int i=0; i<domn->ngrd; i++) {
            try {
                domn->domc->setGasStateAtPt(i);
            } catch (const CanteraError& e) {
                throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
            }
            MMw.at(i) = domn->gas->thermo()->meanMolecularWeight();
            domn->gas->transport()->getMixDiffCoeffs(&Di.at(0));
            for (int k=0; k<nspc; k++) {
                rhoD.at(k).at(i)      = domn->rho->d.at(i)*Di.at(k);
                rhoDYinvM.at(k).at(i) = rhoD.at(k).at(i)*domn->ysp[k]->d.at(i)/MMw.at(i);
            }
        }
        //for (int k=0; k<nspc; k++) {      // this breaks sometimes due to division issues, use linear instead (below)
        //    interpVarToFacesHarmonic(rhoD.at(k),      rhoD_f.at(k));
        //    interpVarToFacesHarmonic(rhoDYinvM.at(k), rhoDYinvM_f.at(k));
        //}

        for (int k=0; k<nspc; k++) {         // linear interpolation.
            for(int i=0; i<domn->ngrdf; i++){
                rhoD_f.at(k).at(i)      = linearInterpToFace(i, rhoD.at(k));
                rhoDYinvM_f.at(k).at(i) = linearInterpToFace(i, rhoDYinvM.at(k));
                ysp_f.at(k).at(i)       = linearInterpToFace(i, domn->ysp[k]->d);
            }
        }

        dMdx.at(0)          = 0.0;
        dMdx.at(domn->ngrd) = 0.0;
        for (int i=1, im=0; i < domn->ngrd; i++, im++)
            dMdx.at(i) = gf.at(i) * (MMw.at(i) - MMw.at(im));

        //==========================

        //---------- Interior faces

        double jstar;             // correction flux so that all fluxes sum to zero. This is equal to using a correction velocity
                                  // j_i_corrected = j_i - Yi*jstar; jstar = sum(j_i).
                                  // the previous approch just made j_last = -sum(j_all_but_last), where N2 was last.
        for (int i=1, im=0; i < domn->ngrd; i++, im++) {
            jstar = 0.0;
            for(int k=0; k<nspc; k++) {
                domn->ysp[k]->flux.at(i) = -rhoD_f.at(k).at(i)*gf.at(i)*(domn->ysp[k]->d.at(i) - domn->ysp[k]->d.at(im))
                    -rhoDYinvM_f.at(k).at(i)*dMdx.at(i);
                jstar += domn->ysp[k]->flux.at(i);
            }
            for(int k=0; k<nspc; k++)
                domn->ysp[k]->flux.at(i) -= ysp_f.at(k).at(i) * jstar;
        }

        //---------- Boundary faces

        for(int k=0; k<nspc; k++) {
            domn->ysp[k]->flux.at(0)          = 0.0;                          // for wall or outflow; todo: wall flame specifics
            domn->ysp[k]->flux.at(domn->ngrd) = 0.0;
        }
    }
}