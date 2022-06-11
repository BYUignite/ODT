/**
 * @file dv_soot.cc
 * Source file for class dv_soot
 *
 */

#include "dv_soot.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

using namespace soot;

////////////////////////////////////////////////////////////////////////////////
// Static members

vector<double>         *dv_soot::yi;
int                     dv_soot::nsvar;
int                     dv_soot::N;

int                     dv_soot::i_c2h2;
int                     dv_soot::i_o2;
int                     dv_soot::i_o;
int                     dv_soot::i_h;
int                     dv_soot::i_h2;
int                     dv_soot::i_oh;
int                     dv_soot::i_h2o;
int                     dv_soot::i_co;
int                     dv_soot::i_c10h8;
int                     dv_soot::i_c12h8;
int                     dv_soot::i_c12h10;
int                     dv_soot::i_c14h10;
int                     dv_soot::i_c16h10;
int                     dv_soot::i_c18h10;
//int                     dv_soot::i_elem_c;
//int                     dv_soot::i_elem_h;
//vector<int>             dv_soot::i_pah;
//vector<double>          dv_soot::MW_sp;

////////////////////////////////////////////////////////////////////////////////
/*! dv_soot  constructor function
 *
 * @param
 * @param
 */

dv_soot::dv_soot(domain     *line,
                 const      string s,
                 const bool Lt,
                 const bool Lo) {

    domn               = line;
    var_name           = s;
    L_transported      = Lt;
    L_output           = Lo;
    d                  = vector<double>(domn->ngrd, 0.0);

    kMe         = N++;                    ///< kMe

    if(L_transported){
        rhsMix.resize(domn->ngrd, 0.0);
        rhsSrc.resize(domn->ngrd, 0.0);
    }

    //----------------------------------------

    if (kMe == 0) {

        gasSootSources.resize(domn->gas->thermo()->nSpecies());                     // set first dim of [nspc][ngrd]

        //-------------- populate user-specified parameters

        Cmin             = domn->pram->Cmin;
        rhoSoot          = domn->pram->rhoSoot;
        b_coag           = domn->pram->b_coag;
        nucleation_mech  = soot::str2nucMech(domn->pram->nucleation_mech);      // Nucleation: NONE, LL, LIN, PAH
        growth_mech      = soot::str2grwMech(domn->pram->growth_mech);          // Surface growth: NONE, LL, LIN, HACA
        oxidation_mech   = soot::str2oxiMech(domn->pram->oxidation_mech);       // Oxidation: NONE, LL, LEE_NEOH, NSC_NEOH, HACA
        coagulation_mech = soot::str2coaMech(domn->pram->coagulation_mech);     // Coagulation: NONE, LL, FUCHS, FRENK
        psd_mech         = soot::str2psdMech(domn->pram->PSD_method);           // PSD mechanisms: MONO, LOGN, QMOM, MOMIC
        nsvar            = domn->pram->nsvar;                             // number of soot moments

        //-------------- create soot model and soot state objects

        SM = new sootModel(psd_mech, nsvar, nucleation_mech, growth_mech, oxidation_mech, coagulation_mech);
        S = new state();

        //-------------- populate list of gas and PAH species indices

        int isp;

        isp = domn->gas->thermo()->speciesIndex("C2H2");                                      // Copy/paste to add species to this list
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c2h2");
        i_c2h2 = isp;

        isp = domn->gas->thermo()->speciesIndex("O2");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("o2");
        i_o2 = isp;

        isp = domn->gas->thermo()->speciesIndex("O");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("o");
        i_o = isp;

        isp = domn->gas->thermo()->speciesIndex("H");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("h");
        i_h = isp;

        isp = domn->gas->thermo()->speciesIndex("H2");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("h2");
        i_h2 = isp;

        isp = domn->gas->thermo()->speciesIndex("OH");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("oh");
        i_oh = isp;

        isp = domn->gas->thermo()->speciesIndex("H2O");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("h2o");
        i_h2o = isp;

        isp = domn->gas->thermo()->speciesIndex("CO");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("co");
        i_co = isp;

//        isp = domn->gas->thermo()->elementIndex("C");
//        isp = (isp >= 0) ? isp : domn->gas->thermo()->elementIndex("c");
//        i_elem_c = isp;
//
//        isp = domn->gas->thermo()->elementIndex("H");
//        isp = (isp >= 0) ? isp : domn->gas->thermo()->elementIndex("h");
//        i_elem_h = isp;

        isp = domn->gas->thermo()->speciesIndex("C10H8");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c10h8");
        i_c10h8 = isp;

        isp = domn->gas->thermo()->speciesIndex("C12H8");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c12h8");
        i_c12h8 = isp;

        isp = domn->gas->thermo()->speciesIndex("C12H10");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c12h10");
        i_c12h10 = isp;

        isp = domn->gas->thermo()->speciesIndex("C14H10");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c14h10");
        i_c14h10 = isp;

        isp = domn->gas->thermo()->speciesIndex("C16H10");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c16h10");
        i_c16h10 = isp;

        isp = domn->gas->thermo()->speciesIndex("C18H10");
        isp = (isp >= 0) ? isp : domn->gas->thermo()->speciesIndex("c18h10");
        i_c18h10 = isp;

//        for(int i=0; i<domn->pram->PAH_species.size(); i++) {
//            i_pah.push_back( domn->gas->thermo()->speciesIndex(domn->pram->PAH_species[i]) );
//            if (i_pah[i] < 0) {
//                cout << endl << "ERROR: Invalid PAH species provided: check input file and mechanism." << endl;
//                exit(0);
//            }
//        }

        //-------------- TO DO: test that the species present are sufficient for the desired soot mechanism
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! dv mixing term part of the rhs function
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 * Note: temperature variable should be set already.
 * Note: species mass fluxes should be set already.
 */

void dv_soot::getRhsMix(const vector<double> &gf,
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
/*! dv set face fluxes
 * @param gf  \input grid geometric factor for derivatives
 * @param dxc \input = \abs{\Delta(x^c)}
 */

void dv_soot::setFlux(const vector<double> &gf,
                      const vector<double> &dxc){

    int nVar = domn->io->sootParams["nsvar"].as<int>();

    static vector<vector<double> > dvisc(nVar);
    static vector<vector<double> > dvisc_f(nVar);

    if(kMe==0) {                 // to save cost, compute needed terms for all dv_soot objects using this one.

        for(int k=0; k<nVar; k++) {
            dvisc.at(k).resize(domn->ngrd);
            dvisc_f.at(k).resize(domn->ngrdf);
            domn->svar[k]->flux.resize(domn->ngrdf);
        }
        for(int i=0; i<domn->ngrd; i++) {
            try {
                domn->domc->setGasStateAtPt(i);   // keep this commented for decoupled soot
            } catch (const odtCanteraError& e) {
                throw odtCanteraError(STR_TRACE, "setGasStateAtPt",e);
            }

            for (int k=0; k<nVar; k++) {
                dvisc.at(k).at(i) = domn->dvisc->d.at(i);
            }
        }

        //for (int k=0; k<nVar; k++) {      // this breaks sometimes due to division issues, use linear instead (below)
        //    interpVarToFacesHarmonic(dvisc.at(k),      dvisc_f.at(k));
        //}

        for (int k=0; k<nVar; k++) {         // linear interpolation.
            for(int i=0; i<domn->ngrdf; i++){
                dvisc_f.at(k).at(i)     = linearInterpToFace(i, dvisc.at(k));
            }
        }

        //==========================

        //---------- Interior faces

        for (int i=1, im=0; i < domn->ngrd; i++, im++) {

            //-------------- Thermophoretic transport is "hyperbolic", flux looks like f=v*phi, so upwind it for stability

            for(int k=0; k<nVar; k++){
                double vel = -0.55415 * dvisc_f.at(k).at(i)*(log(domn->temp->d.at(i))-log(domn->temp->d.at(im)))*gf.at(i);
                double ii = (vel > 0) ? im : i;
                domn->svar[k]->flux.at(i) = vel * domn->svar[k]->d.at(ii);
            }
        }

        //---------- Boundary faces

        for(int k=0; k<nVar; k++) {
            domn->svar[k]->flux.at(0)          = 0.0;                          // for wall or outflow; todo: wall flame specifics
            domn->svar[k]->flux.at(domn->ngrd) = 0.0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! dv source term part of the rhs function
 *
 *  @param ipt \input optional point to compute source at.
 *
 *  adapted from dv_ygas.cc
 */

void dv_soot::getRhsSrc(const int ipt) {

    if(!L_transported)
        return;

    if(L_source_done)                                    //  return if this routine was already called (by dv_ygas)
        return;

    rhsSrc.resize(domn->ngrd, 0.0);

    for(int k=0; k<gasSootSources.size(); k++)           // resize gasSootSources for each gas species k
        gasSootSources[k].resize(domn->ngrd,0.0);

    static vector<vector<double> > rrSvar(nsvar);                  // temp storage for moment rates [nsvar][ngrd]
    static vector<double>          yGas(8);                     // number of gas species involved with soot = 8
    static vector<double>          yPAH(6);                     // number of PAH species involved with soot = 6
    static vector<double>          ySootVar(domn->pram->nsvar);    // number of soot moments

    int iS, iE;
    if(ipt==-1) {
        iS = 0;
        iE = (domn->ngrd)-1;
    }
    else {
        iS = ipt;
        iE = ipt;
    }

    if(kMe==0) {                                       // to save cost, compute needed terms for all dv_soot objects using this one.

        for(int k=0; k<nsvar; k++)                     // for each mom/section,
            rrSvar.at(k).resize(domn->ngrd);           // resize to match grid size

        for(int k=0; k<gasSootSources.size(); k++)         // resize gasSootSources for each gas species k
            gasSootSources[k].resize(domn->ngrd,0.0);

        for(int i=iS; i<=iE; i++) {                    // loop over grid points

            yGas[0] = (i_h < 0)    ? 0 : domn->ysp[i_h]->d[i];          // gas species mass fractions [H, H2, O, O2, OH, H2O, CO, C2H2]
            yGas[1] = (i_h2 < 0)   ? 0 : domn->ysp[i_h2]->d[i];
            yGas[2] = (i_o < 0)    ? 0 : domn->ysp[i_o]->d[i];
            yGas[3] = (i_o2 < 0)   ? 0 : domn->ysp[i_o2]->d[i];
            yGas[4] = (i_oh < 0)   ? 0 : domn->ysp[i_oh]->d[i];
            yGas[5] = (i_h2o < 0)  ? 0 : domn->ysp[i_h2o]->d[i];
            yGas[6] = (i_co < 0)   ? 0 : domn->ysp[i_co]->d[i];
            yGas[7] = (i_c2h2 < 0) ? 0 : domn->ysp[i_c2h2]->d[i];

            yPAH[0] = (i_c10h8 < 0)  ? 0 : domn->ysp[i_c10h8]->d[i];    // PAH species mass fractions [C10H8, C12H8, C12H10, C14H10, C16H10, C18H10]
            yPAH[1] = (i_c12h8 < 0)  ? 0 : domn->ysp[i_c12h8]->d[i];
            yPAH[2] = (i_c12h10 < 0) ? 0 : domn->ysp[i_c12h10]->d[i];
            yPAH[3] = (i_c14h10 < 0) ? 0 : domn->ysp[i_c14h10]->d[i];
            yPAH[4] = (i_c16h10 < 0) ? 0 : domn->ysp[i_c16h10]->d[i];
            yPAH[5] = (i_c18h10 < 0) ? 0 : domn->ysp[i_c18h10]->d[i];

            domn->domc->enforceSootMom();
            for(int k=0; k<ySootVar.size(); k++) 
                ySootVar[k] = domn->svar[k]->d[i];  // soot moment values [M0, M1, M2, M3]

            // set the thermodynamic state
            S->setState(domn->temp->d[i], domn->pram->pres, domn->rho->d[i], domn->dvisc->d[i],
                       domn->gas->thermo()->meanMolecularWeight(), yGas, yPAH, ySootVar, ySootVar.size());

            // calculate source terms
            SM->calcSourceTerms(*S);

            // retrieve soot moment values
            for (int k=0; k<nsvar; k++) {
                rrSvar[k][i] = SM->sourceTerms->sootSourceTerms.at(k);
            }

            // retrieve gas source term values
            if (i_h > 0)    gasSootSources[i_h][i]    = SM->sourceTerms->gasSourceTerms.at(gasSp::H);
            if (i_h2 > 0)   gasSootSources[i_h2][i]   = SM->sourceTerms->gasSourceTerms.at(gasSp::H2);
            if (i_o > 0)    gasSootSources[i_o][i]    = SM->sourceTerms->gasSourceTerms.at(gasSp::O);
            if (i_o2 > 0)   gasSootSources[i_o2][i]   = SM->sourceTerms->gasSourceTerms.at(gasSp::O2);
            if (i_oh > 0)   gasSootSources[i_oh][i]   = SM->sourceTerms->gasSourceTerms.at(gasSp::OH);
            if (i_h2o > 0)  gasSootSources[i_h2o][i]  = SM->sourceTerms->gasSourceTerms.at(gasSp::H2O);
            if (i_co > 0)   gasSootSources[i_co][i]   = SM->sourceTerms->gasSourceTerms.at(gasSp::CO);
            if (i_c2h2 > 0) gasSootSources[i_c2h2][i] = SM->sourceTerms->gasSourceTerms.at(gasSp::C2H2);
        }
    }

    for(int i=iS; i<=iE; i++)
        rhsSrc.at(i) = rrSvar.at(kMe).at(i);

    if(domn->pram->Lspatial)
        for(int i=iS; i<=iE; i++)
            rhsSrc.at(i) /= domn->uvel->d.at(i);

    L_source_done = true;

}