/**
 * @file dv_soot.cc
 * Header file for class dv_soot
 *
 * @author Victoria B. Lansinger
 */

#include "dv_soot.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// Static members

int                     dv_soot::N = 0;
int                     dv_soot::nsvar;

double                  dv_soot::Cmin;
double                  dv_soot::MW_c;
double                  dv_soot::MW_h;
double                  dv_soot::rhoSoot;
double                  dv_soot::b_coag;
string                  dv_soot::nucleation_mech;
string                  dv_soot::growth_mech;
string                  dv_soot::oxidation_mech;
string                  dv_soot::coagulation_mech;

double                  dv_soot::DIMER;
double                  dv_soot::m_dimer;

double                  dv_soot::T;
double                  dv_soot::P;
double                  dv_soot::rho;
double                  dv_soot::MW;
double                  dv_soot::mu;
vector<double>         *dv_soot::yi;

int                     dv_soot::i_c2h2;
int                     dv_soot::i_o2;
int                     dv_soot::i_h;
int                     dv_soot::i_h2;
int                     dv_soot::i_oh;
int                     dv_soot::i_h2o;
int                     dv_soot::i_co;
int                     dv_soot::i_elem_c;
int                     dv_soot::i_elem_h;
vector<int>             dv_soot::i_pah;
vector<double>          dv_soot::MW_sp;

double                  dv_soot::rC2H2_rSoot_n = 0.0;
double                  dv_soot::rH2_rSoot_ncnd = 0.0;
vector<double>          dv_soot::rPAH_rSoot_ncnd;
double                  dv_soot::rO2_rSoot_go = 0.0;
double                  dv_soot::rOH_rSoot_go = 0.0;
double                  dv_soot::rH_rSoot_go = 0.0;
double                  dv_soot::rCO_rSoot_go = 0.0;
double                  dv_soot::rH2_rSoot_go = 0.0;
double                  dv_soot::rC2H2_rSoot_go = 0.0;


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

    //----------------------------------------

    if (kMe == 0) {

        nsvar = domn->pram->nsvar;

        gasSootSources.resize(domn->gas->nSpecies());                          // set first dim of [nspc][ngrd]

        //-------------- populate list of gas species indices

        int isp;

        isp = domn->gas->speciesIndex("C2H2");                                      // Copy/paste to add species to this list
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("c2h2");
        i_c2h2 = isp;

        isp = domn->gas->speciesIndex("O2");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("o2");
        i_o2 = isp;

        isp = domn->gas->speciesIndex("H");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("h");
        i_h = isp;

        isp = domn->gas->speciesIndex("H2");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("h2");
        i_h2 = isp;

        isp = domn->gas->speciesIndex("OH");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("oh");
        i_oh = isp;

        isp = domn->gas->speciesIndex("H2O");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("h2o");
        i_h2o = isp;

        isp = domn->gas->speciesIndex("CO");
        isp = (isp >= 0) ? isp : domn->gas->speciesIndex("co");
        i_co = isp;

        isp = domn->gas->elementIndex("C");
        isp = (isp >= 0) ? isp : domn->gas->elementIndex("c");
        i_elem_c = isp;

        isp = domn->gas->elementIndex("H");
        isp = (isp >= 0) ? isp : domn->gas->elementIndex("h");
        i_elem_h = isp;

        for(int i=0; i<domn->pram->PAH_species.size(); i++) {
            i_pah.push_back( domn->gas->speciesIndex(domn->pram->PAH_species[i]) );
            if (i_pah[i] < 0) {
                cout << endl << "ERROR: Invalid PAH species provided: check input file and mechanism." << endl;
                exit(0);
            }
        }
        rPAH_rSoot_ncnd.resize(i_pah.size(),0.0);

        //--------------

        MW_sp = domn->gas->molecularWeights();
        MW_c  = domn->gas->atomicWeight(i_elem_c);
        MW_h  = domn->gas->atomicWeight(i_elem_h);

        //--------------

        Cmin             = domn->pram->Cmin;
        rhoSoot          = domn->pram->rhoSoot;
        b_coag           = domn->pram->b_coag;
        nucleation_mech  = domn->pram->nucleation_mech;
        growth_mech      = domn->pram->growth_mech;
        oxidation_mech   = domn->pram->oxidation_mech;
        coagulation_mech = domn->pram->coagulation_mech;

        //-------------- TO DO: test that the species present are sufficient for the desired soot mechanism
    }

    //----------------------------------------


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
            // domn->domc->setGasStateAtPt(i);   // keep this commented for decoupled soot
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
/*! set_gas_state function
 *
 *      Sets gas state properties. These will be implied by functions in this
 *         class that depend on the gas state.
 *
 *      @param T_p      /input (K)
 *      @param P_p      /input (Pa)
 *      @param rho_p    /input (gas density, kg/m3)
 *      @param MW_p     /input (mean molecular weight, kg/kmol)
 *      @param mu_p     /input (dynamic viscosity kg/m*s)
 *      @param y_p      /input mass fractions
 */

void dv_soot::set_gas_state_vars(const double &T_p,
                                 const double &P_p,
                                 const double &rho_p,
                                 const double &MW_p,
                                 const double &mu_p,
                                 vector<double> &y_p){
    T   = T_p;
    P   = P_p;
    rho = rho_p;
    MW  = MW_p;
    mu  = mu_p;
    yi  = &y_p;
}

////////////////////////////////////////////////////////////////////////////////
/*! getNucleationRate function
 *
 *      Calls appropriate function for nucleation chemistry
 *      based on the nucleation_mech flag.
 *      Returns soot nucleation rate in #/m3*s.
 *
 *      @param mi    /input vector of soot particle sizes    (optional, has default)
 *      @param wi    /input vector of soot particle weights  (optional, has default)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::getNucleationRate(const vector<double> &mi, const vector<double> &wi) {

    if (nucleation_mech      == "NONE")
        return 0;
    else if (nucleation_mech == "LL")
        return nucleation_LL();
    else if (nucleation_mech == "LIN")
        return nucleation_Linstedt();
    else if (nucleation_mech == "PAH")
        return nucleation_PAH(mi, wi);
    else {
        cout << endl << "ERROR: Invalid soot nucleation mechanism." << endl;
        exit(0);
    }

    return 0;

}

////////////////////////////////////////////////////////////////////////////////
/*! getGrowthRate function
 *
 *      Calls appropriate function for surface growth chemistry
 *      based on the growth_mech flag.
 *      Returns chemical surface growth rate in kg/m2*s.
 *
 *      @param M0     / input moment 0 (#/m3),       (optional, has a default)
 *      @param M1     / input moment 1 (kg-soot/m3), (optional, has a default)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::getGrowthRate(const double &M0, const double &M1) {

    if (growth_mech      == "NONE")
        return 0;
    else if (growth_mech == "LIN")
        return growth_Lindstedt();
    else if (growth_mech == "LL")
        return growth_LL(M0, M1);
    else if (growth_mech == "HACA")
        return growth_HACA(M0, M1);
    else {
        cout << endl << "ERROR: Invalid soot growth mechanism." << endl;
        exit(0);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/*! getOxidationRate function
 *
 *      Calls appropriate function for oxidation chemistry
 *      based on the oxidation_mech flag.
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      @param M0     / input moment 0 (#/m3),       (optional, has a default)
 *      @param M1     / input moment 1 (kg-soot/m3), (optional, has a default)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::getOxidationRate(const double &M0, const double &M1) {

    if (domn->pram->oxidation_mech      == "NONE")
        return 0;
    else if (domn->pram->oxidation_mech == "LL")
        return oxidation_LL();
    else if (domn->pram->oxidation_mech == "LEE_NEOH")
        return oxidation_Lee_Neoh();
    else if (domn->pram->oxidation_mech == "NSC_NEOH")
        return oxidation_NSC_Neoh();
    else if (domn->pram->oxidation_mech == "HACA")
        return oxidation_HACA(M0, M1);
    else {
        cout << endl << "ERROR: Invalid soot oxidation mechanism." << endl;
        exit(0);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/*! getCoagulationRate function
 *
 *      Calls appropriate function for coagulation chemistry
 *      based on the coagulation_mech flag.
 *      Returns the value of the collision rate function beta in m3/#*s.
 *
 *      @param m1      /input  mass of particle 1
 *      @param m2      /input  mass of particle 2
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::getCoagulationRate(const double &m1, const double &m2) {

    if (domn->pram->coagulation_mech      == "NONE")
        return 0;
    else if (domn->pram->coagulation_mech == "LL")
        return coagulation_LL(m1, m2);
    else if (domn->pram->coagulation_mech == "FUCHS")
        return coagulation_Fuchs(m1, m2);
    else if (domn->pram->coagulation_mech == "FRENK")
        return coagulation_Frenk(m1, m2);
    else {
        cout << endl << "ERROR: Invalid soot coagulation mechanism." << endl;
        exit(0);
    }

    return 0;

}

////////////////////////////////////////////////////////////////////////////////
/*! Nucleation by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical nucleation rate in #/m3*s.
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::nucleation_LL() {

    double cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];   // kmol/m3
    double Rnuc  =  0.1E5 * exp(-21100/T) * cC2H2;        // kmol/m^3*s

    rC2H2_rSoot_n  = -MW_sp[i_c2h2]/(2*MW_c);              // kg C2H2 / kg Soot
    rH2_rSoot_ncnd =  MW_sp[i_h2]  /(2*MW_c);              // kg H2   / kg Soot

    return Rnuc * 2 * Na / domn->pram->Cmin;              // #/m3*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Nucleation by Lindstedt (2005)
 *
 *      Rate from Lindstedt (2005), Proc. Comb. Inst. 30:775
 *      Uses Cmin = 10 for Naphthalene.
 *      Returns chemical nucleation rate in #/m3*s.
 *
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::nucleation_Linstedt() {

    double cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];   // kmol/m3
    double Rnuc  =  0.63E4 * exp(-21100/T) * cC2H2;       // kmol/m^3*s

    rC2H2_rSoot_n  = -MW_sp[i_c2h2]/(2*MW_c);              // kg C2H2 / kg Soot
    rH2_rSoot_ncnd =  MW_sp[i_h2]  /(2*MW_c);              // kg H2   / kg Soot

    return Rnuc * 2 * Na / domn->pram->Cmin;              // #/m3*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Helper function for PAH nucleation, and condensation
 *
 *      Sets the dimer mass m_dimer.
 *
 *      Rate from Blanquart & Pitsch (2009) article "A joint
 *      volume-surface-hydrogen multi-variate model for soot formation," ch. 27
 *      in Combustion Generated Fine Carbonaceous Particles ed. Bockhorn et al.
 *      Returns wdotD;
 *
 *      Call set_gas_state_vars first.
 *
 *      @param mi    /input vector of soot particle sizes
 *      @param wi    /input vector of soot particle weights
 *
 *      Note, Cmin is reset. (Some mechanisms have Cmin as an input).
 *      Note, the "preFac" is (1/2) the F.M. rate in coagulation_Frenk without eps_c.
 *
 */

double dv_soot::set_m_dimer() {

    //------------ compute wdotD, the dimer self collision rate

    double preFac = sqrt(4*M_PI*kb*T)*pow(6/(M_PI*rhoSoot), 2.0/3.0);
    double wdotD = 0.0;                      // dimer self collision rate (formation rate: #/m3*s)
    double gamma_i;                          // sticking coefficient
    double m_ipah;                           // PAH species mass per molecule
    double N_i;                              // PAH species number density: molecules / m3
    double wdoti = 0.0;                      // convenience variable
    m_dimer = 0.0;                           // dimer mass kg/part.
    for(int i=0; i<i_pah.size(); i++) {
        m_ipah  = MW_sp[i_pah[i]]/Na;
        gamma_i = m_ipah > 153 ? 1.501E-11*pow(m_ipah,4) : 1.501E-11*pow(m_ipah,4) / 3.0;
        N_i     = rho * (*yi)[i_pah[i]] / MW_sp[i_pah[i]] * Na;
        wdoti   = abs(gamma_i * preFac * pow(m_ipah, 1.0/6.0) * N_i*N_i);
        wdotD   += wdoti;
        m_dimer += wdoti*m_ipah;
        Cmin    += wdoti*domn->gas->nAtoms(i_pah[i],i_elem_c);
        rPAH_rSoot_ncnd[i] = wdoti*m_ipah;
    }
        for(int i=0; i<i_pah.size(); i++)
            rPAH_rSoot_ncnd[i] /= m_dimer;      // now mdot_i_pah = pah_relative_rates[i]*mdot, where mdot is a total gas rate
    m_dimer *= 2/wdotD;
    Cmin    *= 4/wdotD;                        // This is reset here. Some mechanisms have this as an input

    for(int i=0; i<i_pah.size(); i++)
        rPAH_rSoot_ncnd[i] *= -2.0*m_dimer/(Cmin*MW_c/Na);
    rH2_rSoot_ncnd =  2.0*m_dimer/(Cmin*MW_c/Na) - 1.0;

    return wdotD;
}

////////////////////////////////////////////////////////////////////////////////
/*! Helper function for PAH nucleation, and condensation
 *
 *      Sets the dimer number density DIMER (#/m3).
 *
 *      Rate from Blanquart & Pitsch (2009) article "A joint
 *      volume-surface-hydrogen multi-variate model for soot formation," ch. 27
 *      in Combustion Generated Fine Carbonaceous Particles ed. Bockhorn et al.
 *
 *      Call set_gas_state_vars first.
 *
 *      @param mi    /input vector of soot particle sizes
 *      @param wi    /input vector of soot particle weights
 *
 */

void dv_soot::set_Ndimer(const vector<double> &mi, const vector<double> &wi) {

    double wdotD = set_m_dimer();

    //------------- compute the dimer concentration as solution to quadratic
    // Steady state approximation.
    // Dimer creation rate = dimer destruction from self collision + from soot collision
    // wdotD = beta_DD*[D]^2 + sum(beta_DS*w_i)*[D]

    double beta_DD = coagulation_Frenk(m_dimer, m_dimer);          // dimer self-collision rate
    double I_beta_DS = 0.0;                                        // sum of dimer-soot collision rates
    for(int i=0; i<mi.size(); i++)                                 // loop over soot "particles" (abscissas)
        I_beta_DS += abs(wi[i]) * coagulation_Frenk(m_dimer, mi[i]);

    //------------- solve quadratic for D: beta_DD*(D^2) + I_beta_DS*(D) - wdotD = 0
    // See numerical recipies 3rd ed. sec 5.6 page 227.
    // Choosing the positive root.

    DIMER = 2.0*wdotD/(I_beta_DS + sqrt(I_beta_DS*I_beta_DS + 4*beta_DD*wdotD));   // #/m3

}

////////////////////////////////////////////////////////////////////////////////
/*! PAH nucleation by Blanquart et al. (2009)
 *
 *      Rate from Blanquart & Pitsch (2009) article "A joint
 *      volume-surface-hydrogen multi-variate model for soot formation," ch. 27
 *      in Combustion Generated Fine Carbonaceous Particles ed. Bockhorn et al.
 *      Returns chemical nucleation rate in #/m3*s.
 *
 *      Call set_gas_state_vars first.
 *
 *      @param mi    /input vector of soot particle sizes
 *      @param wi    /input vector of soot particle weights
 *
 *
 */

double dv_soot::nucleation_PAH(const vector<double> &mi, const vector<double> &wi) {

    set_Ndimer(mi, wi);
    double beta_DD = coagulation_Frenk(m_dimer, m_dimer);          // dimer self-collision rate

    return 0.5*beta_DD*DIMER*DIMER;                                // Jnuc (=) #/m3*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Growth by Lindstedt (1994)
 *
 *      Rate from Bockhorn (1994) pg. 417, "Simplified Soot Nucleation and Surface Growth Steps..."
 *      Equation (27.35).
 *      Returns chemical surface growth rate in kg/m2*s.
 *
 *      Call set_gas_state_vars first.
 *
 */

double dv_soot::growth_Lindstedt() {

    double cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];        // kmol/m3
    double rSoot = 750.0 * exp(-12100.0/T) * cC2H2 * 2.0*MW_c; // kg/m^2*s

    rC2H2_rSoot_go = -MW_sp[i_c2h2]/(2*MW_c);                      // kg C2H2 / kg Soot
    rH2_rSoot_go   =  MW_sp[i_h2]  /(2*MW_c);                      // kg H2   / kg Soot

    return rSoot;

}

////////////////////////////////////////////////////////////////////////////////
/*! Growth by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical surface growth rate in kg/m2*s.
 *
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::growth_LL(const double &M0, const double &M1) {

    double Am2m3 = 0.0;
    double cC2H2 = 0.0;
    double rSoot = 0.0;

    if (M0 > 0.0)
        Am2m3 = M_PI * pow(abs(6/(M_PI*rhoSoot)*M1/M0),2.0/3.0) * abs(M0);    // m^2_soot / m^3_total = pi*di^2*M0

    cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];                          // kmol/m3

    if (Am2m3 > 0)
        rSoot = 0.6E4 * exp(-12100.0/T) * cC2H2/sqrt(Am2m3) * 2.0*MW_c;       // kg/m^2*s

    rC2H2_rSoot_go = -MW_sp[i_c2h2]/(2*MW_c);                      // kg C2H2 / kg Soot
    rH2_rSoot_go   =  MW_sp[i_h2]  /(2*MW_c);                      // kg H2   / kg Soot

    return rSoot;

}

//////////////////////////////////////////////////////////////////////////////////
/*! Growth by HACA
 *
 *      See Appel, Bockhorn, & Frenklach (2000), Comb. & Flame 121:122-136.
 *      For details, see Franklach and Wang (1990), 23rd Symposium, pp. 1559-1566.
 *
 *      Parameters for steric factor alpha updated to those given in Balthasar
 *      and Franklach (2005) Comb. & Flame 140:130-145.
 *
 *      Call set_gas_state_vars first.
 *
 *      Returns the chemical soot growth rate in kg/m2*s.
 *
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 */

double dv_soot::growth_HACA(const double &M0, const double &M1) {

    double cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];      // kmol/m3
    double cO2   = rho * (*yi)[i_o2]   / MW_sp[i_o2];        // kmol/m3
    double cH    = rho * (*yi)[i_h]    / MW_sp[i_h];         // kmol/m3
    double cH2   = rho * (*yi)[i_h2]   / MW_sp[i_h2];        // kmol/m3
    double cOH   = rho * (*yi)[i_oh]   / MW_sp[i_oh];        // kmol/m3
    double cH2O  = rho * (*yi)[i_h2o]  / MW_sp[i_h2o];       // kmol/m3

    //---------- calculate alpha, other constants
    double RT       = 1.9872036E-3 * T;         // R (=) kcal/mol
    double chi_soot = 2.3E15;                   // (=) sites/cm^2
    double a_param  = 33.167 - 0.0154 * T;      // a parameter for calculating alpha
    double b_param  = -2.5786 + 0.00112 * T;    // b parameter for calculating alpha

    //---------- calculate raw HACA reaction rates
    double fR1 = 4.2E13 * exp(-13.0 / RT) * cH / 1000;
    double rR1 = 3.9E12 * exp(-11.0 / RT) * cH2 / 1000;
    double fR2 = 1.0E10 * pow(T,0.734) * exp(-1.43 / RT) * cOH / 1000;
    double rR2 = 3.68E8 * pow(T,1.139) * exp(-17.1 / RT) * cH2O /1000;
    double fR3 = 2.0E13 * cH / 1000;
    double fR4 = 8.00E7 * pow(T,1.56) * exp(-3.8 / RT) * cC2H2 / 1000;
    double fR5 = 2.2E12 * exp(-7.5 / RT) * cO2 / 1000;
    double fR6 = 1290.0 * 0.13 * P * (cOH/rho*MW_sp[i_oh]) / sqrt(T);    // gamma = 0.13 from Neoh et al.

    //---------- Steady state calculation of chi for soot radical; see Frenklach 1990 pg. 1561
    double denom = rR1 + rR2 + fR3 + fR4 + fR5;
    double chi_rad = 0.0;
    if(denom != 0.0)
        chi_rad = 2 * chi_soot * (fR1 + fR2 + fR6) / denom;        // sites/cm^2

    double alpha = 1.0;                                            // alpha = fraction of available surface sites
    if (M0 > 0.0)
        alpha = tanh(a_param/log10(M1/M0)+b_param);
    if (alpha < 0.0)
        alpha = 1.0;

    double c_soot_H   = alpha * chi_soot * 1E4;                    // sites/m2-mixture
    double c_soot_rad = alpha * chi_rad  * 1E4;                    // sites/m2-mixture

    return (fR5*c_soot_rad + fR6*c_soot_H) / Na * 2 * MW_c;        // kg/m2*s
}

////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by Leung_Lindstedt (1991)
 *
 *      Rate from Leung & Lindstedt (1991), Comb. & Flame 87:289-305.
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      C + 0.5 O2 --> CO
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::oxidation_LL() {

    double cO2 = rho * (*yi)[i_o2] / MW_sp[i_o2];             // kmol/m3
    double rSoot = 0.1E5 * sqrt(T) * exp(-19680.0/T) * cO2 * MW_c;    // kg/m^2*s

    rO2_rSoot_go = -0.5*MW_sp[i_o2]/MW_c;                        // kg O2 / kg Soot
    rCO_rSoot_go =      MW_sp[i_co]/MW_c;                        // kg CO / kg Soot

    return rSoot;

}

////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by Lee et al. + Neoh
 *
 *      Rates from Lee et al. (1962) Comb. & Flame 6:137-145 and Neoh (1981)
 *      "Soot oxidation in flames" in Particulate Carbon Formation During
 *      Combustion book
 *      C + 0.5 O2 --> CO
 *      C + OH     --> CO + H
 *
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::oxidation_Lee_Neoh() {

    double pO2 = (*yi)[i_o2] * MW / MW_sp[i_o2] * P / 101325.0;      // partial pressure of O2 (atm)
    double pOH = (*yi)[i_oh] * MW / MW_sp[i_oh] * P / 101325.0;      // partial pressure of OH (atm)

    double rSootO2 = 1.085E4*pO2/sqrt(T)*exp(-1.977824E4/T)/1000.0;  // kg/m^2*s
    double rSootOH = 1290.0*0.13*pOH/sqrt(T);                        // kg/m^2*s

    rO2_rSoot_go = -0.5*MW_sp[i_o2]/MW_c * rSootO2/(rSootO2+rSootOH);    // kg O2 / kg Soot
    rOH_rSoot_go =     -MW_sp[i_oh]/MW_c * rSootOH/(rSootO2+rSootOH);    // kg OH / kg Soot
    rH_rSoot_go  =      MW_sp[i_h] /MW_c * rSootOH/(rSootO2+rSootOH);    // kg H  / kg Soot
    rCO_rSoot_go =      MW_sp[i_co]/MW_c;                                // kg CO / kg Soot

    return rSootO2 + rSootOH;

}

////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by NSC + Neoh
 *
 *      Rates from Nagle and Strickland-Constable (1961) and Neoh (1981) "Soot
 *      oxidation in flames" in Particulate Carbon Formation During Combustion
 *      book
 *      C + 0.5 O2 --> CO
 *      C + OH     --> CO + H
 *
 *      Returns chemical soot oxidation rate in kg/m2*s.
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::oxidation_NSC_Neoh() {

    double pO2 = (*yi)[i_o2] * MW / MW_sp[i_o2] * P / 101325.0; // partial pressure of O2 (atm)
    double pOH = (*yi)[i_oh] * MW / MW_sp[i_oh] * P / 101325.0; // partial pressure of OH (atm)

    double kA = 20.0     * exp(-15098.0/T);                     // rate constants
    double kB = 4.46E-3  * exp(-7650.0/T);
    double kT = 1.51E5   * exp(-48817.0/T);
    double kz = 21.3     * exp(2063.0/T);

    double x  = 1.0/(1.0+kT/(kB*pO2));                          // x = unitless fraction
    double NSC_rate = kA*pO2*x/(1.0+kz*pO2) + kB*pO2*(1.0-x);   // kmol/m^2*s
    double rSootO2 = NSC_rate*rhoSoot;                          // kg/m2*s
    double rSootOH = 1290.0*0.13*pOH/sqrt(T);                   // kg/m2*s

    rO2_rSoot_go = -0.5*MW_sp[i_o2]/MW_c * rSootO2/(rSootO2+rSootOH);    // kg O2 / kg Soot
    rOH_rSoot_go =     -MW_sp[i_oh]/MW_c * rSootOH/(rSootO2+rSootOH);    // kg OH / kg Soot
    rH_rSoot_go  =      MW_sp[i_h] /MW_c * rSootOH/(rSootO2+rSootOH);    // kg H  / kg Soot
    rCO_rSoot_go =      MW_sp[i_co]/MW_c;                                // kg CO / kg Soot

    return rSootO2 + rSootOH;                                   // kg/m2*s

}

//////////////////////////////////////////////////////////////////////////////////
/*! Oxidation by HACA
 *
 *      See Appel, Bockhorn, & Frenklach (2000), Comb. & Flame 121:122-136.
 *      For details, see Franklach and Wang (1990), 23rd Symposium, pp. 1559-1566.
 *
 *      Parameters for steric factor alpha updated to those given in Balthasar
 *      and Franklach (2005) Comb. & Flame 140:130-145.
 *
 *      Returns the chemical soot oxidation rate in kg/m2*s.
 *
 *      Call set_gas_state_vars first.
 *
 *      @param M0       /input  local soot number density (#/m3)
 *      @param M1       /input  local soot mass density (kg/m3)
 */

double dv_soot::oxidation_HACA(const double &M0, const double &M1) {

    double cC2H2 = rho * (*yi)[i_c2h2] / MW_sp[i_c2h2];      // kmol/m3
    double cO2   = rho * (*yi)[i_o2]   / MW_sp[i_o2];        // kmol/m3
    double cH    = rho * (*yi)[i_h]    / MW_sp[i_h];         // kmol/m3
    double cH2   = rho * (*yi)[i_h2]   / MW_sp[i_h2];        // kmol/m3
    double cOH   = rho * (*yi)[i_oh]   / MW_sp[i_oh];        // kmol/m3
    double cH2O  = rho * (*yi)[i_h2o]  / MW_sp[i_h2o];       // kmol/m3

    //---------- calculate alpha, other constants
    double RT       = 1.9872036E-3 * T;         // R (=) kcal/mol
    double chi_soot = 2.3E15;                   // (=) sites/cm^2
    double a_param  = 33.167 - 0.0154 * T;      // a parameter for calculating alpha
    double b_param  = -2.5786 + 0.00112 * T;    // b parameter for calculating alpha

    //---------- calculate raw HACA reaction rates
    double fR1 = 4.2E13 * exp(-13.0 / RT) * cH / 1000;
    double rR1 = 3.9E12 * exp(-11.0 / RT) * cH2 / 1000;
    double fR2 = 1.0E10 * pow(T,0.734) * exp(-1.43 / RT) * cOH / 1000;
    double rR2 = 3.68E8 * pow(T,1.139) * exp(-17.1 / RT) * cH2O /1000;
    double fR3 = 2.0E13 * cH / 1000;
    double fR4 = 8.00E7 * pow(T,1.56) * exp(-3.8 / RT) * cC2H2 / 1000;
    double fR5 = 2.2E12 * exp(-7.5 / RT) * cO2 / 1000;
    double fR6 = 1290.0 * 0.13 * P * (cOH/rho*MW_sp[i_oh]) / sqrt(T);    // gamma = 0.13 from Neoh et al.

    //---------- Steady state calculation of chi for soot radical; see Frenklach 1990 pg. 1561
    double denom = rR1 + rR2 + fR3 + fR4 + fR5;
    double chi_rad = 0.0;
    if(denom != 0.0)
        chi_rad = 2 * chi_soot * (fR1 + fR2 + fR6) / denom;        // sites/cm^2

    double alpha = 1.0;                                            // alpha = fraction of available surface sites
    if (M0 > 0.0)
        alpha = tanh(a_param/log10(M1/M0)+b_param);
    if (alpha < 0.0)
        alpha = 1.0;

    double c_soot_H   = alpha * chi_soot * 1E4;                    // sites/m2-mixture
    double c_soot_rad = alpha * chi_rad  * 1E4;                    // sites/m2-mixture

    double Roxi = -fR1*c_soot_H + rR1*c_soot_rad - fR2*c_soot_H + rR2*c_soot_rad +
                   fR3*c_soot_rad + fR4*c_soot_rad - fR6*c_soot_H; // #-available-sites/m2-mix*s
    return Roxi / Na * MW_c;                                       // kg/m2*s

}

////////////////////////////////////////////////////////////////////////////////
/*! Coagulation by Leung_Lindstedt
 *
 *      @param m1   \input  first particle size (kg)
 *      @param m2   \input  second particle size (kg)
 *
 *      Returns the value of the collision rate function beta in m3/#*s.
 *
 *      Note, this assumes the free molecular regime.
 *      The original LL model is for monodispersed and has the form
 *         2*Ca*sqrt(dp)*sqrt(6*kb*T/rhoSoot)
 *         This is Eq. (4) in LL but LL is missing the 1/2 power on (6*kb*T/rhoSoot)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::coagulation_LL (const double &m1, const double &m2) {

    const double Ca = 9.0;

    //--------- Free molecular form from Fuchs and/or Frenklack below.
    //double Dp1 = pow(6.0*abs(m1)/M_PI/rhoSoot, 1.0/3.0);
    //double Dp2 = pow(6.0*abs(m2)/M_PI/rhoSoot, 1.0/3.0);
    //double m12 = abs(m1*m2/(m1+m2));
    //return Ca/2.0*sqrt(M_PI*kb*T*0.5/m12) * pow(Dp1+Dp2, 2.0);

    //--------- Equivalent L&L form assuming m1 = m2
    double Dp1 = pow(6.0*abs(m1)/M_PI/rhoSoot, 1.0/3.0);
    return 2.0*Ca*sqrt(Dp1*6*kb*T/rhoSoot);

}

////////////////////////////////////////////////////////////////////////////////
/*! Coagulation by Fuchs
 *
 *      Rate comes from Seinfeld and Pandis Atmospheric Chemistry book (2016), pg. 548, chp 13.
 *      See also chapter 9.
 *      Details and clarification in Fuchs' Mechanics of Aerosols book (1964)
 *      Seinfeld is missing the sqrt(2) in the final term for g. This is needed to reproduce his plot.
 *      Fuchs' book has the sqrt(2).
 *      (I've seen another book https://authors.library.caltech.edu/25069/7/AirPollution88-Ch5.pdf
 *      with 1.0 in place of both sqrt(2) factors. This gives 5% max error in the D=D curve
 *
 *      Returns the value of the collision rate function beta in m3/#*s.
 *
 *      @param m1       \input  first particle size (kg)
 *      @param m2       \input  second particle size (kg)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::coagulation_Fuchs(const double &m1, const double &m2) {

    double Dp1 = pow(6.0*abs(m1)/M_PI/rhoSoot, 1.0/3.0);
    double Dp2 = pow(6.0*abs(m2)/M_PI/rhoSoot, 1.0/3.0);

    double c1 = sqrt(8.0*kb*T/M_PI/m1);
    double c2 = sqrt(8.0*kb*T/M_PI/m2);

    double mfp_g = get_gas_mean_free_path();

    double Kn1 = 2.0*mfp_g/Dp1;
    double Kn2 = 2.0*mfp_g/Dp2;

    double Cc1 = 1 + Kn1*(1.257 + 0.4*exp(-1.1/Kn1));    // Seinfeld p. 372 eq. 9.34. This is for air at 298 K, 1 atm
    double Cc2 = 1 + Kn2*(1.257 + 0.4*exp(-1.1/Kn2));    // for D<<mfp_g, Cc = 1 + 1.657*Kn; Seinfeld p. 380: 10% error at Kn=1, 0% at Kn=0.01, 100

    double D1 = kb*T*Cc1/(3.0*M_PI*mu*Dp1);
    double D2 = kb*T*Cc2/(3.0*M_PI*mu*Dp2);

    double l1 = 8.0*D1/M_PI/c1;
    double l2 = 8.0*D2/M_PI/c2;

    double g1 = sqrt(2.0)/3.0/Dp1/l1*( pow(Dp1+l1,3.0) - pow(Dp1*Dp1 + l1*l1, 3.0/2.0) ) - sqrt(2.0)*Dp1;
    double g2 = sqrt(2.0)/3.0/Dp2/l2*( pow(Dp2+l2,3.0) - pow(Dp2*Dp2 + l2*l2, 3.0/2.0) ) - sqrt(2.0)*Dp2;

    return 2.0*M_PI*(D1+D2)*(Dp1+Dp2) / ((Dp1+Dp2)/(Dp1+Dp2+2.0*sqrt(g1*g1+g2*g2)) + 8.0/eps_c*(D1+D2)/sqrt(c1*c1+c2*c2)/(Dp1+Dp2));

}

////////////////////////////////////////////////////////////////////////////////
/*! Coagulation by Frenklach
 *
 *      Returns the value of the collision rate function beta in m3/#*s.
 *
 *      @param m1       \input  first particle size (kg)
 *      @param m2       \input  second particle size (kg)
 *
 *      Call set_gas_state_vars first.
 */

double dv_soot::coagulation_Frenk(const double &m1, const double &m2) {

    double Dp1 = pow(6.0*abs(m1)/M_PI/rhoSoot, 1.0/3.0);
    double Dp2 = pow(6.0*abs(m2)/M_PI/rhoSoot, 1.0/3.0);

    //------------ free molecular rate

    double m12 = abs(m1*m2/(m1+m2));

    double beta_12_FM = eps_c*sqrt(M_PI*kb*T*0.5/m12) * pow(Dp1+Dp2, 2.0);

    //------------ continuum rate

    double mfp_g = get_gas_mean_free_path();

    double Kn1 = 2.0*mfp_g/Dp1;
    double Kn2 = 2.0*mfp_g/Dp2;

    double Cc1 = 1 + Kn1*(1.257 + 0.4*exp(-1.1/Kn1));    // Seinfeld p. 372 eq. 9.34. This is for air at 298 K, 1 atm
    double Cc2 = 1 + Kn2*(1.257 + 0.4*exp(-1.1/Kn2));    // for D<<mfp_g, Cc = 1 + 1.657*Kn; Seinfeld p. 380: 10% error at Kn=1, 0% at Kn=0.01, 100

    double beta_12_C = 2*kb*T/(3*mu)*(Cc1/Dp1 + Cc2/Dp2)*(Dp1 + Dp2);

    //------------ return harmonic mean

    return beta_12_FM * beta_12_C / (beta_12_FM + beta_12_C);

}

////////////////////////////////////////////////////////////////////////////////
/*! Gas mean free path
 *      Returns the value of the collision rate function beta in m3/#*s.
 *      Call set_gas_state_vars first.
 */

double dv_soot::get_gas_mean_free_path() {
    return mu/rho*sqrt(M_PI*MW/(2.0*Rg*T));
}

////////////////////////////////////////////////////////////////////////////////
/*! Kc
 *      Returns continuum coagulation coefficient Kc
 *      Call set_gas_state_vars first.
 */

double dv_soot::get_Kc() {
    return 2.0*kb*T/(3.0*mu);
}

////////////////////////////////////////////////////////////////////////////////
/*! Kcp
 *      Returns continuum coagulation coefficient Kc prime
 *      Call set_gas_state_vars first.
 */

double dv_soot::get_Kcp() {
    return 2.0*1.657*get_gas_mean_free_path()*pow(M_PI/6*rhoSoot,1./3.);
}

////////////////////////////////////////////////////////////////////////////////
/*! Kfm
 *      Returns continuum coagulation coefficient Kc prime
 *      Call set_gas_state_vars first.
 */

double dv_soot::get_Kfm() {
    return eps_c*sqrt(M_PI*kb*T/2)*pow(6./M_PI/rhoSoot,2./3.);
}

////////////////////////////////////////////////////////////////////////////////
/*! set_gasSootSources
 *      Sets the gas source terms from reaction with soot
 *      Call set_gas_state_vars first.
 *
 *      @param N1       \input  soot M1 source term (kg/m3*s) from nucleation
 *      @param Cnd1     \input  soot M1 source term (kg/m3*s) from condensation
 *      @param G1       \input  soot M1 source term (kg/m3*s) from growth
 *      @param X1       \input  soot M1 source term (kg/m3*s) from oxidation
 *      @param igrd     \input  grid point
 *
 */

void dv_soot::set_gasSootSources(const double &N1, const double &Cnd1, const double &G1, const double &X1, const int igrd) {


    //---nucleation: see dv_soot.cc for rC2H2_rSoot_n, etc.
    gasSootSources[i_c2h2][igrd] = N1 * rC2H2_rSoot_n  / rho;  // some of these terms might be 0 dep. on mechanisms used.
    gasSootSources[i_h2][igrd]   = N1 * rH2_rSoot_ncnd / rho;  // signs are ebedded in terms
    for(int i=0; i<i_pah.size(); i++)
        gasSootSources[i_pah[i]][igrd] = N1 * rPAH_rSoot_ncnd[i] / rho;

    //---growth

    gasSootSources[i_c2h2][igrd] += G1 * rC2H2_rSoot_go / rho;
    gasSootSources[i_h2][igrd]   += G1 * rH2_rSoot_go   / rho;

    //---oxidation

    gasSootSources[i_o2][igrd]   = X1 * rO2_rSoot_go   / rho;
    gasSootSources[i_oh][igrd]   = X1 * rOH_rSoot_go   / rho;
    gasSootSources[i_h ][igrd]   = X1 * rH_rSoot_go    / rho;
    gasSootSources[i_co][igrd]   = X1 * rCO_rSoot_go   / rho;

    //---PAH condensation

    gasSootSources[i_h2][igrd]   += Cnd1 * rH2_rSoot_ncnd  / rho;
    for(int i=0; i<i_pah.size(); i++)
        gasSootSources[i_pah[i]][igrd] += Cnd1 * rPAH_rSoot_ncnd[i] / rho;

    //---coagulation: NA

}

