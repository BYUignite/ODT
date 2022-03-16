/**
 * @file domain.cc
 * @brief Source file for class \ref domain
 */

#include "domain.h"
#include "processor.h"
#include "dv.h"
#include "dv_uvw.h"
#include "dv_pos.h"
#include "dv_posf.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_sca.h"
#include "dv_aDL.h"
#include "domaincase_odt_channel.h"
#include "domaincase_odt_channelScalar.h"
#include "domaincase_odt_isothermalWall.h"
#include "domaincase_odt_jetMixlRxn.h"
#include "domaincase_odt_jetFlame.h"
#include "domaincase_odt_MFjetFlame.h"
#include "domaincase_odt_coldPropaneJet.h"
#include "domaincase_odt_coldJet.h"
#include "domaincase_odt_RT.h"
#include "chemical_mechanisms/onestep_ch4.h"
#include "chemical_mechanisms/fourstep_ch4.h"
#include "chemical_mechanisms/onestep_c2h4.h"
#include "chemical_mechanisms/simple_dlr.h"
#include "chemical_mechanisms/ch4red.h"
#include "chemical_mechanisms/c2h4red.h"
#include "chemical_mechanisms/canteraRR.h"
#include "chemical_mechanisms/chemNone.h"
#include <cmath>
#include <iomanip>

extern processor proc;

/////////////////////////////////////////////////////////////////////
/** Constructor
 */

domain::domain(domain *p_domn, param *p_pram) {

    domn = p_domn;
    pram = p_pram;
    domc = 0;               // initialize for destruction of eddy domains

 }

/////////////////////////////////////////////////////////////////////
/** Initializer
 */

void domain::init(inputoutput     *p_io,
                  meshManager     *p_mesher,
                  streams         *p_strm,
//                  IdealGasPhase   *p_gas,
//                  Transport       *p_tran,
                  shared_ptr<Solution> p_gas,
                  micromixer      *p_mimx,
                  eddy            *p_ed,
                  domain          *p_eddl,
                  solver          *p_solv,
                  randomGenerator *p_rand,
                  chemMech        *p_chem,
                  bool             LisEddyDomain) {

    //----------------------

    io     = p_io;
    mesher = p_mesher;
    gas    = p_gas;
//    tran   = p_tran;
    strm   = p_strm;
    mimx   = p_mimx;
    ed     = p_ed;
    eddl   = p_eddl;
    solv   = p_solv;
    rand   = p_rand;
    chem   = p_chem;

    //----------------------

    ngrd    = pram->ngrd0;
    ngrdf   = ngrd + 1;

    //----------------------

    if(LisEddyDomain) {        // eddy domain needs less data
        initEddyDomain();
        return;
    }

    //----------------------
    io->init(this);
    pram->init(this);
    ed->init(this, eddl);
    solv->init(this);
    // mesher is init below in caseinit for phi
    // strm is init below in caseinit  (domc), (if needed)
    // mimx is init below since it needs v[] set for cvode

    //---------------------- Continue setting up the case using the case_somecase class.
    // Adds to the above variable list, and initializes solution for the run

     if(pram->probType == "CHANNEL")
         domc = new domaincase_odt_channel();    // cold channel flow
     else if(pram->probType == "CHANNEL_SCALAR")
         domc = new domaincase_odt_channelScalar();  // cold channel flow with passive scalar
     else if(pram->probType == "JETMIXL_RXN")
         domc = new domaincase_odt_jetMixlRxn(); // jet, wake, mixing layer with gaseous reaction
     else if(pram->probType == "COLDPROPANEJET")
         domc = new domaincase_odt_coldPropaneJet(); // TNF jet
     else if(pram->probType == "COLDJET")
         domc = new domaincase_odt_coldJet(); // Hussein 1994
     else if(pram->probType == "JETFLAME")
         domc = new domaincase_odt_jetFlame(); // Shaddix jet
     else if(pram->probType == "MF_JETFLAME")
         domc = new domaincase_odt_MFjetFlame(); // jet flame w/ mixt frac density profile
     else if(pram->probType == "ISOTHERMAL_WALL")
         domc = new domaincase_odt_isothermalWall(); // isothermal wall
     else if(pram->probType == "RT")
         domc = new domaincase_odt_RT();      // simple Rayleigh Taylor flow
     else {
         cout << endl << "ERROR, probType UNKNOWN" << endl;
         exit(0);
     }

    domc->init(this);

    //----------------------

    if (pram->chemMech == "simple_dlr")
        chem = new simple_dlr(this);
    else if (pram->chemMech == "onestep_ch4")
        chem = new onestep_ch4(this);
    else if (pram->chemMech == "onestep_c2h4")
        chem = new onestep_c2h4(this);
    else if (pram->chemMech == "fourstep_ch4")
        chem = new fourstep_ch4(this);
    else if (pram->chemMech == "c2h4red")
        chem = new c2h4red(this);
    else if (pram->chemMech == "ch4red")
        chem = new ch4red(this);
    else if (pram->chemMech == "gri30" || pram->chemMech == "gri30_highT" || pram->chemMech == "air")
        chem = new canteraRR(this);
    else if (pram->chemMech == "none")
        chem = new chemNone(this);
    else
        throw domain_error("Unknown chemical mechanism requested");

    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);

    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------

    mimx->init(this);

    //----------------------

    if(pram->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(pram->trst);
    }

}

/////////////////////////////////////////////////////////////////////
/** Compute size of domain based on faces.
 */

double domain::Ldomain() {
     return posf->d.at(ngrd) - posf->d.at(0);
}

/////////////////////////////////////////////////////////////////////
/** Initialize data members of the eddy domain.
 *  Note, none of the other members of this domain should be used (like random).
 *  Note, all variables here should have corresponding variables (by var_name) in the
 *     main domn. This is needed for using the eddl object.
 */

void domain::initEddyDomain() {

    v.push_back(new dv_pos(  this, "pos",   false, true));
    v.push_back(new dv_posf( this, "posf",  false, true));
    v.push_back(new dv_uvw(  this, "uvel",  true,  true));   // last are: L_transported, L_output
    v.push_back(new dv_uvw(  this, "vvel",  true,  true));
    v.push_back(new dv_uvw(  this, "wvel",  true,  true));
    v.push_back(new dv_rho(  this, "rho",   false, false));
    v.push_back(new dv_dvisc(this, "dvisc", false, false));
    if(domn->pram->LdoDL)
       v.push_back(new dv_aDL(this, "aDL",   false, false));

    int k = 0;
    pos   = v.at(k++);
    posf  = v.at(k++);
    uvel  = v.at(k++);
    vvel  = v.at(k++);
    wvel  = v.at(k++);
    rho   = v.at(k++);
    dvisc = v.at(k++);
    if(domn->pram->LdoDL)
        aDL   = v.at(k++);

}

/////////////////////////////////////////////////////////////////////
/** Set the domain from a region of the domn.  Normally called by eddy domain.
 *  @param i1 \input index of starting cell of domn to build from
 *  @param i2 \input index of ending cell of domn to build from
 *  If i2 < i1, we have a periodic region (wrap around the domain).
 *     This only happens in planar cases, not cylindrical or sphericial.
 *  nonwrap: |   | * | * | * | * | * | * |   |   |
 *                i1                  i2
 *  new domain consists of *'d cells
 *
 *  Wrap: | 4 | 5 |   |   |   |   | 1 | 2 | 3 |
 *             i2                  i1
 *  New domain consists of #'d cells: 1 2 3 4 5}
 */

void domain::setDomainFromRegion(const int i1, const int i2) {

    ngrd  = i2-i1+1;
    ngrdf = ngrd+1;

    for(int k=0; k<v.size(); k++)
        v.at(k)->setDvFromRegion(i1,i2);
}

///////////////////////////////////////////////////////////////////////////////
/** Find index of cell for given position (residing in cell).
 *  Start search assuming a uniform grid,
 *  then search forward or back till hit the cell index.
 *  If position is on cell face j, then if LowSide true, return j, else j-1.         \n
 * For start of eddy region, set LowSide to true                                     \n
 * For end of eddy region, set LowSide to false                                      \n
 * (This is so triplet maps don't overlap cells)
 *                                                                                   <pre><code>
 * e.g., usual:   | { | | | | } |    5 pts, eddy pos between cell faces
 *       okay:    {   | | | |   }    5 pts, eddy pos on cell faces (1 or both)
 *       bad:     |   { | | }   |    5 pts, eddy pos on internal faces (1 or both)
 *                                                                                   </code></pre>
 * @param position \input position to find corresponding index.
 * @param LowSide  \input flag true, then return j if position is on cell face j, else j-1.
 * @return index of position.
 */

int domain::domainPositionToIndex(double position, const bool LowSide, int dbg) {

    if(abs(position-posf->d.at(0)) < 1.0E-14)
        return 0;
    if(abs(position-posf->d.at(ngrd)) < 1.0E-14)
        return ngrd-1;

    //if(position < posf->d.at(0))         // for periodic (from eddies only)
    //    position += Ldomain();
    //if(position > posf->d.at(ngrd))
    //    position -= Ldomain();

    if(position < posf->d.at(0) || position > posf->d.at(ngrd)) {
       *io->ostrm << "\ndbg = " << dbg << endl; //doldb
       *io->ostrm << scientific;
       *io->ostrm << setprecision(14);
       *io->ostrm << "\n ERROR odt_grid::domainPositionToIndex position < posf->d.at(0) or > posf->d.at(ngrd) \n"
               " and at processor's id---> " << proc.myid
               <<" Value of position is---> "<<position << " and values of posf->d.at(0) and posf->d.at(ngrd) are "
               <<posf->d.at(0)<< " and "<<posf->d.at(ngrd) <<" respectively "<< endl;
       //io->outputProperties("dbg.dat", 0.0); //doldb
       exit(0);
    }

    int i;
    int ipos = static_cast<int>((position-posf->d.at(0))/Ldomain()*ngrd);

    if(posf->d.at(ipos+1) > position) {      // case 1: grd skewed more pts on right half
        for(i=ipos+1; i>=0; i--)  {
            if(posf->d.at(i) <= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i;
            }
        }
    }

    else  {                           // case 2: grd skewed more pts on left half
        for(i=ipos+1; i<=ngrdf; i++) {
            if(posf->d.at(i) >= position) {
                if(position == posf->d.at(i)) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i-1;
            }
        }
    }

    *io->ostrm << "\n\n******** ERROR IN odt_grid::domainPositionToIndex "
         << position << '\t' << posf->d.at(0) << '\t' << posf->d.at(ngrd) << '\t' << endl << endl;

    return -1;
}

/////////////////////////////////////////////////////////////////////
/** Cycle domain for periodic flows.
 *  @param icycle \input move all cells before and including this one
 *   to the end of the domain.
 *  @return the cycle distance (used for backcycling).
 */

double domain::cyclePeriodicDomain(const int icycle) {

    double cycleDistance = posf->d.at(icycle+1)-posf->d.at(0);

    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
        v.at(k)->d.insert(v.at(k)->d.end(),   v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
        v.at(k)->d.erase( v.at(k)->d.begin(), v.at(k)->d.begin()+icycle+1);
    }

    //---------- now do posf, and pos

    double xend = posf->d.at(ngrd);
    for(int i=1; i<=icycle+1; i++)
        posf->d.push_back(xend+(posf->d.at(i)-posf->d.at(0)));
    posf->d.erase(posf->d.begin(), posf->d.begin()+icycle+1);

    pos->setVar();     // does a little extra work (whole domain) but doesn't happen that often
                       //    only when periodic eddies are accepted.

    return cycleDistance;
}

/////////////////////////////////////////////////////////////////////
/** Back cycle domain for periodic flows. Intended to be called some time
 *  after cyclePeriodicDomain is called.
 *  Splits the cell at posf.at(ngrd) - backCycleDistace, then moves end cells
 *     after the split to the beginning of the domain.
 *  @param \input distance from the end to split and move the domain.
 */

void domain::backCyclePeriodicDomain(const double backCycleDistance) {

    double xend = posf->d.at(ngrd) - backCycleDistance;     // end loc.
    double icycle = domainPositionToIndex(xend, true, 1);  // cycle cells greater than this to beginning

    //------------ split the cell where the back cycle happens

    vector<double> interPos(3);
    if(abs(posf->d.at(icycle) - xend) > 1.0e-15) {
        interPos.at(0) = posf->d.at(icycle);
        interPos.at(1) = xend;
        interPos.at(0) = posf->d.at(icycle+1);
        mesher->splitCell(icycle, 1, interPos);
        icycle++;
    }

    //------------ now move the cells

    int nmove = ngrd-icycle+1;

    for(int k=0; k<v.size(); k++) {
        if (v.at(k)->var_name=="pos" || v.at(k)->var_name=="posf")
            continue;
        v.at(k)->d.insert(v.at(k)->d.begin(), v.at(k)->d.begin()+icycle, v.at(k)->d.end());
        v.at(k)->d.erase( v.at(k)->d.begin()+icycle+nmove, v.at(k)->d.end() );
    }

    //---------- now do posf, and pos

    double xstart_orig = posf->d.at(0);

    posf->d.insert(posf->d.begin(), nmove, 0.0);
    icycle += nmove;
    for(int i=0; i<nmove; i++)
        posf->d.at(i) = xstart_orig - (posf->d.at(posf->d.size()-1-i) - xend);

    posf->d.erase(posf->d.begin()+icycle+1, posf->d.end());

    pos->setVar();     // does a little extra work (whole domain) but doesn't happen that often
                       //    only when periodic eddies are accepted.
}