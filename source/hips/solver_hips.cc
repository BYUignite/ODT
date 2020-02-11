
#include "solver_hips.h"
#include "domain.h"
#include <cmath>     // pow, log, ceil
#include <algorithm> // copy, sort
#include <sstream>
#include <string>
#include <iomanip>   // precision

///////////////////////////////////////////////////////////////////////////////

bool sortFunc(pair<double,int> &a, pair<double,int> &b){
    return a.first < b.first;
}

///////////////////////////////////////////////////////////////////////////////
/** Initializer
    @param p_domn \input pointer to domain object
 */

void solver_hips::init(domain *p_domn) {

    domn = p_domn;

    //------------------- Set number of parcels, and level lengthscales, timescales, and rates

    iEta = domn->pram->nLevels - 3;
    if(domn->pram->LScHips){
        int maxSc=1.0;
        for(int i=0; i<domn->io->scalarSc.size(); i++)
            maxSc = domn->io->scalarSc[i].as<double>() > maxSc ? domn->io->scalarSc[i].as<double>() : maxSc;
        if(maxSc > 1.0) {
            domn->pram->nLevels += ceil(log(maxSc)/log(4));
            cout << endl << endl << "maxSc > 1.0, increasing nLevels to " << domn->pram->nLevels << endl << endl;
        }
    }

    Nm1 = domn->pram->nLevels - 1;
    Nm3 = domn->pram->nLevels - 3;

    domn->ngrd = static_cast<int>(pow(2, Nm1));

    vector<double> levelLengths(domn->pram->nLevels);    // including all levels, but last 2 don't count:
    vector<double> levelTaus(domn->pram->nLevels);       //     smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(domn->pram->nLevels);

    for(int i=0; i<domn->pram->nLevels; i++){
        levelLengths[i] = domn->pram->domainLength * pow(domn->pram->Afac,i);
        levelTaus[i]    = domn->pram->tau0 *
                          pow(levelLengths[i]/domn->pram->domainLength, 2.0/3.0) /
                          domn->pram->C_param;
        levelRates[i]   = 1.0/levelTaus[i] * pow(2.0,i);
    }
    if(domn->pram->LScHips){     // correct levels for high Sc (levels above Kolmogorov)
        for(int i=iEta+1; i<domn->pram->nLevels; i++) {
            levelTaus[i]    = domn->pram->tau0 *
                              pow(levelLengths[iEta]/domn->pram->domainLength, 2.0/3.0) /
                              domn->pram->C_param;
            levelRates[i]   = 1.0/levelTaus[i] * pow(2.0,i);
        }
    }

    //-------------------

    totalEddyRate = 0.0;
    for(int i=0; i<=Nm3; i++)
        totalEddyRate += levelRates[i];
    eddyLevelCDF.resize(domn->pram->nLevels - 2);
    eddyLevelCDF[0] = levelRates[0]/totalEddyRate;
    for(int i=1; i<=Nm3; i++)
        eddyLevelCDF[i] = eddyLevelCDF[i-1] + levelRates[i]/totalEddyRate;

    //-------------------

    tMix = domn->pram->fmix * levelTaus[domn->pram->nLevels-2];

    //------------------- Set the parcel addresses (index array)

    pLoc.resize(domn->ngrd);
    for(int i=0; i<domn->ngrd; i++)
        pLoc[i] = i;

    //------------------- set the eddy event times at all levels

    setEddyEventTimes();

    //------------------- verify settings

    if(domn->pram->LScHips && !domn->pram->LsimpleMix){
        cout << endl << "\nERROR: LScHips requires LsimpleMix\n" << endl;
        exit(0);
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Sets arrays s.eTimes and s.eLevels which hold the eddy event times and the corresponding tree base level.
    First make an array of arrays for times at each level.
    These are from a Poisson process with the given rate at each level.
    Then collapse these into a single array.
    Then sort these times, along with the corresponding level array.
 */

void solver_hips::setEddyEventTimes(){

    double dmb;

    for(int i=0; i<levelRates.size()-2; i++) {
        double rate = levelRates[i];
        double time = 0.0;
        for(;;) {
            double r = domn->rand->getRand();
            time += -log(r)/rate;
            if(time > domn->pram->tEnd)
                break;
            eTL.push_back(make_pair(time, i));
        }
    }

    //--------- now sort the list of times and the associated levels

    sort(eTL.begin(), eTL.end(), sortFunc);

}
///////////////////////////////////////////////////////////////////////////////

void solver_hips::sample_hips_eddy(double &dt, double &iLevel) {

    //--------------- time to next eddy

    double r = domn->rand->getRand();
    dt = -log(r)/totalEddyRate;

    //--------------- eddy level: option 1

    static double a = 1/totalEddyRate/domn->pram->tau0;
    static double b = 2.0/pow(domn->pram->Afac, 2.0/3.0);
    static double coef1 = (1-b)/a;
    static double coef2 = 1.0/log(b);

    r = domn->rand->getRand();
    iLevel = ceil(log(1.0-r*coef1)*coef2 - 1);

    //--------------- eddy level: option 2

    //r = domn->rand->getRand();
    //iLevel = Nm3;
    //for(int i=0; i<Nm3; i++)
    //    if(r <= eddyLevelCDF[i]){
    //        iLevel = i;
    //        break;
    //    }

    return;

}

///////////////////////////////////////////////////////////////////////////////
/** Function performs eddy events: parcel swaps.
    @param iLevel \input  level of the tree for the base of the swap.
    @param Qstart \output starting index for the Q-tree.
    @param Rstart \output starting index for the R-tree.
    @param nPswap \output number of parcels swapped.

    Randomly select a node on iLevel.
    Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
    Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
    Then swap the cells.

    For a 6 level tree: 0, 1, 2, 3, 4, 5:
    If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
    Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
       or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **

    We use bitwise shifts for easy powers of 2.
    The swap is done by adding or subtracting a value (shift),
        which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
                                                                                                              Level
                                                                                                            ---------
                                                    *                                                           0
                                                 /     \
                                              /           \
                                           /                 \
                                        /                       \
                                     /                             \
                                  /                                   \
                               /                                         \
                            /                                               \
                           *                                                (*)  01|0000                        1
                          / \                                               / \
                        /     \                                           /     \
                      /         \                                       /         \
                    /             \                                   /             \
                  /                 \                               /                 \
                /                     \                           /                     \
               *                       *                         *                       *                      2
              / \                     / \                       / \                     / \
            /     \                 /     \                   /     \                 /     \
          /         \             /         \               /         \             /         \
         *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
        / \         / \         / \         / \           / \         / \         / \         / \
       /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
      *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
     / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
    00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
                                                      ^^^^^^^^^^^                         ^^^^^^^^^^^

*/
void solver_hips::selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap){

    int iTree = domn->rand->getRandInt((1 << iLevel)-1);          // base node (at iLevel) for the swap

    int zero_q = domn->rand->getRandInt(1);                       // 0q where q is 0 or 1
    int one_r  = 2+domn->rand->getRandInt(1);                     // 1r where r is 0 or 1

    Qstart = (zero_q << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));  // starting index of Q parcels
    Rstart = (one_r  << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));  // starting index of R parcels
    nPswap = 1 << (Nm3-iLevel);                                   // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                      // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                      // inclusive indices are Rstart to Rend-1

    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);  // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                    // python: pLoc[Rstart:Rend]=aa


}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

void solver_hips::calculateSolution() {

    domn->io->writeDataFile("hips_init.dat", 0.0);

    //--------------------

    int QS, RS, nPs;

    //---------- advance to the first eddy

    domn->mimx->advanceOdt(0.0, eTL[0].first);     // advance to the first eddy

    //---------- do all eddies but the last: EE, then advance to next, repeat...

    int i;
    for(i=0; i<eTL.size()-1; i++) {

        *domn->io->ostrm << endl << "eddy #, time, level: "
                         << i << "  " << eTL[i].first
                         << "  " << eTL[i].second;
        selectAndSwapTwoSubtrees(eTL[i].second, QS, RS, nPs);
        domn->mimx->advanceOdt(eTL[i].first, eTL[i+1].first, eTL[i].second);

        if(i % domn->pram->modDump == 0){        /*{{{*/
            stringstream ss1; string s1;
            ss1.clear();  ss1 << setfill('0') << setw(4) << i; ss1 >> s1;
            domn->io->writeDataFile("hips_eddy_"+s1+".dat", eTL[i].first);
        }/*}}}*/

    }

    //---------- do the last eddy

    i=eTL.size()-1;
    *domn->io->ostrm << endl << "eddy #, time, level: "
                     << i << "  " << eTL[i].first
                     << "  " << eTL[i].second;
    selectAndSwapTwoSubtrees(eTL[i].second, QS, RS, nPs);
    domn->mimx->advanceOdt(eTL[i].first, domn->pram->tEnd, eTL[i].second);

    *domn->io->ostrm << endl;

    //--------------------

    domn->io->writeDataFile("hips_final.dat", 0.0);
}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 */

//void solver_hips::calculateSolution() {
//
//    domn->io->writeDataFile("hips_init.dat", 0.0);
//
//    //--------------------
//
//    int    iEE = 0;                // eddy counter
//    int    QS, RS, nPs;
//    double time;                   // current simulation time
//    double dt;                     // time to next eddy event (EE)
//    double iLevel;                 // current level of EE
//    double iLevel_p;               // previous level of EE
//
//    time = 0.0;                          // init time
//    iLevel_p  = -1;                      // init init prev level
//    sample_hips_eddy(dt, iLevel);        // init next EE
//
//    while(time+dt <= domn->pram->tEnd) {
//        domn->mimx->advanceOdt(time, time+dt, iLevel_p);       //--- ADVANCE to sampled EE ---
//        selectAndSwapTwoSubtrees(iLevel, QS, RS, nPs);         //--- IMPLEMENT  sampled EE ---
//
//        *domn->io->ostrm << endl << "eddy #, time, level: " << ++iEE << "  " << time+dt << "  " << iLevel;/*{{{*/
//        if(iEE % domn->pram->modDump == 0){
//            stringstream ss1; string s1;
//            ss1.clear();  ss1 << setfill('0') << setw(4) << iEE; ss1 >> s1;
//            domn->io->writeDataFile("hips_eddy_"+s1+".dat", time+dt);
//        }/*}}}*/
//
//        time += dt;                      // init time
//        iLevel_p = iLevel;               // init prev level
//        sample_hips_eddy(dt, iLevel);    // init next EE
//
//
//    }
//    domn->mimx->advanceOdt(time, domn->pram->tEnd, iLevel_p);    //--- ADVANCE ---
//
//    //--------------------
//
//    *domn->io->ostrm << endl;
//    domn->io->writeDataFile("hips_final.dat", 0.0);
//}

